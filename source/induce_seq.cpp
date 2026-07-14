#include "induce_seq.h"
#include "art_framework.h"
#include "fastafile_handler.h"
#include "fileio.h"
#include "random_generator.h"

#include <fstream>
#include <iostream>
#include <map>
#include <algorithm>
#include <cmath>
#include <omp.h>
#include <sys/mman.h>



InduceSeq::InduceSeq(NGSsdd& sddData) : sdd_data(sddData) {
}

InduceSeq::InduceSeq(NGSsdd& sddData, NGSParameters parameters, std::string tempFolderPath) : sdd_data(sddData) {
    parameter = parameters;
    set_genome_data(tempFolderPath);
    set_fragment_size_distribution_from_file();
}

// Reads the DNA fragment size distribution file (same "length count" format as fragment_size_distribution_path)
// from the path set by the induce_seq_fragment_size_distribution_path parameter, and converts the returned
// (min_fragment_length, normalized_counts) pair into a cumulative-probability -> length lookup map, replacing
// the hard-coded default fragment_size_distribution. get_random_fragment_length samples from this map.
void InduceSeq::set_fragment_size_distribution_from_file() {
    auto fragmentData = readFragmentSizeDist(parameter.get_induce_seq_fragment_size_distribution_path());
    int min_fragment_length = fragmentData.first;
    const std::vector<double>& normalized_counts = fragmentData.second;

    fragment_size_distribution.clear();
    double cumulative_probability = 0.0;
    for (size_t i = 0; i < normalized_counts.size(); i++) {
        cumulative_probability += normalized_counts[i];                            // Running sum of the normalized counts turns the PMF into a CDF
        fragment_size_distribution[static_cast<float>(cumulative_probability)] = min_fragment_length + static_cast<int>(i);
    }
}


std::vector<std::vector<long>>& InduceSeq::get_dsb_locations(int groupTID) {
    return(dsb_locations[groupTID]);
}

void InduceSeq::set_parameter(NGSParameters param) {
    parameter = param;
}

// Sets the shape for data vectors, making each vector contain nGroupThreads empty vectors. 
// This is because the data held in InduceSeq class is indexed by thread group ID, in the same way as NGSsdd. 
// Each thread group processes one damaged genome that is to be sequenced, and accesses only the data corresponding to it. 
void InduceSeq::init_set_data_holders(int nGroupThreads){
    dsb_locations.clear();
    dsb_locations.resize(nGroupThreads);

    dsb_blunted_ends.clear();
    dsb_blunted_ends.resize(nGroupThreads);

    dsb_fragments_left.clear();
    dsb_fragments_left.resize(nGroupThreads);

    dsb_fragments_right.clear();
    dsb_fragments_right.resize(nGroupThreads);

    dsb_strands_left.clear();
    dsb_strands_left.resize(nGroupThreads);

    dsb_strands_right.clear();
    dsb_strands_right.resize(nGroupThreads);

    base_pair_damages_left.clear();
    base_pair_damages_left.resize(nGroupThreads);

    base_pair_damages_right.clear();
    base_pair_damages_right.resize(nGroupThreads);
}

//clears all data for a given groupTID. 
void InduceSeq::reset_permanent_damage_vecs(int groupTID){                                                 // empty all the permanent vectors before the next exposure
    dsb_locations[groupTID].clear();
    dsb_blunted_ends[groupTID].clear();
    dsb_fragments_left[groupTID].clear();
    dsb_fragments_right[groupTID].clear();
    dsb_strands_left[groupTID].clear();
    dsb_strands_right[groupTID].clear();
    base_pair_damages_left[groupTID].clear();
    base_pair_damages_right[groupTID].clear();
}

void InduceSeq::run_simulation(int cell_number, int groupTID, int threadID, int NumWorkerThreads, int threadIDOffset) {
    find_DSBs(parameter.get_dsb_threshold(), groupTID);
    get_blunted_ends(groupTID);
    get_dsb_fragments(groupTID, threadID);
    filter_fragments_size(groupTID);
    filter_dsb_strands_ssd(groupTID);
    find_base_pair_damages(groupTID);
    generate_simulation_output(cell_number, groupTID, NumWorkerThreads, threadIDOffset);
}

// Creates a fasta file containing the genome data in a particular format, as described in buildUndamagedGenomeTemplate_ForwardOnly_MM
// Saves this file in tempFolderPath. Sets genome_fasta to point to the first character of the memory map of the fasta file. 
// Initializes cum_chrom_header_sizes and chrom_headers (described in the header file) 
void InduceSeq::set_genome_data(std::string& tempFolderPath) {
    long ref_genomeFile_size = fileSize_bytes(*parameter.get_reference_genome());                              
    std::string genomeTemplatePath = tempFolderPath+"/genome_spaceless.fa";                                       
    genome_fasta_size = static_cast<size_t>(ref_genomeFile_size*2);                                           // The size of an Undamaged file is estimated to be 2 times the size of the reference sequence file
    genome_fasta = createMemoryMappedFile(genomeTemplatePath, genome_fasta_size);                                     // Generate a memory-map placeholder to store the memory map of the undamaged fasta file as it gets created later
    buildUndamagedGenomeTemplate_ForwardOnly_MM(genome_fasta, genome_fasta_size, sdd_data.get_num_chrom(), sdd_data.get_chrom_mapping(), parameter.get_reference_genome(), cum_chrom_header_sizes);
    if (calculateCumChromHeaderSizes(cum_chrom_header_sizes, chrom_headers, genome_fasta, genome_fasta_size, *sdd_data.get_chrom_end_loc())) {
        std::cerr<<"\n ERROR: The chromosome sizes listed in the sdd file do not match the chromosome sizes in the genome fasta file " << parameter.get_reference_genome();
        exit(EXIT_FAILURE);
    }
}

// finds DSBs from the ssd data, for a given groupTID. Populates the dsb_locations[groupTID] vector (format described in header file). 
// DSBthreshold is the maximum distance between a pair of strand breaks for them to qualify as a dsb. 
// ssd_data should have its backbone break vectors initialized and sorted least position to greatest position before this function is called
void InduceSeq::find_DSBs(int DSBthreshold, int groupTID){
    const std::vector<long>& chromEnds = *sdd_data.get_chrom_end_loc();                                         // chrom_end_loc is sorted in ascending order

    // Each chromosome is a separate, already-blunted DNA molecule, so both backbones "break" at its two physical
    // ends. chromEnds itself already lists exactly these boundary positions (num_chrom+1 of them): position 0 is
    // the start of the first chromosome, chromEnds.back() is the end of the last chromosome, and each entry in
    // between is simultaneously the end of one chromosome and the start of the next, since global coordinates run
    // contiguously across the chromosome boundary. These are merged into local copies of the backbone break lists
    // (both real ssd breaks and chromEnds are already sorted ascending, so a linear merge keeps the result sorted)
    // so that the dsb-finding logic below, completely unmodified, can also find dsbs formed between a chromosome
    // boundary and a nearby real break. The original break vectors held by sdd_data are left untouched, since
    // later steps (e.g. filter_dsb_strands_ssd) read them again expecting only genuine ssd breaks.
    std::vector<long>& real_backbone1_breaks = sdd_data.get_backbone1_break_loc(groupTID);
    std::vector<long>& real_backbone2_breaks = sdd_data.get_backbone2_break_loc(groupTID);

    std::vector<long> backbone1_breaks;
    backbone1_breaks.reserve(real_backbone1_breaks.size() + chromEnds.size());
    std::merge(real_backbone1_breaks.begin(), real_backbone1_breaks.end(), chromEnds.begin(), chromEnds.end(), std::back_inserter(backbone1_breaks));

    std::vector<long> backbone2_breaks;
    backbone2_breaks.reserve(real_backbone2_breaks.size() + chromEnds.size());
    std::merge(real_backbone2_breaks.begin(), real_backbone2_breaks.end(), chromEnds.begin(), chromEnds.end(), std::back_inserter(backbone2_breaks));

    std::vector<long>::iterator site1 = backbone1_breaks.begin();
	std::vector<long>::iterator site2 = backbone2_breaks.begin();

    // variable to keep track of whether the previous step in the while loop was a dsb
    // used for keeping track of dsbs that are part of interconnected dsb, as explained below
    bool prevStepDSB = 0;
    
    //go through backbone breaks, checking for dsbs.  
    while (site1 != backbone1_breaks.end() && site2 != backbone2_breaks.end()){
        int siteDiff = *site1 - *site2;                                                                 // separation in number of bp
        // #pragma omp critical 
        // {
        //     std::cout << "site 1: " << *site1 << ",  site 2: " << *site2 << "\n";
        // }
		bool isDSB{0};                                                                                  // initiating with zero
        if(abs(siteDiff) <= DSBthreshold){
            // Chromosome identity itself is not kept: it is only needed transiently here, to check that the two
            // sites fall in the same interval (get_chrom_idx is the same binary search used wherever a chromosome
            // index is actually needed downstream, e.g. generate_simulation_output).
            long chromIdx1 = get_chrom_idx(*site1);
            long chromIdx2 = get_chrom_idx(*site2);
            isDSB = (chromIdx1 == chromIdx2);                                                          // same chromosome if both sites fall in the same interval
        }

        if(isDSB){
            // std::cout << "loc: " << *site1 << "\n";
            // save dsb to data vector
            dsb_locations[groupTID].push_back({*site1, *site2, prevStepDSB});


            // One site is incremented, to move onto the next strand break. The choice of which site to increment has the form below for the following reason. 
            // Consider the strand breaks as nodes in a graph, with an edge between any two breaks that satifsy the distance criterion for a dsb. 
            // If dsbs are all isolated, then the graph is a set of 2-node connected components, and no special consideration is needed.
            // However, if many breaks are close together, then the graph can include a connected component with many nodes. 
            // In such a component, we are only interested in the first and the last pair of connected nodes, in terms of position in the genome:
            // these two pairs will form the edges of two dna segments, and all the DNA in between them will get all broken up, and will be irrelevent to induce-seq
            // The simplest way to find and record these pairs is with the incrementing method below.
            // This method is gauranteed to detect the first and last connected pair in any connected component, and on every iteration step in between, it will detect a dsb
            // The fourth data field in recorded dsbs is whether the previous step was a dsb. 
            // Thus, when parsing through dsbs, if consecutive dsbs have 1 in that field, they are part of a set of interconnected dsbs.
            // In that case, the last dsb in that sequence and dsb directly before that, which has a 0 in the fourth field, are the dsbs to consider for InduceSeq
            bool site1_has_next = (site1 + 1) != backbone1_breaks.end();
            bool site2_has_next = (site2 + 1) != backbone2_breaks.end();
            if (site1_has_next && site2_has_next) {
                if ( (*(site1+1) - *site2) < (*(site2+1) - *site1)) {
                    site1++;
                } else {
                    site2++;
                }
            } else if (site1_has_next) {
                // site2 has no more breaks to compare against; advance site1 so its next value can still be checked against the current site2
                site1++;
            } else {
                // Either site2 has a next break to check against the current site1, or neither side has one left (in which case
                // this just advances site1 to end the loop; site2's value has already been fully considered by this point)
                site2++;
            }
        }else{
            // move to the next strand break
            if (siteDiff>0){site2++;}                                                                   // if site1 is after site2; update site2
            else if (siteDiff<0){site1++;}                                                              // if site1 is before site2; update site1
        }
        prevStepDSB = isDSB;
    }

    // const auto& dsbs = dsb_locations[groupTID];
    // std::cout << "DSB locations for groupTID=" << groupTID << " (" << dsbs.size() << " DSBs):\n";
    // for (size_t i = 0; i < dsbs.size(); i++) {
    //     std::cout << "  [" << i << "] backbone1=" << dsbs[i][0]
    //               << "  backbone2=" << dsbs[i][1]
    //               << "  prevDSB="  << dsbs[i][2] << "\n";
    // }
}

// Returns the 0-based index of the chromosome that a global bp position falls within, via binary search over
// chrom_end_loc. Chromosome c owns positions chrom_end_loc[c]+1 .. chrom_end_loc[c+1] inclusive (matching the
// convention used everywhere else, e.g. chrom_start/chrom_end in the fragment-generation code), so this looks
// for the smallest chrom_end_loc entry that is >= position rather than the largest one that is <= position:
// a position exactly on a chromosome boundary is that chromosome's own last base, not the next chromosome's
// first. This is the only place chromosome identity is derived from a position outside of find_DSBs.
int InduceSeq::get_chrom_idx(long position) {
    const std::vector<long>& chromEnds = *sdd_data.get_chrom_end_loc();
    return std::lower_bound(chromEnds.begin(), chromEnds.end(), position) - chromEnds.begin() - 1;
}

void InduceSeq::close() {
    munmap(genome_fasta, genome_fasta_size);
}

void InduceSeq::get_blunted_ends(int groupTID) {
    dsb_blunted_ends[groupTID].clear();
    std::vector<std::vector<long>> dsb_locs = get_dsb_locations(groupTID);
    // dsb_blunted ends is in the form {location of base on the left edge, location of base on the right edge,
    // strand1 (backbone1) location of the dsb that caused the left edge, strand2 (backbone2) location of the dsb that caused the left edge,
    // strand1 (backbone1) location of the dsb that caused the right edge, strand2 (backbone2) location of the dsb that caused the right edge}
    // A cluster of chained/connected dsbs can have the left edge caused by a different dsb than the right edge, so both are tracked separately.
    // Chromosome identity is deliberately not tracked here: find_DSBs treats every chromosome boundary as a break on both strands, so a
    // boundary always shows up as an ordinary entry in dsb_locs, indistinguishable from a genuine dsb as far as this function is concerned
    // (see get_dsb_fragments for how that is enough, on its own, to keep fragments from crossing a chromosome boundary).
    // base locations start at 1
    std::vector<long> new_dsb_blunted_ends = {0, 0, 0, 0, 0, 0};
    for (size_t i = 0; i < dsb_locs.size(); i++) {
        std::vector<long> dsb = dsb_locs[i];
        if (dsb[2] == 0) {
            new_dsb_blunted_ends[0] = dsb[1];
            new_dsb_blunted_ends[2] = dsb[0];                                     // strand1 location of the dsb that starts this cluster (causes the left edge)
            new_dsb_blunted_ends[3] = dsb[1];                                     // strand2 location of the dsb that starts this cluster (causes the left edge)
        }
        if (i + 1 == dsb_locs.size() || dsb_locs[i+1][2] == 0) {
            new_dsb_blunted_ends[1] = dsb_locs[i][0] + 1;
            new_dsb_blunted_ends[4] = dsb_locs[i][0];                             // strand1 location of the dsb that closes this cluster (causes the right edge)
            new_dsb_blunted_ends[5] = dsb_locs[i][1];                             // strand2 location of the dsb that closes this cluster (causes the right edge)
            dsb_blunted_ends[groupTID].push_back(new_dsb_blunted_ends);
        }
    }

    // const auto& blunted = dsb_blunted_ends[groupTID];
    // std::cout << "Blunted ends for groupTID=" << groupTID << " (" << blunted.size() << " entries):\n";
    // for (size_t i = 0; i < blunted.size(); i++) {
    //     std::cout << "  [" << i << "] left=" << blunted[i][0]
    //               << "  right=" << blunted[i][1] << "\n";
    // }
}

//DDDD
void InduceSeq::get_dsb_fragments(int groupTID, int threadID) {
    dsb_fragments_left[groupTID].clear();
    dsb_fragments_right[groupTID].clear();

    // Chromosome identity is never tracked here. find_DSBs treats every chromosome boundary as a break on both
    // strands, so a boundary always shows up as an ordinary entry in dsb_blunted_ends, exactly like a real dsb -
    // the same overlap handling below that keeps two neighbouring fragments from encroaching on each other
    // therefore also keeps a fragment from crossing into the next chromosome, with no need to know where
    // chromosome boundaries fall. The only thing that is special-cased is the two ends of the whole genome: the
    // very first entry always originates (at least in part) from the break at the start of the genome, and the
    // very last entry always originates from the break at its end (both breaks are unconditionally present, and
    // are respectively the smallest and largest possible break position), so those two never grow a fragment on
    // the side that would run off the genome.
    const std::vector<std::vector<long>>& blunted_ends = dsb_blunted_ends[groupTID];

    // State carried over from the previous entry's right fragment, so the current entry's left fragment
    // can be checked against it for overlap.
    long previous_right_start = 0;
    long previous_right_end = 0;
    long previous_right_dsb_strand1 = 0;                                            // The right-causing dsb's strand1/strand2 locations for the previous blunted end,
    long previous_right_dsb_strand2 = 0;                                            // carried over in case the previous dsb's right fragment needs to be re-pushed after an overlap merge
    bool previous_right_appended = false;                                          // Whether the previous entry's right fragment is currently the last element of dsb_fragments_right[groupTID]

    // Set by the right-fragment check below when entry i's right fragment runs into entry i+1's left
    // blunted end, so that entry i+1's left fragment is skipped entirely on the next loop iteration.
    bool skip_left_fragment = false;

    size_t i = 0;
    while (i < blunted_ends.size()) {
        const std::vector<long>& dsb_blunted_end = blunted_ends[i];
        long left_dsb_strand1 = dsb_blunted_end[2];                                // strand1 location of the dsb that caused the left edge
        long left_dsb_strand2 = dsb_blunted_end[3];                                // strand2 location of the dsb that caused the left edge
        long right_dsb_strand1 = dsb_blunted_end[4];                               // strand1 location of the dsb that caused the right edge
        long right_dsb_strand2 = dsb_blunted_end[5];                               // strand2 location of the dsb that caused the right edge

        bool skip_this_left_fragment = skip_left_fragment || i == 0;               // Capture the flag set by the previous iteration; the very start of the genome never grows a left fragment
        skip_left_fragment = false;

        int fragment_length = get_random_fragment_length(threadID) - parameter.get_P5_adapter_length();
        if (fragment_length > 1 && !skip_this_left_fragment) {
            long left_start = dsb_blunted_end[0];
            long left_end = left_start - fragment_length + 1;
            bool drop_current_left = false;
            if (left_end < previous_right_end) {
                if (previous_right_appended) {
                    // The current left fragment overlaps with the previous entry's right fragment.
                    long overlap = previous_right_end - left_end;                      // How far the two fragments overlap
                    long blunted_end_gap = left_start - previous_right_start;          // Raw distance between the two dsb break points (not the fragment ends)
                    if (overlap > parameter.get_maximum_overlap_fragment_generation() * (blunted_end_gap + 2*parameter.get_P5_adapter_length())) {
                        // Overlap is large; the section in between these two DSBs is taken to not have been fragmented. No reads produced from this fragment, since it contains 2 P5 adapters.
                        // right fragment that was already pushed, and skip pushing the current left fragment.
                        dsb_fragments_right[groupTID].pop_back();
                        drop_current_left = true;
                    } else {
                        // Overlap is small enough: both fragments are shortened to meet near the midpoint of
                        // their two ends, but with two consecutive end positions rather than sharing a base pair.
                        long end_sum = left_end + previous_right_end;
                        long new_right_end;                                        // previous dsb's right fragment new end (the smaller of the two)
                        long new_left_end;                                         // current dsb's left fragment new end (the larger of the two)
                        if (end_sum % 2 != 0) {
                            // Mean is not an integer: use the two integers either side of it.
                            new_right_end = end_sum / 2;
                            new_left_end = new_right_end + 1;
                        } else {
                            // Mean is an integer: one fragment gets that exact position, the other is 1bp away.
                            // Which fragment gets the exact position alternates with the parity of the mean, so
                            // neither the left nor the right fragments are systematically longer.
                            long mean = end_sum / 2;
                            if (mean % 2 == 0) {
                                new_right_end = mean;
                                new_left_end = mean + 1;
                            } else {
                                new_left_end = mean;
                                new_right_end = mean - 1;
                            }
                        }
                        left_end = new_left_end;

                        // Re-push the previous right fragment with its new (shortened) end. Size filtering happens
                        // later, in filter_fragments_size, once all fragments have been generated.
                        dsb_fragments_right[groupTID].pop_back();
                        dsb_fragments_right[groupTID].push_back({previous_right_start, new_right_end, previous_right_dsb_strand1, previous_right_dsb_strand2});
                    }
                } else {
                    // No previous fragment was actually kept (e.g. it was too short to be a fragment at all):
                    // its territory still shouldn't be encroached on.
                    left_end = previous_right_end;
                }
            }

            if (!drop_current_left) {
                dsb_fragments_left[groupTID].push_back({left_start, left_end, left_dsb_strand1, left_dsb_strand2});
            }
        }

        fragment_length = get_random_fragment_length(threadID) - parameter.get_P5_adapter_length();
        long right_start = dsb_blunted_end[1];
        long right_end = right_start + fragment_length - 1;
        bool right_appended = false;
        bool is_last = (i + 1 == blunted_ends.size());
        if (fragment_length > 1 && !is_last) {
            // Check whether this right fragment runs into the next entry's left blunted end (plus the adapter length).
            // If so, drop this right fragment and remember to also skip the next entry's left fragment.
            const std::vector<long>& next_dsb_blunted_end = blunted_ends[i+1];
            bool overlaps_next_dsb = right_end > next_dsb_blunted_end[0] + parameter.get_P5_adapter_length();
            if (overlaps_next_dsb) {
                skip_left_fragment = true;
            } else {
                dsb_fragments_right[groupTID].push_back({right_start, right_end, right_dsb_strand1, right_dsb_strand2});
                right_appended = true;
            }
        }

        previous_right_start = right_start;
        previous_right_end = right_end;
        previous_right_dsb_strand1 = right_dsb_strand1;
        previous_right_dsb_strand2 = right_dsb_strand2;
        previous_right_appended = right_appended;

        i++;
    }

    // const auto& frags_left  = dsb_fragments_left[groupTID];
    // const auto& frags_right = dsb_fragments_right[groupTID];
    // std::cout << "DSB fragments left for groupTID=" << groupTID << " (" << frags_left.size() << " entries):\n";
    // for (size_t j = 0; j < frags_left.size(); j++) {
    //     std::cout << "  [" << j << "] start=" << frags_left[j][0]
    //               << "  end="   << frags_left[j][1] << "\n";
    // }
    // std::cout << "DSB fragments right for groupTID=" << groupTID << " (" << frags_right.size() << " entries):\n";
    // for (size_t j = 0; j < frags_right.size(); j++) {
    //     std::cout << "  [" << j << "] start=" << frags_right[j][0]
    //               << "  end="   << frags_right[j][1] << "\n";
    // }
}

// Removes fragments from dsb_fragments_left[groupTID] and dsb_fragments_right[groupTID] that are shorter than
// parameter.get_first_size_filter() once the P5 adapter is accounted for. get_dsb_fragments appends every
// fragment it generates unfiltered, so that overlap handling between neighbouring fragments always sees the
// true previous fragment rather than treating one dropped for being too short as if it were never generated.
void InduceSeq::filter_fragments_size(int groupTID) {
    std::vector<std::vector<long>>& frags_left = dsb_fragments_left[groupTID];
    frags_left.erase(std::remove_if(frags_left.begin(), frags_left.end(), [this](const std::vector<long>& frag) {
        return frag[0] - frag[1] + parameter.get_P5_adapter_length() + 1 < parameter.get_first_size_filter();
    }), frags_left.end());

    std::vector<std::vector<long>>& frags_right = dsb_fragments_right[groupTID];
    frags_right.erase(std::remove_if(frags_right.begin(), frags_right.end(), [this](const std::vector<long>& frag) {
        return frag[1] - frag[0] + parameter.get_P5_adapter_length() + 1 < parameter.get_first_size_filter();
    }), frags_right.end());
}



// void InduceSeq::get_dsb_fragments(int groupTID, int threadID) {
//     dsb_fragments_left[groupTID].clear();
//     dsb_fragments_right[groupTID].clear();
    
//     long previous_right_end = 0;
//     bool previous_right_end_in_next_adaptor = false;
//     size_t i = 0;
//     while (i < dsb_blunted_ends[groupTID].size()) {
//         std::vector<long> dsb_blunted_end = dsb_blunted_ends[groupTID][i]; 
//         int chrom_idx = dsb_blunted_end[2];  
        

//         int fragment_length = get_random_fragment_length(threadID) - parameter.get_P5_adapter_length();
//         if (fragment_length > 1) {
//             long left_start = dsb_blunted_end[0];
//             long left_end = left_start - fragment_length + 1;
//             int chrom_start = (*sdd_data.get_chrom_end_loc())[chrom_idx] + 1;
//             if (left_end < chrom_start) {
//                 left_end = chrom_start;
//             } else if (left_end < previous_right_end) {
//                 left_end = 0.5 * (left_end + previous_right_end);
//                 std::vector<long> prev_dsb_frag = dsb_fragments_right[groupTID].back();
//                 dsb_fragments_right.pop_back();
//                 if (left_end - prev_dsb_frag[0] + parameter.get_P5_adapter_length() + 1 >= parameter.get_first_size_filter()) {
//                     dsb_fragments_right[groupTID].push_back({prev_dsb_frag[0], left_end, chrom_idx});
//                 }
//             }
            
//             if (left_start - left_end + parameter.get_P5_adapter_length() + 1 >= parameter.get_first_size_filter()) {
//                 dsb_fragments_left[groupTID].push_back({left_start, left_end, chrom_idx});
//             }
//         }

//         fragment_length = get_random_fragment_length(threadID) - parameter.get_P5_adapter_length();
//         long right_start = dsb_blunted_end[1];
//         long right_end = right_start + fragment_length - 1;
//         if (fragment_length > 1) {    
//             int chrom_end = (*sdd_data.get_chrom_end_loc())[chrom_idx + 1];
//             // std::cout << "chrom end: " <<chrom_end << "  right_end: " << right_end << "\n";
//             if (right_end > chrom_end) right_end = chrom_end;
//             if (right_end - right_start + parameter.get_P5_adapter_length() + 1 >= parameter.get_first_size_filter()) {
//                 dsb_fragments_right[groupTID].push_back({right_start, right_end, chrom_idx});
//             }
//         }
//         previous_right_end = right_end;
        
//         i++;
//     }

// }



void InduceSeq::filter_dsb_strands_ssd(int groupTID) {
    dsb_strands_left[groupTID].clear();
    dsb_strands_right[groupTID].clear();
    std::vector<long> strand_1_breaks = sdd_data.get_backbone1_break_loc(groupTID);
    std::vector<long> strand_2_breaks = sdd_data.get_backbone2_break_loc(groupTID);

    size_t i = 0;
    size_t ssd_i = 0;
    while (i < dsb_fragments_left[groupTID].size()) {
        std::vector<long> dsb_frag = dsb_fragments_left[groupTID][i];
        long frag_end = dsb_frag[1];
        long causing_break = dsb_frag[2];                                          // backbone1 break that formed this dsb

        bool is_good = true;

        // checks for ssbs 
        // Because of overhang fill-in during blunting, the edge of a fragment can be at a different location than the original dsb edge. 
        // damages with positions in between these two positions would not be on the fragment, since that part of the fragment is remade. 
        // That is why causing_break is used as a bound, instead of the end of the break.  
        
        while (ssd_i < strand_1_breaks.size() && strand_1_breaks[ssd_i] < causing_break) {
            if (strand_1_breaks[ssd_i] >= frag_end) {
                is_good = false;
            }
            ssd_i++;
        }

        if (is_good) {
            dsb_strands_left[groupTID].push_back(dsb_frag);
        }

        
        i++;
    }

    i = 0;
    ssd_i = 0;
    while (i < dsb_fragments_right[groupTID].size()) {
        std::vector<long> dsb_frag = dsb_fragments_right[groupTID][i];
        long frag_start = dsb_frag[0];
        long frag_end = dsb_frag[1];
        long causing_break = dsb_frag[3];                                          // strand2 location of the dsb that caused the right edge

        bool is_good = true;

        // Because of overhang fill-in during blunting, the edge of a fragment can be at a different location than the original dsb edge. 
        // damages with positions in between these two positions would not be on the fragment, since that part of the fragment is remade. This loop filters them out. 
        while (ssd_i < strand_2_breaks.size() && strand_2_breaks[ssd_i] <= causing_break) {
            ssd_i++;
        }
        // Breaks strictly after the causing break are genuinely 
        while (ssd_i < strand_2_breaks.size() && strand_2_breaks[ssd_i] < frag_end) {
            if (strand_2_breaks[ssd_i] >= frag_start) {
               is_good = false;
            }
            ssd_i++;
        }

        if (is_good) {
            dsb_strands_right[groupTID].push_back(dsb_frag);
        }
        
        i++;

    }

    const auto& strands_left  = dsb_strands_left[groupTID];
    const auto& strands_right = dsb_strands_right[groupTID];
    std::cerr << "DSB strands left for groupTID=" << groupTID << " (" << strands_left.size() << " entries):\n";
    for (size_t j = 0; j < strands_left.size(); j++) {
        std::cerr << "  [" << j << "] start=" << strands_left[j][0]
                  << "  end="   << strands_left[j][1] << "\n";
    }
    std::cerr << "DSB strands right for groupTID=" << groupTID << " (" << strands_right.size() << " entries):\n";
    for (size_t j = 0; j < strands_right.size(); j++) {
        std::cerr << "  [" << j << "] start=" << strands_right[j][0]
                  << "  end="   << strands_right[j][1] << "\n";
    }
}

void InduceSeq::find_base_pair_damages(int groupTID) {
    base_pair_damages_left[groupTID].clear();
    base_pair_damages_right[groupTID].clear();

    size_t bp_i = 0;
    std::vector<long> bp_damages = sdd_data.get_basestrand1_damage_loc(groupTID);
    for (std::vector<long> dsb_strand : dsb_strands_left[groupTID]) {
        std::vector<long> bp_damages_in_strand;
        while (bp_i < bp_damages.size() && bp_damages[bp_i] <= dsb_strand[0]) {
            if (bp_damages[bp_i] >= dsb_strand[1]) {
                bp_damages_in_strand.push_back(bp_damages[bp_i]);
            }
            bp_i++;
        }
        base_pair_damages_left[groupTID].push_back(bp_damages_in_strand);
    }

    bp_i = 0;
    bp_damages = sdd_data.get_basestrand2_damage_loc(groupTID);
    for (std::vector<long> dsb_strand : dsb_strands_right[groupTID]) {
        std::vector<long> bp_damages_in_strand;
        while (bp_i < bp_damages.size() && bp_damages[bp_i] <= dsb_strand[1]) {
            if (bp_damages[bp_i] >= dsb_strand[0]) {
                bp_damages_in_strand.push_back(bp_damages[bp_i]);
            }
            bp_i++;
        }
        base_pair_damages_right[groupTID].push_back(bp_damages_in_strand);
    }
}


void InduceSeq::generate_simulation_output(int cell_number, int groupTID, int num_available_threads, int threadIDOffset) {
    
    ART::initiate_read_generation(parameter.get_read_length(), parameter.get_GC_binSize(), parameter.get_fraction_nonFR_read_pairs(), parameter.get_read_artifacts_rate());
    ART::set_read_quality_distribution(*parameter.get_r1_quality_profile(), *parameter.get_r2_quality_profile());

    // ofstream object of the output fastq file for read 1
    std::string output_file_ending;
    std::string output_fastq_R1_filename;
    std::ofstream fastq_R1_file;
    if (parameter.get_compress_output()) {
        output_file_ending = ".fastq.gz";
        output_fastq_R1_filename = (*parameter.get_output_directory())+"/"+(*parameter.get_output_fastq_filename_prefix())+"_"+std::to_string(cell_number)+"_R1" + output_file_ending;    
        fastq_R1_file.open(output_fastq_R1_filename.c_str(),std::ios::binary);
    } else {
        output_file_ending = ".fastq";
        output_fastq_R1_filename = (*parameter.get_output_directory())+"/"+(*parameter.get_output_fastq_filename_prefix())+"_"+std::to_string(cell_number)+"_R1" + output_file_ending;
        fastq_R1_file.open(output_fastq_R1_filename.c_str());
    } 
    
    // ofstream object of the output file listing the DSBs that were sequenced, if requested. Never compressed.
    std::string output_sequenced_dsbs_filename;
    std::ofstream sequenced_dsbs_file;
    if (parameter.get_output_sequenced_dsbs()) {
        output_sequenced_dsbs_filename = (*parameter.get_output_directory())+"/"+(*parameter.get_output_fastq_filename_prefix())+"_"+std::to_string(cell_number)+"_sequenced_dsbs.csv";
        sequenced_dsbs_file.open(output_sequenced_dsbs_filename.c_str());
        sequenced_dsbs_file << "fragment_start,fragment_end,chromosome_index,is_left,strand1_damage_location,strand2_damage_location,read_number\n";
    }
    
    ART read1;                                                                                      // Creating an ART class object and setting the insertion and deletion probability vectors for that read object
    read1.set_read_error_rates(parameter.get_insertion_error_rate_read1(), parameter.get_deletion_error_rate_read1());
    read1.set_read_error_probability(parameter.get_read_length(), parameter.get_insertion_error_rate_read1(), read1.insertion_probability_vec, parameter.get_max_errors_in_read());
    read1.set_read_error_probability(parameter.get_read_length(), parameter.get_deletion_error_rate_read1(), read1.deletion_probability_vec, parameter.get_max_errors_in_read());
    read1.resize_vectors(parameter.get_number_of_threads());                                         // Sized to the full global thread count since threadID below spans that range, not just this call's num_available_threads
    
    std::string chromSegSeq;                                                                        // Temporary variable to hold each chromosome segment sequence from the fasta file one at a time
    std::string chromSegSeq_ID;                                                                     // Temporary variable to hold IDs of each hromosome segment sequence


    const int batchSize{2000};                                                                      // Define a batch size for writing reads to the output file. These much data will be stored in cache before writing it on the file

    int batchSize_thread = std::round(batchSize/num_available_threads);                                     // Devide the total cache size for the buffer equally for all the threads

    std::vector<std::vector<std::string>> batch_buffer(num_available_threads);
    std::vector<std::vector<std::string>> dsb_batch_buffer(num_available_threads);                   // Per-thread buffer for the sequenced-DSBs CSV lines, flushed the same way as batch_buffer

    #pragma omp parallel for num_threads(num_available_threads)
    for (size_t i=0; i<dsb_strands_left[groupTID].size() + dsb_strands_right[groupTID].size(); i++) {
        int localTID = omp_get_thread_num();                                                       // ID local to this call's team (0..num_available_threads-1); safe to index batch_buffer, which is private to this call
        int threadID = threadIDOffset + localTID;                                                  // Globally-unique ID (0..nThreads_User-1) required by rng::local_mt and ART's per-thread buffers, which are shared across all concurrently-running groups
        std::string dna_seq;
        std::vector<long> dsb_strand;
        std::vector<long> bp_damages;
        bool is_left;
        if (i < dsb_strands_left[groupTID].size()) {
            dsb_strand = dsb_strands_left[groupTID][i];
            bp_damages = base_pair_damages_left[groupTID][i];
            is_left = true;
        } else {
            dsb_strand = dsb_strands_right[groupTID][i - dsb_strands_left[groupTID].size()];
            bp_damages = base_pair_damages_right[groupTID][i - dsb_strands_left[groupTID].size()];
            is_left = false;
        }
        // Chromosome identity isn't carried through the pipeline (see get_dsb_fragments); this is the one place
        // it's actually needed, so it's looked up here from dsb_strand[0], the fragment's end nearest to the
        // dsb that caused it, which is guaranteed to sit well within the correct chromosome.
        int chrom_idx = get_chrom_idx(dsb_strand[0]);

        if (parameter.get_output_sequenced_dsbs()) {
            std::string dsb_data = std::to_string(dsb_strand[0])+","+std::to_string(dsb_strand[1])+","+std::to_string(chrom_idx)+","+std::to_string(is_left)+","+std::to_string(dsb_strand[2])+","+std::to_string(dsb_strand[3])+","+std::to_string(i)+"\n";
            dsb_batch_buffer[localTID].push_back(dsb_data);                         // Add the DSB data to the buffer vector of the respective thread

            if (dsb_batch_buffer[localTID].size() >= static_cast<size_t>(batchSize_thread)) {// Check if the batch buffer is full, and write it to the file if needed.
                #pragma omp critical(section2)
                {
                    writeBatchToFile(dsb_batch_buffer[localTID], sequenced_dsbs_file, false);
                }
            }
        }

        get_dna_sequence(dna_seq, bp_damages, dsb_strand, is_left, chrom_idx);
        // std::cout << dna_seq << "\n";
        read1.generate_read_with_indel_from_frag(dna_seq, threadID);                                   // Make a read with random indel errors
        std::vector<short> read1_quality_score_vec;                                 // Vector to hold the quality scores for read 1
        read1.get_read_quality(read1_quality_score_vec, 1, threadID);               // Get the read quality scores for the read positions
        read1.add_baseCall_error(read1_quality_score_vec, threadID);                // Add base call errors to the read based on the quality scores

        //BBBB
        std::string chromID = chrom_headers[chrom_idx];
        std::string read_data = "@"+chromID+"_read"+std::to_string(i)+"\n";                           // @readID
        read_data += (*read1.get_final_read_sequence(threadID))+ "\n+\n";           // read sequence and +
        for(size_t k=0; k<(*read1.get_final_read_sequence(threadID)).size(); k++){  // read quality scores; insert only as many quality values as with the length of sequence
            read_data += static_cast<char>(read1_quality_score_vec[k]+32);          // +33 to get the phred score
        }
        read_data += "\n";
        // std::cout<< read_data;

        batch_buffer[localTID].push_back(read_data);                                // Add the read data to the buffer vector of the respective thread

        if (batch_buffer[localTID].size() >= static_cast<size_t>(batchSize_thread)) {// Check if the batch buffer is full, and write it to the file if needed.
            if (parameter.get_compress_output()) {
                // compression can be done in paralell, since there is no shared memory.
                std::string compressed_batch = getCompressedBatch(batch_buffer[localTID]);
                // writing can only be done by one thread at a time.
                // writeBatchToFile with compression = true is not used so that compression and writing can be done in separate blocks
                #pragma omp critical(section1)
                {
                    fastq_R1_file.write(compressed_batch.c_str(), compressed_batch.size());
                }
            } else {
                #pragma omp critical(section1)
                {
                    writeBatchToFile(batch_buffer[localTID], fastq_R1_file, false);
                }
            }
        }
    }
    for (size_t l=0;l<batch_buffer.size();l++){
        writeBatchToFile(batch_buffer[l], fastq_R1_file, parameter.get_compress_output());      // If there are unwritten data in batch buffer, write that too when the loop ends
    }

    fastq_R1_file.close();
    if (parameter.get_output_sequenced_dsbs()) {
        for (size_t l=0;l<dsb_batch_buffer.size();l++){
            writeBatchToFile(dsb_batch_buffer[l], sequenced_dsbs_file, false);      // If there is unwritten data in the dsb batch buffer, write that too when the loop ends
        }
        sequenced_dsbs_file.close();
    }
}

void InduceSeq::get_dna_sequence(std::string& dna_seq, std::vector<long>& bp_damages, std::vector<long>& dsb_strand, bool is_left, int chrom_idx) {
    long start_char_i = dsb_strand[0] + cum_chrom_header_sizes[chrom_idx] - 1;
    long end_char_i = dsb_strand[1] + cum_chrom_header_sizes[chrom_idx] - 1;

    if (is_left) {
        dna_seq = std::string(genome_fasta + end_char_i, genome_fasta + start_char_i + 1);
        for (long bp_damage : bp_damages) {
            int dam_i = bp_damage - dsb_strand[1];
            dna_seq[dam_i] = 'N';
        }
        std::reverse(dna_seq.begin(), dna_seq.end());
    } else {
        dna_seq = std::string(genome_fasta + start_char_i, genome_fasta + end_char_i + 1);
        for (long bp_damage : bp_damages) {
            int dam_i = bp_damage - dsb_strand[0];
            dna_seq[dam_i] = 'N';
        }
        for (char& base : dna_seq) {
            if      (base == 'A') base = 'T';
            else if (base == 'T') base = 'A';
            else if (base == 'C') base = 'G';
            else if (base == 'G') base = 'C';
        }
    }
}

int InduceSeq::get_random_fragment_length(int threadID) {
    float rand_val = rng::rand_float(0.0f, 1.0f, threadID);
    int fragment_length = fragment_size_distribution.upper_bound(rand_val)->second;
    return fragment_length;
}
