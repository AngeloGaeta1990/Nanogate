Nanogate pseudocode
=========================

**proedure nanogate**
kmer-size <- find_kmer_size(parts_library)

reads_hashes <- convert_read_to_hashed(reads, kmer_size)

parts_hashes <- convert_parts_to_hashes(reads, kmer_size)

jaccards <- evaluate_jaccards(parts_hashes, read_hashes)

filtered_jaccard <- filter(jaccards) # for each read take only parts hashes with the highes jaccard, until sum(len(parts sequence))>= len(read sequence)




