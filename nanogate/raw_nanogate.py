from nanogate.read import Read
from nanogate.read import Part
from nanogate.io import write_raw_json
from nanogate.io import write_raw_jasonlog
from nanogate.utils import *
import time


def raw_nanogate(input_reads,
                 input_parts_library,
                 output_json_reads,
                 output_json_log,
                 construct_parts=10,
                 kmer_size=11):
    """
    :param input_reads: fasta file including reads
    :param input_parts_library: fasta file including the library of all parts used to assembly constructs
    :param input_verification_table: .csv table matching each construct wirh a list of parts
    :param kmer_size: size of each kmer to consider
    :param scale: value used to filter kmer after conversion into hashes
    :param output_txt: .txt file including the accuracy score
    :return: a score reads assigned correctly/number of jaccard scores evaluated


    Takes as input a fasta file -input_reads containing reads, split each read in k-mer of length -kmer_size,
    each k-mer is converted into hash, the hashes can be filtered applying the scale value.
    The library of part is taken as second input -input_parts_library, each part is converted into k-mer of -kmer_size
    kmers are converted into hashes, hashes can be filtered as well as for reads.
    For ecach read jaccard containment coefficient is evaluated between the filtered hashes of the read and the filtered
    ashes of each part
    To verify that each read matches with the right part, a .csv table is uploaded -input_verification_table
    the accuracy score is returned reads assigned correctly/number of jaccard scores evaluated
    A .txt file -output_txt is returned including the accuracy score
    """
    start_time = time.perf_counter()
    #kmer_size = find_kmer_length(input_parts_library)

    # initialising parts
    parts = []
    for record in screed.open(input_parts_library):
        part = Part()
        part.initialise_part(record, kmer_size)
        parts.append(part)

    #building a parts_length dictionary
    part_length = build_parts_length_dictionary(parts)

    # Theorical reads length given by combining parts
    poisson_constructs = poisson_distribution(parts, construct_parts)

    #Initialising reads
    reads = []
    filtered_out = 0
    for record in screed.open(input_reads):
        read = Read()
        read.initialise(record, kmer_size)
        outliars = evalute_likelihood(read, poisson_constructs, reads)
        filtered_out += outliars

    

    # Evaluating jaccard
    for read in reads:
        read.get_json()
        for part in parts:
            jaccard = part.hashes.contained_by(read.hashes)
            read.add_jaccard_to_json(int(part.name),jaccard)
            #read.add_length_to_json(int(part.name),int(part.length))

    end_time =time.perf_counter()
    computational_time = end_time - start_time

    # Write unfiltered result
    write_raw_json(output_json_reads, reads)
    write_raw_jasonlog(computational_time,part_length, filtered_out, output_json_log)


    # # Filtering parts
    # for read in reads:
    #     read.parts_jaccard.sort(key=lambda x: x.jaccard, reverse=True)
    #     sum_parts = 0
    #     parts_list_index = 0
    #     while sum_parts < read.length:
    #         if read.parts_jaccard[parts_list_index].jaccard > threshold:
    #             read.parts.append(read.parts_jaccard[parts_list_index])
    #         sum_parts += read.parts_jaccard[parts_list_index].length
    #         parts_list_index += 1
    #
    #     write_output(output_filename, reads)
