from nanogate.io import read_raw_json
from nanogate.io import write_filtered_output_csv
from nanogate.io import read_part_length
from nanogate.read import Read
from nanogate.read import Part
from nanogate.utils import *


def filtered_nanogate(raw_data: 'nanogate raw data',
                      raw_log: 'raw log file',
                      filtered_output: 'filtered output',
                      threshold="00"):

    """

    :param filtered_output:
    :param threshold:
    :return:
    """
    #Building reads
    reads = []
    reads_dict_list = read_raw_json(raw_data)
    for read_dict in reads_dict_list:
        read = Read()
        read.build_read_from_dictionary(read_dict)
        reads.append(read)


    #Reading threshold jacard
    threshold = read_theshold(threshold)

    #Reading part length
    part_length = read_part_length(raw_log)


    #Building parts
    for read in reads:
         for part_name in list(read.parts_jaccard.keys()):
             part = Part()
             part.build_part_from_json(part_name,part_length,read)
             read.unfiltered_parts.append(part)

    # Filtering parts
    for read in reads:
        sum_parts = 0
        parts_list_index = 0
        read.unfiltered_parts.sort(key=lambda x: x.jaccard, reverse=True)
        while sum_parts < read.length:
            if read.unfiltered_parts[parts_list_index].jaccard > threshold:
                read.parts.append(read.unfiltered_parts[parts_list_index])
            sum_parts += read.unfiltered_parts[parts_list_index].length
            parts_list_index += 1

    write_filtered_output_csv(filtered_output, reads)



