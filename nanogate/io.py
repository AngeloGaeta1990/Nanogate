import pandas as pd
import json


def write_raw_json(output_filename,reads):
    """
    Method to write raw jason file
    mapping constructs to parts jaccard
    :param output_filename:
    :param reads:
    :return:
    """
    jsons = []
    for read in reads:
        jsons.append(read.json)
    with open(output_filename, 'w') as outfile:
        json.dump(jsons, outfile)


def read_raw_json(raw_json):
    """
    Read raw json file
    :param raw_json:
    :return: list of dictionaries representing reads
    """
    with open(raw_json, 'r') as raw_file:
        json_reads = raw_file.read()

    reads = json.loads(json_reads)
    return reads

def read_part_length(json_log):
    with open(json_log) as json_file:
        json_data = json.load(json_file)
        part_length = json_data["parts_length"]
    return part_length

def write_raw_jasonlog(time, parts_length, filtered_out, log_json_filename):
    """
    write a json log file including time:sec, parts:{id:lenght.....}, junk reads: int
    :param time:
    :param parts_length:
    :param filtered_out:
    :return:
    """
    json_log = {}
    json_log["time"] = time
    json_log["parts_length"] = parts_length
    json_log["junk reads"] = filtered_out
    with open(log_json_filename, 'w') as outfile:
        json.dump( json_log, outfile)


def write_filtered_output_csv(output_filename, reads):
    """
    write output table
    """
    reads_dict = []
    for read in reads:
        read.get_parts_name()
        read.get_identity()
        read_pandas = read.to_dict()
        reads_dict.append(read_pandas)
    output_table = pd.DataFrame(reads_dict)
    output_table.to_csv(output_filename)




