
from nanogate.read import Read
from nanogate.read import Part
# from nanogate.read import PartJaccard
from nanogate.utils import *


# def find_kmer_length(input_parts_library, maximum_range=20):
#     for kmer_length in range(3, maximum_range):
#         kmers = []
#         parts_kmers = []
#         for record in screed.open(input_parts_library):
#             part = Part(record)
#             part.get_sequence()
#             part.kmers_from_read(kmer_length)
#             parts_kmers.append(part.kmers)
#             for kmer in part.kmers:
#                 if kmer not in kmers:
#                     kmers.append(kmer)
#         duplicates_list = []
#         for kmer in kmers:
#             duplicates = 0
#             for part in parts_kmers:
#                 copies = part.count(kmer)
#                 duplicates += copies
#             if duplicates not in duplicates_list:
#                 duplicates_list.append(duplicates)
#         if len(duplicates_list) == 1:
#             return kmer_length
