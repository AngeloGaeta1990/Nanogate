import sourmash
import screed

"""
Class to handle reads object
"""


class Read():
    def __init__(self):
        # record contains all the information contatined in a fasta or fastaq file
        self.record = None
        # sequence of a read
        self.sequence = None
        # header of the fasta or fastaq file
        self.features = None
        # id related to read obtained from header
        self.id = None
        # name of the source construct sequenced to obtain a read obtained from header
        self.name = None
        # strand related to the read  obtained from header
        self.strand = None
        # length of sequence contained in the read obtained from header
        self.length = None
        # identity with source sequence obtained from header
        self.readidentity = None
        # length of each k-mer
        self.kmer_length = None
        # list of hashes
        self.hashes = None
        # list of scaled and filtered ashes
        self.filtered_hashes = []
        # list of all unfiltered parts object object
        self.unfiltered_parts = []
        # list of all parts_jaccard objects belonging to the construct
        self.parts = []
        # list of parts name found by the algorithm
        self.parts_name = None
        # list jaccard cofficient related to parts found by the algorithm
        self.jaccards = None
        # dictionary used to print the json file
        self.json = { }
        # dictionary parts:jaccard
        self.parts_jaccard = None
        # dictionary parts:length
        self.parts_length = None

    def get_sequence(self, record):
        """
        method to get the sequence of a read object
        """
        self.record = record
        self.sequence = record.sequence

    def __get_features_for_reads(self):
        """
        Method features: id,
        source construct, strand, sequence length
        and read identity from record of a fasta file
        can be applied only to reads
        """
        self.features = self.record.name
        splitname = self.features.split(" ")
        self.id = splitname[0]
        doublesplit = splitname[1].split(",")
        self.name = int(doublesplit[0])
        self.strand = doublesplit[1]
        self.length = len(self.sequence)
        percentage = splitname[4].replace("read_identity=", "")
        self.readidentity = float(percentage.replace("%", "")) / 100

    # def __get_parts(self, table):
    #     """
    #     Method to get all parts included into the construct sequenced to
    #     generate this read applied only to read object
    #     """
    #     row = table.loc[table['Construct'] == self.construct_name]
    #     parts = row["Parts"].values[0]
    #     parts_list = parts.split(" ")
    #     self.parts_verification = parts_list

    # def kmers_from_read(self, ksize):
    #     """
    #     method to get all k-mers of read given read sequence
    #     and the k-mer length
    #     """
    #     self.kmer_length = ksize
    #     n_kmers = len(self.sequence) - ksize + 1
    #
    #     for i in range(n_kmers):
    #         kmer = self.sequence[i:i + ksize]
    #         self.kmers.append(kmer)

    def hash_from_kmers(self, kmer_size):
        """
        Method to turn all kmers into hashes code
        """
        max_kmer_number = (self.length - kmer_size + 1)
        if max_kmer_number > 0:
            hash_object = sourmash.MinHash(n=max_kmer_number, ksize=kmer_size)
            hash_object.add_sequence(self.sequence)
            self.hashes = hash_object

    def get_parts_name(self):
        """
        Method to go trough list of parts associated to the construct and get name
        """
        parts_name_list = []
        for part in self.parts:
                parts_name_list.append(part.name)
        self.parts_name = str(parts_name_list).strip('[]')
    #
    def get_identity(self):
        """
        Method to go trough list of parts associated to the construct and get the jaccard coefficient
        """
        part_jaccards_list = []
        for part in self.parts:
            part_jaccards_list .append(part.jaccard)
        self.jaccards = str(part_jaccards_list).strip('[]')

    def to_dict(self):
        """
        method to return a dictionary representation of the class

        :return:
        """
        return {
            'Construct': self.name,
            'Parts': self.parts_name,
            'Identity': self.jaccards
        }

    def get_json(self):
        """
        Method to initialise json content
        :return:
        """

        self.json["Construct id"] = self.name
        self.json["read length"] = self.length
        self.json["part_jaccard"] = {}
        self.json["part_length"] = {}

    def add_jaccard_to_json(self, part, jaccard):
        self.json["part_jaccard"][part] = jaccard

    def add_length_to_json(self,part,length):
        self.json["part_length"][part] = length


    def build_read_from_dictionary(self, read_dictionary):
        """
        from json dictionrary to read object
        :param read_dictionary:
        :return:
        """
        self.name = read_dictionary["Construct id"]
        self.length = read_dictionary["read length"]
        self.parts_jaccard = read_dictionary["part_jaccard"]
        self.parts_length = read_dictionary["part_length"]



    def initialise(self, record, kmer_size):
        """
        Method to initialise Read object
        """
        self.get_sequence(record)
        self.__get_features_for_reads()
        self.hash_from_kmers(kmer_size)

    # def filter_hashes(self, scale):
    #     """
    #     method to filter kmers on the basis of hash_maxvalue/scale
    #     only hashes with a value <  hash_maxvalue/scale are kept
    #     """
    #     max_hash = 2 ** 64
    #     keep_below = max_hash / scale
    #     for hash in self.hashes:
    #         if hash < keep_below:
    #             # otherwise, discard
    #             self.filtered_hashes.append(hash)


class Part(Read):
    def __init__(self):
        Read.__init__(self)

        # jaccard containment value evaluted after comparison with the construct
        self.jaccard = None

    def get_features_for_part(self):
        """
        Method to get the ID of a part
        applied only on parts
        """
        self.features = self.record.name
        splitname = self.features.split(" ")
        self.name = splitname[0]

    def get_length(self):
        """
        get the length of a sequence
        """
        self.length = len(self.sequence)

    def build_part_from_json(self,part_name,part_length,read):
        self.name = int(part_name)
        self.jaccard = read.parts_jaccard[part_name]
        self.length = part_length[part_name]



    def initialise_part(self, record, kmer_size):
        """
        Method to initialise Part object
        """
        self.get_sequence(record)
        self.get_features_for_part()
        self.get_length()
        self.hash_from_kmers(kmer_size)


# """
# Class to store in memory only part name and jaccard value
# related to a construct
# """
#
#
# class PartJaccard():
#     def __init__(self, part_name, length, jaccard):
#         # name of the part
#         self.part_name = part_name
#         # length of a part
#         self.length = length
#         # jaccard cofficient evaluated between the part and the construct
#         self.jaccard = jaccard
