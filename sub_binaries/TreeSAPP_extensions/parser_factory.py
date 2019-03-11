import sys
import time
import re
import logging
from distutils.core import setup,Extension
from parsers import parsefile


def match_file_to_dict(file_handler, key_dict, sep="\t", join_by=0):
    """
    Generator function for mapping a particular field in a file, separated by 'sep', to dictionary keys
    :param file_handler: Opened file object
    :param key_dict: Dictionary to map the selected field to
    :param sep: Field separator
    :param join_by: The field number to search the dictionary for
    :return: Line that matches a key
    """
    for line in file_handler:
        if line.split(sep)[join_by] in key_dict:
            yield line


class ReferenceSequence:
    def __init__(self):
        self.accession = ""
        self.description = ""
        self.organism = ""
        self.lineage = ""
        self.short_id = ""
        self.sequence = ""
        self.locus = ""
        self.cluster_rep = False
        self.cluster_rep_similarity = 0
        self.cluster_lca = None

    def get_info(self):
        """
        Returns a string with the ReferenceSequence instance's current fields

        :return: str
        """
        info_string = ""
        info_string += "accession = " + self.accession + ", " + "mltree_id = " + self.short_id + "\n"
        info_string += "description = " + self.description + ", " + "locus = " + self.locus + "\n"
        info_string += "organism = " + self.organism + "\n"
        info_string += "lineage = " + self.lineage + "\n"
        return info_string

class EntrezRecord(ReferenceSequence):
    def __init__(self, acc, ver):
        super().__init__()
        self.accession = acc
        self.versioned = ver
        self.ncbi_tax = ""
        self.bitflag = 0  # For monitoring progress during download stage


def c_map_accession2taxid(query_accession_list, accession2taxid_list):
    er_acc_dict = dict()
    unmapped_queries = list()

    result_list = list()
    # Create a dictionary for O(1) look-ups, load all the query accessions into unmapped queries
    for acc in query_accession_list:  # type: str
        if acc.find('.') >= 0:
            ver = acc
            # Strip off any version numbers from the accessions so we only need to check for one item
            acc = '.'.join(acc.split('.')[0:-1])
        else:
            ver = ""
        er_acc_dict[acc] = EntrezRecord(acc, ver)
        unmapped_queries.append(acc)

    logging.info("Mapping query accessions to NCBI taxonomy IDs... ")
    for accession2taxid in accession2taxid_list.split(','):
        init_qlen = len(unmapped_queries)
        final_qlen = len(unmapped_queries)

        start = time.time()

        ## Call C extension
        result_list = parsefile(accession2taxid, query_accession_list)
        
        end = time.time()
        for i in range(len(result_list[0])):
            try:
                # Update the EntrezRecord elements
                accession = result_list[0][i]
                record = er_acc_dict[accession]
                record.versioned = result_list[1][i]
                record.ncbi_tax = result_list[2][i]
                record.bitflag = 3  # Necessary for downstream filters - indicates taxid has been found
                # Remove accession from unmapped queries
                i = 0
                while i < final_qlen:
                    if unmapped_queries[i] == accession:
                        unmapped_queries.pop(i)
                        final_qlen -= 1
                        break
                    i += 1
                if final_qlen == 0:
                    break
            except KeyError:
                logging.error("Bad key returned by generator.\n")
                sys.exit(13)
                
        print("Time required to parse '" + accession2taxid + "': " + str(end - start) + "s.\n")
        # Report the number percentage of query accessions mapped
        print(
            str(round(((init_qlen - final_qlen) * 100 / len(query_accession_list)), 2)) +
            "% of query accessions mapped by " + accession2taxid + ".\n")
    logging.info("done.\n")

    return er_acc_dict

def map_accession2taxid(query_accession_list, accession2taxid_list):
    er_acc_dict = dict()
    unmapped_queries = list()

    # Create a dictionary for O(1) look-ups, load all the query accessions into unmapped queries
    for acc in query_accession_list:  # type: str
        if acc.find('.') >= 0:
            ver = acc
            # Strip off any version numbers from the accessions so we only need to check for one item
            acc = '.'.join(acc.split('.')[0:-1])
        else:
            ver = ""
        er_acc_dict[acc] = EntrezRecord(acc, ver)
        unmapped_queries.append(acc)

    logging.info("Mapping query accessions to NCBI taxonomy IDs... ")
    for accession2taxid in accession2taxid_list.split(','):
        init_qlen = len(unmapped_queries)
        final_qlen = len(unmapped_queries)
        start = time.time()
        try:
            rosetta_handler = open(accession2taxid, 'r')
        except IOError:
            logging.error("Unable to open '" + accession2taxid + "' for reading.\n")
            sys.exit(13)

        for line_match in match_file_to_dict(rosetta_handler, er_acc_dict):
            try:
                accession, ver, taxid, _ = line_match.strip().split("\t")
            except (ValueError, IndexError):
                logging.warning("Parsing '" + accession2taxid + "' failed.\n")
                break

            try:
                # Update the EntrezRecord elements
                record = er_acc_dict[accession]
                record.versioned = ver
                record.ncbi_tax = taxid
                record.bitflag = 3  # Necessary for downstream filters - indicates taxid has been found
                # Remove accession from unmapped queries
                i = 0
                while i < final_qlen:
                    if unmapped_queries[i] == accession:
                        unmapped_queries.pop(i)
                        final_qlen -= 1
                        break
                    i += 1
                if final_qlen == 0:
                    break
            except KeyError:
                logging.error("Bad key returned by generator.\n")
                sys.exit(13)

        rosetta_handler.close()
        end = time.time()
        print("Time required to parse '" + accession2taxid + "': " + str(end - start) + "s.\n")
        # Report the number percentage of query accessions mapped
        print(
            str(round(((init_qlen - final_qlen) * 100 / len(query_accession_list)), 2)) +
            "% of query accessions mapped by " + accession2taxid + ".\n")
    logging.info("done.\n")

    return er_acc_dict

if __name__ == "__main__":
    print("starting...")
        
    query_accession_list = ["X53814", "X63317", "X63826","X61113", "X51426", "X51425", "X60058", "X60062", "X58230", "X58229", "X62343", "X64519", "X15224"]
    accession2taxid_list = "accession2taxid_data.txt"

    py_results = map_accession2taxid(query_accession_list, accession2taxid_list)
    print(py_results)
    
    results = c_map_accession2taxid(query_accession_list, accession2taxid_list)
    print("C RESULTS: \n")
    print(results)




