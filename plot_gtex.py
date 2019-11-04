import sys
import time
import argparse as ap
import gzip
import matplotlib
import matplotlib.pyplot as plt
import data_viz as dv
sys.path.insert(1, "./hash-tables-qyang13/")
import hash_tables as ht
import hash_functions as hf

matplotlib.use('Agg')

def linear_search(key, L):
    '''Linear search to find the key in a list'''
    hit = -1
    for i in range(len(L)):
        curr = L[i]
        if key == curr:
            return i
    return -1


def binary_search(key, D):
    '''Binary search to find the key in a SORTED list'''
    lo = -1
    hi = len(D)
    while(hi - lo > 1):
        mid = (hi + lo) // 2

        if key == D[mid][0]:
            return D[mid][1]

        if(key < D[mid][0]):
            hi = mid
        else:
            lo = mid

    return -1


def parse_args():
    '''    Argument Parser    '''
    parser = ap.ArgumentParser(description="correct way to parse",
                               prog='Plot GTEX')

    parser.add_argument('-o',
                        '--out_file',
                        type=str,
                        help="Output filename",
                        required=False)

    parser.add_argument('-a',
                        '--sample_attributes',
                        type=str,
                        help="Input meta data filename",
                        required=False)

    parser.add_argument('-c',
                        '--gene_reads',
                        type=str,
                        help="Input read counts filename",
                        required=False)

    parser.add_argument('-t',
                        '--group_type',
                        type=str,
                        help="Group type: SMTS or SMTSD",
                        required=False)

    parser.add_argument('-g',
                        '--gene',
                        type=str,
                        help="Enter a gene name",
                        required=False)

    parser.add_argument('-k',
                        '--hash',
                        action='store_true',
                        help="Enable hash",
                        required=False)

    return parser.parse_args()


def main():
    main_start = time.time()
    args = parse_args()
    meta_data_file_name = args.sample_attributes
    rna_data_file_name = args.gene_reads
    target_type = args.group_type
    target_gene_name = args.gene
    out_file_name = args.out_file

    # Debug mode
    meta_data_file_name = 'GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt'
    rna_data_file_name = 'GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.acmg_59.gct.gz'
    target_type = 'SMTS'
    target_gene_name = 'ACTA2'
    out_file_name = 'ACTA2.png'

    # For unit test, no need to be random
    # Test if the return values match, and the length
    # Also make sure the headers are there
    # print(SAMPID[0],SMTS[0],SMTSD[0])

    SAMPID = []
    SMTS = []
    metadata_header = None
    tissue_group = []
    xlabels = []

    if args.hash!=True:
        # This is a 2D array, first layer stores tissue groups, second stores
        # all the sample IDs that belong to that group
        categoraized_ids = []
        sample_info = None

        # Processing the meta data file
        for l in open(meta_data_file_name):
            # Store the first row i.e. header in a separate array
            if metadata_header is None:
                metadata_header = l.rstrip().split("\t")
                SAMPID_idx = linear_search("SAMPID", metadata_header)
                SMTS_idx = linear_search(target_type, metadata_header)
            else:
                # Split lines by tab and stored them in array
                # rstrip() returns a copy of the string
                # with trailing characters removed
                sample_info = l.rstrip().split("\t")
                sample_id = sample_info[SAMPID_idx]
                tissue_type = sample_info[SMTS_idx]
                # Try to find the tissue type in the existing list
                tissue_idxs = linear_search(tissue_type, tissue_group)
                # If not found, add that tissue type to the list
                if tissue_idxs == -1:
                    tissue_idxs = len(tissue_group)
                    tissue_group.append(tissue_type)
                    categoraized_ids.append([])
                categoraized_ids[tissue_idxs].append(sample_id)
        counts = [[] for i in range(len(tissue_group))]

    # Hash implementation of plot_gtex
    else:  # Hash enabled
        # Hash stores tissue type as key and sample ids as value
        categorized_hash = ht.ChainHashTable(100, hf.h_sedgwicks)
        sample_info = []
        # Processing the meta data file
        for l in open(meta_data_file_name):
            # Store the first row i.e. header in a separate array
            if metadata_header is None:
                metadata_header = l.rstrip().split("\t")
                SAMPID_idx = linear_search("SAMPID", metadata_header)
                SMTS_idx = linear_search(target_type, metadata_header)
            else:
                sample_info.append(l.rstrip().split("\t"))
        for i in sample_info:
            if categorized_hash.search(i[SMTS_idx]) is None:
                categorized_hash.add(i[SMTS_idx],i[SAMPID_idx])
            else:
                hash_slot = categorized_hash.search_slot(i[SMTS_idx])
                categorized_hash.T[hash_slot][0][1].append(i[SAMPID_idx])
        counts = [[] for i in range(len(categorized_hash.K))]


    version = None
    dim = None
    rna_header = None

    if args.hash!=True:
        # Processing the count file
        # use gzip.open to read gzip file
        for l in gzip.open(rna_data_file_name, "rt"):
            # Assume first line stores the version
            if version is None:
                version = l
                continue
            # Assume second line stores the dimension
            if dim is None:
                dim = l
                continue
            # Assume thrid line stores the header for counts
            if rna_header is None:
                rna_header = l.rstrip().split("\t")
                # Use tuple to store the original index
                # and then sort based on the first element in tuple
                rna_header_plus_index = []
                sort_start = time.time()
                for i in range(len(rna_header)):
                    rna_header_plus_index.append([rna_header[i], i])
                rna_header_plus_index.sort(key=lambda pair: pair[0])
                sort_end = time.time()
                # Store the index of description
                # Description_idx = linear_search("Description",
                #                                 rna_header)
                Description_idx = binary_search("Description",
                                                rna_header_plus_index)
                if Description_idx == -1:
                    sys.exit('Description not found in header.')

            else:
                rna_counts = l.rstrip().split("\t")
                if rna_counts[Description_idx] == target_gene_name:
                    search_start = time.time()
                    # For each tissue type
                    for tissue_idx in range(len(tissue_group)):
                        # For each individual in the same tissue type
                        for sample in categoraized_ids[tissue_idx]:
                            # rna_header_idx = linear_search(sample,
                            #                                rna_header)
                            rna_header_idx = binary_search(sample,
                                                           rna_header_plus_index)
                            # Some of the sampeles do not have rna counts
                            if rna_header_idx != -1:
                                count = int(rna_counts[rna_header_idx])
                                counts[tissue_idx].append(count)
                    search_end = time.time()
                    break
        dv.boxplot(counts, tissue_group, target_type, target_gene_name, out_file_name)

    # Hash implementation
    else:
        # Second hash table used to store counts data for each SAMPID
        counts_hash = ht.ChainHashTable(100000, hf.h_sedgwicks)

        for l in gzip.open(rna_data_file_name, "rt"):
            # Assume first line stores the version
            if version is None:
                version = l
                continue
            # Assume second line stores the dimension
            if dim is None:
                dim = l
                continue
            # Assume thrid line stores the header for counts
            if rna_header is None:
                rna_header = l.rstrip().split("\t")
                # Use tuple to store the original index
                # and then sort based on the first element in tuple
                rna_header_plus_index = []
                # sort_start = time.time()
                for i in range(len(rna_header)):
                    rna_header_plus_index.append([rna_header[i], i])
                rna_header_plus_index.sort(key=lambda pair: pair[0])
                # sort_end = time.time()
                # Store the index of description
                # Description_idx = linear_search("Description",
                #                                 rna_header)
                Description_idx = binary_search("Description",
                                                rna_header_plus_index)
                if Description_idx == -1:
                    sys.exit('Description not found in header.')

            else:
                rna_counts = l.rstrip().split("\t")
                sampid_idx = Description_idx + 1
                if rna_counts[Description_idx] == target_gene_name:
                    for count_header, count in\
                     zip(rna_header[sampid_idx:], rna_counts[sampid_idx:]):
                        # Adding samples and its corresponding gene counts
                        counts_hash.add(count_header, int(count))
        # For each tissue type
        for i in range(len(counts)):
            # For each individual in the same tissue type
            search_start = time.time()
            for sample_id in categorized_hash.search(categorized_hash.K[i]):
                sample_count = counts_hash.search(sample_id)
                if sample_count is not None:
                    counts[i].append(sample_count[0])
            search_end = time.time()
        dv.boxplot(counts, categorized_hash.K, target_type, target_gene_name, out_file_name)
    main_end = time.time()

    # print('Sorting time : ' + str(sort_end - sort_start))
    print('Searching time : ' + str(search_end - search_start))
    print('Main time: ' + str(main_end - main_start))







if __name__ == '__main__':
    main()
