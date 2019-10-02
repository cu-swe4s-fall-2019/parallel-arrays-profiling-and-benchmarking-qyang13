import sys
import gzip
import matplotlib
import matplotlib.pyplot as plt
import data_viz as dv
matplotlib.use('Agg')

def linear_search(key, L):
    hit = -1
    for i  in range(len(L)):
        curr =  L[i]
        if key == curr:
            return i
    return -1


def binary_search(key, D):
    lo = -1
    hi = len(D)
    while (hi - lo > 1):
        mid = (hi + lo) // 2

        if key == D[mid][0]:
            return D[mid][1]

        if ( key < D[mid][0] ):
            hi = mid
        else:
            lo = mid

    return -1


def linear_search_all_hits (key, L):
    '''
    Gives indexs not values
    '''
    hit = []
    for i in range(len(L)):
        if key == L[i]:
            hit.append(i)
    return hit


def main():
    meta_data_file_name = "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
    rna_data_file_name = "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.acmg_59.gct.gz"
    target_type = "SMTS"
    target_gene_name = "ACTA2"
    out_file_name = target_gene_name + "_" + target_type + "_boxplot.png"


    #meta_data_file_name = args.meta
    #rna_data_file_name = args.counts
    #tissue_type = args.tissue
    #target_gene_name = args.gene

    SAMPID = []
    SMTS = []
    metadata_header = None
    sample_info = None
    tissue_group = []
    xlabels = []
    # This is a 2D array, first layer stores tissue groups, second stores
    # all the sample IDs that belong to that group
    categoraized_ids = []
    # Processing the meta data file
    for l in open(meta_data_file_name):
        # Store the first row i.e. header in a separate array
        if metadata_header is None:
            metadata_header = l.rstrip().split("\t")
            SAMPID_idx = linear_search("SAMPID", metadata_header)
            SMTS_idx = linear_search(target_type, metadata_header)
        else:
            # Split lines by tab and stored them in array
            # rstrip() returns a copy of the string with trailing characters removed
            sample_info = l.rstrip().split("\t")
            sample_id = sample_info[SAMPID_idx]
            tissue_type = sample_info[SMTS_idx]
            # Try to find the tissue type in the existing list
            tissue_idxs = linear_search(tissue_type, tissue_group)
            # If not found, add that tissue type to the list
            if tissue_idxs == -1:
                xlabels.append(tissue_type)
                tissue_idxs = len(tissue_group)
                tissue_group.append(tissue_type)
                categoraized_ids.append([])
            categoraized_ids[tissue_idxs].append(sample_id)


    # For unit test, no need to be random
    # Test if the return values match, and the length
    # Also make sure the headers are there
    # print(SAMPID[0],SMTS[0],SMTSD[0])

    version = None
    dim = None
    rna_header = None
    counts = [[] for i in range(len(tissue_group))]

    # Processing the count file
    # use gzip.open to read gzip file
    for l in gzip.open(rna_data_file_name,"rt"):
        # Assume first line stores the version
        if version == None:
            version = l
            continue
        # Assume second line stores the dimension
        if dim == None:
            dim = l
            continue
        # Assume thrid line stores the header for counts
        if rna_header == None:
            rna_header = l.rstrip().split("\t")
            # Use tuple to store the original index
            # and then sort based on the first element in tuple
            rna_header_plus_index = []
            for i in range(len(rna_header)):
                rna_header_plus_index.append([rna_header[i],i])
            rna_header_plus_index.sort(key=lambda pair: pair[0])
            # Store the index of description
            Description_idx = binary_search("Description",
                                            rna_header_plus_index)
            if Description_idx == -1:
                sys.write('Description not found in header.')

        else :
            rna_counts = l.rstrip().split("\t")

            if rna_counts[Description_idx] == target_gene_name:
                # For each tissue type
                for tissue_idx in range(len(tissue_group)):
                    # For each individual in the same tissue type
                    for sample in categoraized_ids[tissue_idx]:
                        rna_header_idx = binary_search(sample, rna_header_plus_index)
                        # Some of the sampeles do not have rna counts
                        if rna_header_idx != -1:
                            count = int(rna_counts[rna_header_idx])
                            counts[tissue_idx].append(count)
                break

    dv.boxplot(counts, xlabels, target_type, target_gene_name, out_file_name)

if __name__ == '__main__':
    main()
