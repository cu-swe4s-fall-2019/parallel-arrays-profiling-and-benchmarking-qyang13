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
    rna_data_file_name="GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.acmg_59.gct.gz"

    SAMPID = []
    SMTS = []
    SMTSD = []

    metadata_header = None

    for l in open(meta_data_file_name):
    	# Split lines by tab and stored them in array
    	# rstrip() returns a copy of the string with trailing characters removed
    	sample_info = l.rstrip().split("\t")

    	# Store the first row i.e. header in a separate array
    	if metadata_header is None:
    		metadata_header = sample_info
    		continue

    	SAMPID_idx = linear_search("SAMPID", metadata_header)
    	SMTS_idx = linear_search("SMTS", metadata_header)
    	SMTSD_idx = linear_search("SMTSD", metadata_header)

    	#print (sample_info[SAMPID_idx],sample_info[SMTS_idx],sample_info[SMTSD_idx])

    	SAMPID.append(sample_info[SAMPID_idx])
    	SMTS.append(sample_info[SMTS_idx])
    	SMTSD.append(sample_info[SMTSD_idx])

    # For unit test, no need to be random
    # Test if the return values match, and the length
    # Also make sure the headers are there
    # print(SAMPID[0],SMTS[0],SMTSD[0])

    version = None
    dim = None
    rna_header = None
    target_gene_name = "SDHC"

    # use gzip.open to read gzip file
    for l in gzip.open(rna_data_file_name,"rt"):
    	if version == None:
    		version = l
    		continue

    	if dim == None:
    		dim = l
    		continue

    	if rna_header == None:
    		rna_header = l.rstrip().split("\t")
    		continue

    	# Use tuple to store the original index and then sort based on the first element in tuple, ie the header
    	rna_header_plus_index = []
    	for i in range(len(rna_header)):
    		rna_header_plus_index.append([rna_header[i],i])
    	rna_header_plus_index.sort(key=lambda pair: pair[0])

    	rna_counts = l.rstrip().split("\t")

    	Description_idx = linear_search("Description", rna_header)

    	if Description_idx == -1:
    		sys.write('Description not found in header.')

    	if rna_counts[Description_idx] == target_gene_name:
    		# Obtain the index in the metadata to get all the names
    		blood_idxs = linear_search_all_hits('Blood',SMTS)

    		# Extract the read counts from rna_counts
    		blood_counts = []
    		for blood_idxs in blood_idxs:
    			#rna_header_idx = linear_search(SAMPID[blood_idxs],rna_header)
    			rna_header_idx = binary_search(SAMPID[blood_idxs], rna_header_plus_index)
    			# Some of the sampeles do not have rna counts
    			if rna_header_idx != -1:
    				count = int(rna_counts[rna_header_idx])
    				blood_counts.append(count)
    		break


    dv.boxplot(blood_counts, 'test.png')

if __name__ == '__main__':
    main()
