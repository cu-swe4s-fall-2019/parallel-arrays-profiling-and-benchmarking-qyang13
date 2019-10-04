# Parallel Arrays, Profiling, and Benchmarking

# Usage
Before start, you must download the following files for analysis:
- https://github.com/swe4s/lectures/blob/master/data_integration/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.acmg_59.gct.gz?raw=true
- https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt

Since the `.gz` file is named weird, be sure to remove the additional suffix:
```
mv GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.acmg_59.gct.gz?raw=true GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.acmg_59.gct.gz
```
Also, make sure the matplotlib is installed:
```
conda install --yes matplotlib
```
## Testing
To run unit test for individual functions implemented in plot_gtex, run:
```
python test_plot_gtex.py
```
To run function test to ensure everything is installed properly, run:
```
bash test_plot_gtex.sh
```

## Run the program
To plot boxplot of expression distribution across tissue types for a specific gene, you can run the python script `plot_gtex.py`:

```
python  plot_gtex.py \
--gene_reads \
GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.acmg_59.gct.gz \
--sample_attributes GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt \
--gene ACTA2 \
--group_type SMTS \
--out_file ACTA2.png
```

Where the required arguments are:

| Argument   | Definition   |
| --- | ---|
| --gene_reads | A read count files where columns are samples and rows are genes |
| --sample_attributes | The meta data file |
| --gene | Your gene of interest |
| --group_type | Your tissue group of interest <`SMTS` or `SMTSD`> |
| --out_file | The filename for the output plot |

The resulting plot should look like this:
![alt text](https://github.com/cu-swe4s-fall-2019/parallel-arrays-profiling-and-benchmarking-qyang13/blob/documentation/ACTA2.png "example plot")


# Profiling
The program is profiled using python cProfile, the output is sorted by the total runtime of each process or function call.

As shown in `plot_gtex.linear_serach.txt`, the most consuming function call is linear search, which took a total of `16.273` seconds:
```
ncalls  tottime  percall  cumtime  percall filename:lineno(function)
 45905   16.273    0.000   16.277    0.000 plot_gtex.py:10(linear_search)
```
Thus to speed up the process, binary search is implemented. As a result, the program ran much faster, wich binary search step taking up a total of `0.085` seconds.
```
ncalls  tottime  percall  cumtime  percall filename:lineno(function)
     2    0.266    0.133    0.266    0.133 {built-in method matplotlib._png.write_png}
   250    0.100    0.000    0.105    0.000 <frozen importlib._bootstrap_external>:914(get_data)
 22952    0.085    0.000    0.086    0.000 plot_gtex.py:20(binary_search)
```
