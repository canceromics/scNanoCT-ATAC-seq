# scNanoCT-ATAC-seq
We developed Nanopore sequencing based single-cell CUT&Tag and ATAC co-profiling (scNanoCT&ATAC-seq), which allows for the simultaneous mapping of genome mutations, histone modifications and chromatin accessibility within the same individual cells.


# Requirements
- Please install the following software before running the scripts.  
- It is recommended to create a new virtual environment (e.g., with Conda) for project analysis.
```
python  3.8.5
GNU parallel  20250322
R  4.2.0
ArchR  1.0.1
seqkit  2.6.1
nanoplexer      0.1
cutadapt        4.1
minimap2        2.26
samtools        1.13
bedtools        2.30.0
cuteSV  1.0.10
macs2   2.2.9.1
SEACR   1.3
FREEC   11.6b
bedGraphToBigWig
deepTools  3.5.5
```



# Inputs
- Please ensure that the input files and information are correctly provided.
- PASS FASTQ data (QS>7) should be merged into one file and placed in the `raw_data` directory, and other files should be placed in the `root_dir` directory.

| Name | Description |
| --- | --- |
| ONT_lib | A custom library name (e.g., test). Demultiplexing will be performed in a folder with that name.  |
| threads | Number of cores used in the parallel analysis step. |
| root_dir | Project root directory (the path containing scripts, data and analysis results). |
| bc_index | A file records the barcode indexes used in the library. |
| his_cell | A file records the custom group name (e.g., Rep1/Rep2) and their corresponding barcodes indexes. All analysis outputs will be saved in a folder named accordingly. |
| IBC | A file records inner barcode index and sequence. |
| OBC | A file records outer barcode index and sequence. |



# Usage
- Each directory corresponds to an analysis module. Please follow the step numbers and run them sequentially.
- It is recommended to use a multi-core server or HPC cluster for parallel processing, as most analyses after demultiplexing are carried out per cell.

## 1_Data_preprocessing
The first part of the code handles data demultiplexing and comprises four steps:
1. `load_dual`: generates a standard FASTA file from the barcode file (e.g., barcode/test/Findex.fa).
2. `outer_barcode` and `inner_barcode`: perform two rounds of demultiplexing to produce single‑cell FASTQ files (e.g., barcode/test/Rep1/1_10.fastq).
3. `dual_index`: identifies and trims the ATAC (I7) and histone‑modification (I5) adapters, producing cleaned FASTQ files (e.g., trim/test/Rep1/1_10.fastq).
4. Adapter information for each cell and each read is recorded in files named cell_list and index_list.  

```
bash $root_dir/1_pipeline.sh $ONT_lib $threads $root_dir $bc_index $his_cell
```

**The analysis outputs will be saved in the** `barcode` **and** `trim` **directory.**


## 2_Signal_extraction
The second part of the code performs data alignment and omics signal extraction, consisting of five steps:
1. `hg38`: indexes generation for the reference genome of human or other species.
2. `align`, `sort_index`, and `bam_rmdup`: aligns reads to the reference genome, sorts alignments by coordinate, and removes PCR duplicates, respectively.
3. `bam2bed_fragment` and `flank_fragment`: extract alignment endpoints in genomic coordinates and assign them to distinct omics modalities based on index_list.
4. `pool_fragment`: merges genomic fragments from all single cells for each modality to generate pseudo‑bulk omics data.
5. `add_barcode`: adds the cell barcode and adapter information to the alignment files in the ‘CB:Z’ tag.

### Run for individual cell
```
bash $root_dir/2_pipeline.sh $library $cell $threads $root_dir
```
### Parallel run in multi-core server
```
awk '{print "bash '$root_dir'/2_pipeline.sh '$library' "$1" 1 '$root_dir'"}' $root_dir/barcode/$library/cell_list > cell_task
cat cell_task | parallel -j $threads
```
### Parallel run in HPC cluster
```
for cell in `cut -f1 $root_dir/barcode/$library/cell_list`;do
  id2=$( sbatch -J $cell -p cpuPartition -c 1 --parsable $root_dir/2_pipeline.sh $library $cell 1 $root_dir )
done
```

**The analysis outputs will be saved in the** `alignment`, `fragment`, **and** `pool` **directory.**


## 3_Peak_calling

## 4_SV_detection

## 5_CNV_detection

# Contact
Please go to the issues page for help or contact with us directly by Xiaoying Fan (fan_xiaoying@gzlab.ac.cn) and Zhixiang Zuo (zuozhx@sysucc.org.cn).
