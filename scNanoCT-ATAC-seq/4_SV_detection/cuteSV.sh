#!/bin/bash

library=$1
sample=$2
threads=$3
root_dir=$4

trim_dir=$root_dir/trim/$library/
fragment_dir=$root_dir/fragment/$library/
pool_dir=$root_dir/pool/$library/

index_list=$root_dir/barcode/$library/index_list

FILT_SIZE=100
#FILT_SIZE=1000
export PATH=/home/linzhuobin/software/minimap2-2.26_x64-linux:$PATH
eval "$(conda shell.bash hook)"
conda activate scNanoATAC-env

run(){
  hg38

  call_cutesv 1
}

# human reference
hg38(){
  genome_ref=/home/linzhuobin/reference/human/GRCh38_v44/GRCh38.primary_assembly.genome.mmi
  genome_chr_size=/home/linzhuobin/reference/human/GRCh38_v44/hg38.chrom.sizes
  ref_genome_fa=/home/linzhuobin/reference/human/GRCh38_v44/GRCh38.primary_assembly.genome.fa
  g_size=hs
  genome_name=hg38

  alignment_dir=$root_dir/alignment/$library/
}

init_dir(){
  mkdir -p $1/log/
  cd $1
}

# run by single cell
call_cutesv(){
  local r=$1
  
  cutesv_dir=$root_dir/cutesv/${r}r/$library/
  init_dir $cutesv_dir
  
  if [ ! -s $library.$sample.cuteSV.ONT.vcf ] ; then
    echo "call sv $library.$sample"
    init_dir $cutesv_dir/$sample/
    cuteSV \
      $alignment_dir/$library.$sample.mapQ30.sorted.bam \
      $ref_genome_fa                          \
      $library.$sample.cuteSV.ONT.vcf         \
      .                                       \
      -t $threads                             \
      --max_cluster_bias_INS	100             \
      --diff_ratio_merging_INS	0.3           \
      --max_cluster_bias_DEL  	100           \
      --diff_ratio_merging_DEL	0.3           \
      --min_support $r                        \
      --sample ${library}.${sample} 
    
    cd $cutesv_dir
    mv ./$sample/$library.$sample.cuteSV.ONT.vcf ./
    rm -rf ./$sample/

    cat $library.$sample.cuteSV.ONT.vcf \
    | sed 's/chr//g' \
    > $library.$sample.cuteSV.ONT.vcf.clean

    mv \
      $library.$sample.cuteSV.ONT.vcf.clean \
      $library.$sample.cuteSV.ONT.vcf
  fi
}

run
