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

run2(){
  hg38

  align
  sort_index
  bam_rmdup
  bam2bed_fragment
  flank_fragment
  pool_fragment
  add_barcode

  #call_cutesv 1
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


align(){
  init_dir $alignment_dir
  
  if [ ! -s $library.$sample.mapQ30.bam ] || \
     [[ `wc -c $library.$sample.mapQ30.bam | grep -o ^[0-9]*` -lt 30 ]] 
  then
    minimap2 \
      --MD \
      -ax map-ont \
      -t $threads \
      $genome_ref \
      $trim_dir/$library.$sample.trimed.fq.gz \
      2> ./log/$library.$sample.log \
    | samtools view -bS -q 30 - \
    > $library.$sample.mapQ30.bam 
  fi
}

sort_index(){
  cd $alignment_dir
  
  if [ ! -s $library.$sample.mapQ30.sorted.bam ] && \
     [ -s $library.$sample.mapQ30.bam ]
  then
    samtools sort \
      -@ $threads \
      $library.$sample.mapQ30.bam \
      -o $library.$sample.mapQ30.sorted.bam
  fi

  if [ ! -s $library.$sample.mapQ30.sorted.bam.bai ]
  then
    samtools index \
      -@ $threads \
      $library.$sample.mapQ30.sorted.bam
  fi
}

bam_rmdup(){
  cd $alignment_dir

  if [ ! -s $library.$sample.mapQ30.rmdup.sorted.bam ] && \
     [ -s $library.$sample.mapQ30.sorted.bam ]
  then
    samtools rmdup -s \
      $library.$sample.mapQ30.sorted.bam \
      $library.$sample.mapQ30.rmdup.sorted.bam
    
    samtools index \
      -@ $threads \
      $library.$sample.mapQ30.rmdup.sorted.bam
  fi
}

bam2bed_fragment(){
  init_dir $fragment_dir
  
  if [ ! -s $library.$sample.bed ]
  then
    bedtools bamtobed \
      -i $alignment_dir/$library.$sample.mapQ30.rmdup.sorted.bam \
    | awk -vOFS='\t' \
      '{ if($3-$2>=fs) print }' \
      fs=$FILT_SIZE \
    | sort -k1,1 -k2,2n \
    > $library.$sample.bed

    ## add index info
    awk 'NR==FNR{a[$3]=$2;next}{if($4 in a){print $0"\t"a[$4]}}' $index_list $library.$sample.bed > $sample.tmp.bed
    mv $sample.tmp.bed $library.$sample.bed
  fi
}

# flank fragment with read end bias
## add index info
flank_fragment(){
  cd $fragment_dir

  if [ ! -s $library.$sample.flank.bed ]
  then
    cat \
    <(cat $library.$sample.bed \
      | grep -E "^[0-9]|^X|^Y|^chr[0-9]|^chrX|^chrY" \
      | bedtools flank \
          -i - \
          -r 1 \
          -l 0 \
          -g $genome_chr_size \
      | awk -vOFS='\t' '{split($7,a,"_");if($6=="+") {print $1,$2,$3,$4,$5,"+",a[2]} else {print $1,$2,$3,$4,$5,"+",a[1]}}') \
    <(cat $library.$sample.bed \
      | grep -E "^[0-9]|^X|^Y|^chr[0-9]|^chrX|^chrY" \
      | bedtools flank \
          -i - \
          -r 0 \
          -l 1 \
          -g $genome_chr_size \
      | awk -vOFS='\t' '{split($7,a,"_");if($6=="+") {print $1,$2,$3,$4,$5,"-",a[1]} else {print $1,$2,$3,$4,$5,"-",a[2]}}') \
    > $library.$sample.flank.bed
  fi
}

## pesudo bulk
pool_fragment(){
  init_dir $pool_dir
  cd $fragment_dir

  if [ -s $library.$sample.flank.bed ]
  then
    awk -vOFS='\t' '{if($6=="+") {print $1,$2,$2+50,"N","1000",$6,$7} else {print $1,$3-50,$3,"N","1000",$6,$7}}' \
    $library.$sample.flank.bed > $pool_dir/$library.$sample.pool.bed
  fi
}

# run by single cell
add_barcode(){
  cd $alignment_dir

  if [[ -s $library.$sample.mapQ30.rmdup.sorted.bam ]] && \
     ([[ ! -s $library.$sample.barcode.bam ]] || \
      [[ `wc -c $library.$sample.barcode.bam | grep -o ^[0-9]*` -lt 30 ]])
  then
    cat <( samtools view $library.$sample.mapQ30.rmdup.sorted.bam -H ) \
        <( samtools view $library.$sample.mapQ30.rmdup.sorted.bam \
           | awk \
              -vOFS='\t' \
              '{ print $0,"CB:Z:"cell }' \
              cell="${library}#${sample}" ) \
    | samtools view -bS - \
    > $library.$sample.barcode.bam
  fi

  if [[ ! -s $library.$sample.barcode.bam.bai ]]
  then
    samtools index \
      -@ $threads \
      $library.$sample.barcode.bam
  fi
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

run2
