#! /bin/bash

library=$1
threads=$2
root_dir=$3

archr_script=$root_dir/archr.R

fragment_dir=$root_dir/fragment/$library/
archr_dir=$root_dir/archr/$library/
pool_dir=$root_dir/pool/$library/

export PATH=/home/linzhuobin/software/nanoplexer:$PATH
export PATH=/home/linzhuobin/software/minimap2-2.26_x64-linux:$PATH
eval "$(conda shell.bash hook)"
conda activate scNanoATAC-env

run3(){
  hg38

  archr_fragment_file
  create_arrow_file
  pool_fragment_file
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

# run by library
archr_fragment_file(){
  cd $(dirname $fragment_dir)

  cell_list=`ls ./$library/ | grep -v log | \
	awk -F "." '{print $2}' | sort | uniq`

  for cell in $cell_list
  do
    cat ./$library/$library.$cell.flank.bed \
    | awk -vOFS='\t' \
      '{
        if ($1!~/chr/) $1="chr"$1
        print $1,$2,$3,cell,1,$7
      }' \
      cell=$cell
  done \
  | sort -k 1,1 -k 2,2n -k 4,4 \
  > $library.fragments.sorted.bed

  ## I5 fragment
  grep -w 'I5' $library.fragments.sorted.bed > $library.I5.sorted.bed
  ## I7 fragment
  grep -w 'I7' $library.fragments.sorted.bed > $library.I7.sorted.bed

  bgzip -f $library.fragments.sorted.bed
  tabix -f $library.fragments.sorted.bed.gz
  bgzip -f $library.I5.sorted.bed
  tabix -f $library.I5.sorted.bed.gz
  bgzip -f $library.I7.sorted.bed
  tabix -f $library.I7.sorted.bed.gz
}

# run by library
# activate r40 first
create_arrow_file(){
  conda deactivate
  conda activate jupyter-env  

  init_dir $archr_dir
  
  Rscript $archr_script \
    $library $(dirname $fragment_dir) \
    $genome_name 
  
  mv $library.arrow ../
  conda deactivate
  conda activate scNanoATAC-env
}

pool_fragment_file(){
  cd $(dirname $pool_dir)

  cell_list=`ls ./$library/ | grep -v log | \
        awk -F "." '{print $2}' | sort | uniq`

  ## all fragment
  cat ./$library/*pool.bed | sort -k 1,1 -k 2,2n -k 4,4 \
  > $library.fragments.sorted.bed
  rm -rf $library

  ## I5 fragment
  grep -w 'I5' $library.fragments.sorted.bed | \
  awk -vOFS='\t' '{print $1,$2,$3,$4,$5,$6}' > $library.I5.sorted.bed
  ## I7 fragment
  grep -w 'I7' $library.fragments.sorted.bed | \
  awk -vOFS='\t' '{print $1,$2,$3,$4,$5,$6}' > $library.I7.sorted.bed

  bgzip -f $library.fragments.sorted.bed
  tabix -f $library.fragments.sorted.bed.gz
  bgzip -f $library.I5.sorted.bed
  tabix -f $library.I5.sorted.bed.gz
  bgzip -f $library.I7.sorted.bed
  tabix -f $library.I7.sorted.bed.gz
}

run3
