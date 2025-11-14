#! /bin/bash

library=$1
threads=$2
root_dir=$3
bc_index=$4
his_cell=$5

FILT_SIZE=100
#FILT_SIZE=1000
export PATH=/home/linzhuobin/software:$PATH
export PATH=/home/linzhuobin/software/nanoplexer:$PATH
eval "$(conda shell.bash hook)"
conda activate scNanoATAC-env

raw_dir=$root_dir/raw_data/$library/
barcode_dir=$root_dir/barcode/$library/
trim_dir=$root_dir/trim/$library/

dual_bc_script=$root_dir/barcode.sh
barcode=$barcode_dir/barcode.fa
outcode=$barcode_dir/outer_barcode.fa
inncode=$barcode_dir/inner_barcode.fa
cell_list=$barcode_dir/cell_list
index_list=$barcode_dir/index_list

run1(){
  load_dual

  outer_barcode
  inner_barcode

  dual_index
}

init_dir(){
  mkdir -p $1/log/
  cd $1
}

load_dual(){
  ## output barcode, index, cell_list
  init_dir $barcode_dir
  bash $dual_bc_script $bc_index $root_dir
  seqkit seq -i -n barcode.fa > cell_list
}

dual_barcode(){
  init_dir $barcode_dir

  nanoplexer \
  -b $barcode \
  -t $threads \
  -p .\
  $root_dir/raw_data/$library/*.pass.f*q.gz 2> $barcode_dir/log/nanoplexer.log
  #ls [0-9]*.fastq |sed 's/.fastq//g' > $cell_list
}

## all histone library
outer_barcode(){
  init_dir $barcode_dir

  nanoplexer \
  -b $outcode \
  -t $threads \
  -p .\
  $root_dir/raw_data/$library/*.pass.f*q.gz 2> $barcode_dir/log/nanoplexer_outer.log

  for his in `cut -f1 $his_cell|uniq`
  do
      echo $his
      his_dir=$root_dir/barcode/$his
      init_dir $his_dir
      cd $barcode_dir

      #for oi in `grep -w $his $his_cell| cut -f2`
      for oi in `awk '{if($1=="'"$his"'") {print $2}}' $his_cell`
      do
          if [ -s ${oi}.fastq ];then
              echo $oi
              mv ${oi}.fastq $his_dir
          fi
      done
  done
}

inner_barcode(){
  for his in `cut -f1 $his_cell|uniq`
  do
      echo $his
      his_dir=$root_dir/barcode/$his
      cd $his_dir
      echo $his > $his_dir/log/nanoplexer_inner.log
      
      #for oi in `grep -w $his $his_cell| cut -f2`
      for oi in `awk '{if($1=="'"$his"'") {print $2}}' $his_cell`
      do
          echo  $oi >> $his_dir/log/nanoplexer_inner.log
          mkdir $oi && cd $oi
          nanoplexer -b $inncode -t $threads -p ./ ../${oi}.fastq 2>> $his_dir/log/nanoplexer_inner.log
          for fq in `ls *.fastq|grep -v "unclassified"`;do mv $fq ../${oi}_$fq;done
          cat unclassified.fastq >> ../unclassified.fastq
          cd ../
          rm -r ${oi} && rm ${oi}.fastq
      done
      ls *.fastq|grep -v "unclassified"|sed 's/.fastq//g' > cell_list
      cp $barcode_dir/*index.fa .
  done
}

## single histone library
dual_index1(){
  ## ATAC(I7)/CUT&TAG(I5) index
  library=$1
  trim_dir=$root_dir/trim/$library/
  barcode_dir=$root_dir/barcode/$library/
  index_list=$barcode_dir/index_list
  cell_list=$barcode_dir/cell_list

  init_dir $trim_dir
  cd $barcode_dir

  if [ -s $index_list ];then rm $index_list;fi
  for bc in `cat $cell_list`
  do
    if [ -d $bc ]
    then
        rm -r $bc
    else
        mkdir -p $bc && cd $bc
        cp $barcode_dir/*index.fa .
        ## 5' index
        cutadapt -e 0.2 -O 12 -g file:Findex.fa -o {name}.fastq ../${bc}.fastq \
        -j $threads > ${trim_dir}/log/${library}.${bc}.txt 2>&1
        lid_list=`ls *.fastq | sed 's/.fastq//g'`
        ## 3' index
        for lid in ${lid_list[@]}
        do
            cutadapt -e 0.2 -O 12 -a file:Rindex.fa -o ${lid}_{name}.fastq ${lid}.fastq \
            --minimum-length $FILT_SIZE -j $threads >> ${trim_dir}/log/${library}.${bc}.txt 2>&1
            rm ${lid}.fastq
        done
        ## read_id barcode index
        rename 'unknown' 'UN' * && rename 'unknown' 'UN' *
        id_list=`ls *.fastq | sed 's/.fastq//g'`
	for id in ${id_list[@]}
        do
            seqkit seq -i -n -j $threads ${id}.fastq| \
            awk '{print "'$bc'","'$id'",$1}'| tr ' '  '\t' >> $index_list
        done
        cat *.fastq | gzip > ${trim_dir}/${library}.${bc}.all_trimed.fq.gz
        rm UN_UN.fastq # reads without index
        cat *.fastq | gzip > ${trim_dir}/${library}.${bc}.trimed.fq.gz
        cd $barcode_dir && rm -r $bc
    fi
  done
}

dual_index(){
  for his in `cut -f1 $his_cell|uniq`
  do
      dual_index1 $his
  done
}


run1
