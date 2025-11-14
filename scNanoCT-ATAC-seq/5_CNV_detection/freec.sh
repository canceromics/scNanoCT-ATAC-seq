#! /bin/bash

library=$1
sample=$2
root_dir=$3

freec_dir=$root_dir/freec/$library
cell_list=$root_dir/barcode/${library}/cell_list

eval "$(conda shell.bash hook)"
conda activate lh3
PathToFREEC=/hwdata/home/linzhuobin/software/FREEC-11.6b

run(){
  hg38

  freec
}

# human reference
hg38(){
    ref_genome_fa=/home/linzhuobin/reference/human/GRCh38_v44/GRCh38.primary_assembly.genome.fa
    ref_prefix=/home/linzhuobin/reference/human/GRCh38_v44/GRCh38.bowtie2
    CHRLIST=/home/linzhuobin/reference/human/GRCh38_v44/chr.list
    ## for peak calling
    GENOMESIZE=hs
    CHRSIZEFILE=/home/linzhuobin/reference/human/GRCh38_v44/hg38.chrom.sizes

    if [ -s ${ref_prefix}.1.bt2 ];then
        echo "Found bowtie2 indexed reference genome."
    else
        echo "Creating bowtie2 indexed reference genome..."
        $BT2_bin/bowtie2-build -f $ref_genome_fa --threads $threads GRCh38.bowtie2
    fi
}

init_dir(){
  mkdir -p $1/log/
  cd $1
}

freec(){
    init_dir $freec_dir
    
    echo $root_dir $library $sample
    if [ -s $library.$sample.bam_sample.cpn ];then rm $library.$sample.bam_sample.cpn ;fi
    samtools view -H $root_dir/alignment/$library/$library.$sample.mapQ30.rmdup.sorted.bam | grep -v "SN:GL\|SN:KI" > $library.$sample.sam
    samtools view -@4 $root_dir/alignment/$library/$library.$sample.mapQ30.rmdup.sorted.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
	   chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM >> $library.$sample.sam
    samtools view -b $library.$sample.sam -o $library.$sample.bam && rm $library.$sample.sam
    samtools index $library.$sample.bam
    cat $PathToFREEC/single_cell.conf | sed -e "s|root_dir|$root_dir|g" -e "s|library|$library|g" -e "s|input.bam|$library.$sample.bam|g" > "$sample.conf"
    #cat $sample.conf
    $PathToFREEC/src/freec -conf $sample.conf > log/$library.$sample.log 2>&1
}

run
