#! /bin/bash

library=$1
root_dir=$2

## MACS2 and SEACR peak calling
ATAC_SIZE=$3
QVAL=$4
## SEACR peak calling
CTAG_SIZE=$5
SIZE=$6

I7=$root_dir/pool/$library.I7.sorted.bed.gz
atac_dir=$root_dir/pseudo_ATAC/MACS2/$library/
I5=$root_dir/pool/$library.I5.sorted.bed.gz
ctag_dir=$root_dir/pseudo_CUT/SEACR/$library/

MACS=/home/linzhuobin/anaconda3/envs/scNanoATAC-env/bin/macs2
SEACR=/home/linzhuobin/software/SEACR/SEACR_1.3.sh
BT2_bin=/home/linzhuobin/anaconda3/envs/scNanoATAC-env/bin

NPEAK=300000

eval "$(conda shell.bash hook)"
conda activate lh3

run4(){
  hg38

  MACS_default
  SEACR_call

}

# human reference
hg38(){
    ref_genome_fa=/home/linzhuobin/reference/human/GRCh38_v44/GRCh38.primary_assembly.genome.fa
    ref_prefix=/home/linzhuobin/reference/human/GRCh38_v44/GRCh38.bowtie2
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

# run by library
MACS_default(){
    init_dir $atac_dir
    mkdir -p MACS_default

    echo "Peak calling - MACS (default/narrow)"
    ## Peak calling - MACS
    EXT_SIZE=$ATAC_SIZE
    NAME=$atac_dir/MACS_default/all
    $MACS callpeak -t $I7 \
    -f BED -n $NAME -g ${GENOMESIZE} \
    -q $QVAL --nomodel --shift -$(($EXT_SIZE/2)) --extsize $EXT_SIZE --keep-dup all -B --SPMR 2> $atac_dir/MACS_default/run.log

    echo "Output MACS results..."
    ## Filter out peaks and rename
    atac_raw=$atac_dir/MACS_default/all_peaks.narrowPeak
    atac_flt=$atac_dir/MACS_default/chr_peaks.narrowPeak
    sort -k 8gr,8gr $atac_raw | grep -P 'chr[\dXY]+[ \t]'| awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}' | \
    head -n ${NPEAK} | bedtools sort > $atac_flt
    gzip -c $atac_flt > $atac_flt.gz
    echo "MACS done."

    EXT_SIZE1=$(($ATAC_SIZE/2))
    zcat $I7 |awk -vOFS='\t' '{if($6=="+") {print $1,$2-"'"$EXT_SIZE1"'",$2+"'"$EXT_SIZE1"'","N","1000",$6,$7} else {print $1,$3-"'"$EXT_SIZE1"'",$3+"'"$EXT_SIZE1"'","N","1000",$6,$7}}' |tr ' ' '\t' | bedtools sort | bgzip -f > $library.I7.sorted.bed.gz
    I7=$atac_dir/$library.I7.sorted.bed.gz
}

SEACR_call(){
    init_dir $ctag_dir
    mkdir -p stringent

    EXT_SIZE=$CTAG_SIZE
    zcat $I5 |awk -vOFS='\t' '{if($6=="+") {print $1,$2,$2+"'"$EXT_SIZE"'","N","1000",$6,$7} else {print $1,$3-"'"$EXT_SIZE"'",$3,"N","1000",$6,$7}}' |tr ' ' '\t' | bedtools sort | bgzip -f > $library.I5.sorted.bed.gz
    I5=$ctag_dir/$library.I5.sorted.bed.gz

    cd stringent
    echo "Peak calling - SEACR"
    ## Peak calling - SEACR
    bedtools genomecov -bg -i ../$library.I5.sorted.bed.gz -g $CHRSIZEFILE > $library.bedgraph
    bash $SEACR $library.bedgraph $SIZE non stringent all_SEACR.peak 2> $ctag_dir/stringent/run.log

    echo "Output SEACR results..."
    ## Filter out peaks and rename
    ctag_raw=$ctag_dir/stringent/all_SEACR.peak.stringent.bed
    ctag_flt=$ctag_dir/stringent/chr_SEACR.peak.stringent.bed
    sort -k 4gr,4gr $ctag_raw | grep -P 'chr[\dXY]+[ \t]'| awk 'BEGIN{OFS="\t"}{$7="Peak_"NR ; print $0}' | \
    head -n ${NPEAK} | bedtools sort > $ctag_flt
    gzip -c $ctag_flt > $ctag_flt.gz
    echo "SEACR done."
}

run4
