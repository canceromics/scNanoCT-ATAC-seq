#!/bin/bash
input=$1
root_dir=$2
# input contents:
# column1: start outer_barcode
# column2: end outer_barcode
# column3: start inner_barcode
# column4: end inner_barcode
# column5: cell name prefix (optional): cell line/date/replicate, letters only

#root_dir=/home/linzhuobin/scNanoCTA/basic
IBC=$root_dir/inner_barcode.txt
OBC=$root_dir/outer_barcode.txt
seq1=CAAGCAGAAGACGGCATACGAGAT
#seq2=AGATGTGTATAAGAGACAG
#seqI5=TCGTCGGCAGCGTC
#seqI7=GTCTCGTGGGCTCGG
FI5=TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
FI7=GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
RI5=CTGTCTCTTATACACATCTGACGCTGCCGACGA
RI7=CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
# Structure: xxx[outer_barcode][seq1][inner_barcode][seqI5/seqI7][seq2]xxx

# I5/I7 dual barcode pair
if [ -s index.fa ];then
        echo "Existed index.fa was removed!"
        rm index.fa
fi
echo ">FI5" >> index.fa
echo "${FI5}" >> index.fa
echo ">FI7" >> index.fa
echo "${FI7}" >> index.fa
echo ">RI5" >> index.fa
echo "${RI5}" >> index.fa
echo ">RI7" >> index.fa
echo "${RI7}" >> index.fa

if [ -s Findex.fa ];then
        echo "Existed Findex.fa was removed!"
        rm Findex.fa
fi
echo ">I5" >> Findex.fa
echo "${FI5}" >> Findex.fa
echo ">I7" >> Findex.fa
echo "${FI7}" >> Findex.fa

if [ -s Rindex.fa ];then
        echo "Existed Rindex.fa was removed!"
        rm Rindex.fa
fi
echo ">I5" >> Rindex.fa
echo "${RI5}" >> Rindex.fa
echo ">I7" >> Rindex.fa
echo "${RI7}" >> Rindex.fa

# cell barcode
if [ -s barcode.fa ];then
	echo "Existed barcode.fa was removed!"
	rm barcode.fa
fi

if [ -s outer_barcode.fa ];then
        echo "Existed outer_barcode.fa was removed!"
        rm outer_barcode.fa
fi

if [ -s inner_barcode.fa ];then
        echo "Existed inner_barcode.fa was removed!"
        rm inner_barcode.fa
fi

while read i
do
    a=($i)
    if [ ${a[0]} -ge 0 ];then
        OBC_start=${a[0]}
        OBC_end=${a[1]}
        IBC_start=${a[2]}
        IBC_end=${a[3]}
        CName=${a[4]}
        [ $OBC_start ] || OBC_start=1
        [ $OBC_end ] || OBC_end=8
        [ $IBC_start ] || IBC_start=1
        [ $IBC_end ] || IBC_end=12
        [ $CName ] || CName=N
        
        for oi in `seq ${OBC_start} ${OBC_end}`
        do
            OBC_seq=`cat $OBC | awk -v i=${oi} '$1==i {print $2}'`
            echo ">${oi}"   >> outer_barcode.fa
            echo "$OBC_seq" >> outer_barcode.fa
            for ii in `seq ${IBC_start} ${IBC_end}`
            do
            IBC_seq=`cat $IBC | awk -v i=${ii} '$1==i {print $2}'`
            Name=${oi}"_"${ii}
            BARCODE=${OBC_seq}${seq1}${IBC_seq}
            #echo ${Name}":"${BARCODE}
            echo ">${Name}" >> barcode.fa
            echo "$BARCODE" >> barcode.fa
            done
        done

        for ii in `seq ${IBC_start} ${IBC_end}`
        do
        IBC_seq=`cat $IBC | awk -v i=${ii} '$1==i {print $2}'`
        echo ">${ii}"  >> inner_barcode.fa
        echo "$IBC_seq" >> inner_barcode.fa
        done
    else
        echo "Barcode is no integer!"
        exit 1
    fi
done < $input
