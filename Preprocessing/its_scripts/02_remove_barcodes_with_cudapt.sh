mkdir trimmed
mkdir outputs

cd trimmed
cp ../inputs/samples .

ln -s ../demuxed/merged/*.gz .

for sample in $(cat samples);
do
echo "On sample: $sample";

/Users/imorelan/miniconda3/envs/qiime2-2019.7/bin/cutadapt  -g GTGCCAGCMGCCGCGGTAA -a ATTAGAWACCCBNGTAGTCC -G GGACTACNVGGGTWTCTAAT -A TTACCGCGGCKGCTGGCAC -m 160 -M 235 -o ${sample}_R1_trimmed.fq.gz -p ${sample}_R2_trimmed.fq.gz ${sample}_R1.fq.gz ${sample}_R2.fq.gz

done
