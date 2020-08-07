cat inputs/its_primer_barcodes_sabre.tsv | cut -f2 | cut -f1 -d "_" > inputs/samples


mkdir demuxed
cd demuxed
mkdir round_one
cd round_one

~/sabre-master/sabre pe -f ../../raw_fastq/KSBITS_S1_L001_R1_001.fastq.gz -r ../../raw_fastq/KSBITS_S1_L001_R2_001.fastq.gz -b ../../inputs/its_primer_barcodes_sabre.tsv -u no_bc_match_R1.fq -w no_bc_match_R2.fq

cd ..
mkdir round_two
cd round_two

~/sabre-master/sabre pe -f ../round_one/no_bc_match_R2.fq -r ../round_one/no_bc_match_R1.fq -b ../../inputs/its_primer_barcodes_sabre.tsv -u round2_no_bc_match_R1.fq -w round2_no_bc_match_R2.fq

cd .. 
mkdir merged

for sample in $(cat ../inputs/samples);
do
echo "On sample: $sample";

cat round_one/${sample}_R1.fq round_two/${sample}_R1.fq > merged/${sample}_R1.fq 
cat round_one/${sample}_R2.fq round_two/${sample}_R2.fq > merged/${sample}_R2.fq 

gzip merged/${sample}_R1.fq
gzip merged/${sample}_R2.fq

done
