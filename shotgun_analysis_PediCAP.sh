#!/bin/bash

eval "$(conda shell.bash hook)"

mkdir 02_fastqc_raw 03_fastq_trimmed 04_fastqc_trimmed 05_fastq_unmapped 06_fastqc_unmapped 07_kraken 08_bracken 09_ARGs-OAP

#Assess quality of the reads with FastQC

cd 01_fastq_raw
fastqc *.fastq.gz -o ../02_fastqc_raw/ -t 24
cd ../02_fastqc_raw/
mkdir zips
for d in *.zip; \
	do unzip $d -d zips/; \
done
cd zips/
mkdir txt
for d in *; \
	do cd $d; \
	mv "fastqc_data.txt" ../txt/"$d.txt"; \
	cd ../; \
done
cd txt/
sed -s -n '4p;7p' *_R1_001_fastqc.txt > summary_raw.txt
cd ../../../01_fastq_raw

#Adapter trimming with TrimGalore

for d in *_R1_001.fastq.gz; \
    do trim_galore --output_dir ../03_fastq_trimmed/ --quality 30 --cores 8 --paired $d ${d%%_R1*}_R2_001.fastq.gz; \
    rm $d; \
    rm ${d%%_R1*}_R2_001.fastq.gz; \
done 

#Assess quality of trimmed reads and calculate number of reads removed

cd ../03_fastq_trimmed
fastqc *.fq.gz -o ../04_fastqc_trimmed/ -t 24
cd ../04_fastqc_trimmed/
mkdir zips
for d in *.zip; \
	do unzip $d -d zips/; \
done
cd zips/
mkdir txt
for d in *; \
	do cd $d; \
	mv "fastqc_data.txt" ../txt/"$d.txt"; \
	cd ../; \
done
cd txt/
sed -s -n '4p;7p' *_val_1_fastqc.txt > summary_trimmed.txt
cd ../../../03_fastq_trimmed

#Host reads removal with Hostile

for d in *1.fq.gz; \
    do hostile clean --fastq1 $d --fastq2 ${d%%_R1*}_R2_001_val_2.fq.gz -o ../05_fastq_unmapped -t 24; \
    rm $d; \
    rm ${d%%_R1*}_R2_001_val_2.fq.gz; \
done
cd ../05_fastq_unmapped
for d in *clean_1*; \
    do mv $d ${d%%_S*}_R1.fastq.gz; \
done
for d in *clean_2*; \
    do mv $d ${d%%_S*}_R2.fastq.gz; \
done

#FastQC to calculate number of reads removed in host removal step

fastqc *.f*.gz -o ../06_fastqc_unmapped -t 24
cd ../06_fastqc_unmapped/
mkdir zips
for d in *.zip; \
	do unzip $d -d zips/; \
done
cd zips/
mkdir txt
for d in *; \
	do cd $d; \
	mv "fastqc_data.txt" ../txt/"$d.txt"; \
	cd ../; \
done
cd txt/
sed -s -n '4p;7p' *_R1_fastqc.txt > summary_unmapped.txt
cd ../../../05_fastq_unmapped

#Taxonomical classification of unmapped files

for d in *_R1.fastq.gz; \
    do kraken2 --paired --threads 24 --db ~/reference_databases/uhgg_v2.0.2/ --output ../07_kraken/${d%%_R1*}_uhgg --report ../07_kraken/${d%%_R1*}_uhgg.report $d ${d%%_R1*}_R2.fastq.gz; \
    rm ../07_kraken/${d%%_R1*}_uhgg; \
done
cd ../07_kraken

#Combine kraken results into one file

combine_kreports.py -r *_uhgg.report -o kraken_uhgg_summary.txt

#Run bracken for microbial abundance estimation

for d in *_uhgg.report; \
	do est_abundance.py -i $d -k ~/reference_databases/uhgg_v2.0.2/database150mers.kmer_distrib -o ../08_bracken/${d%%.report}.bracken; \
	   est_abundance.py -i $d -k ~/reference_databases/uhgg_v2.0.2/database150mers.kmer_distrib -l G -o ../08_bracken/${d%%.report}_genus.bracken; \
	   est_abundance.py -i $d -k ~/reference_databases/uhgg_v2.0.2/database150mers.kmer_distrib -l P -o ../08_bracken/${d%%.report}_phylum.bracken; \
done

#Combine bracken output

cd ../08_bracken
combine_bracken_outputs.py --files *_uhgg.bracken -o bracken_uhgg_species.txt
combine_bracken_outputs.py --files *_uhgg_genus.bracken -o bracken_uhgg_genus.txt
combine_bracken_outputs.py --files *_uhgg_phylum.bracken -o bracken_uhgg_phylum.txt

#Run ARGs-OAP for acquired resistance genes detection on the reads

cd ../05_fastq_unmapped
args_oap stage_one -i ./ -o ../09_ARGs-OAP -f fastq.gz -t 24 
cd ../09_ARGs-OAP 
args_oap stage_two -i ./ -t 24
