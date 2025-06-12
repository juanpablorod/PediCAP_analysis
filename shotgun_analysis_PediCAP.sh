#!/bin/bash

eval "$(conda shell.bash hook)"

mkdir 02_fastqc_raw 03_fastq_trimmed 04_fastqc_trimmed 05_fastq_unmapped 06_fastqc_unmapped 07_kraken 08_bracken 09_ARGs-OAP 11_mumame 12_metaSPAdes 13_metaBAT 14_checkm 10_HUMAnN

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

conda activate trimgalore 
for d in *_R1_001.fastq.gz; \
    do trim_galore --output_dir ../03_fastq_trimmed/ --quality 30 --cores 8 --paired $d ${d%%_R1*}_R2_001.fastq.gz; \
    rm $d; \
    rm ${d%%_R1*}_R2_001.fastq.gz; \
done 
conda deactivate

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

conda activate hostile
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
conda deactivate

#Taxonomical classification of unmapped files

conda activate kraken
for d in *_R1.fastq.gz; \
    do kraken2 --paired --threads 24 --db ~/reference_databases/kraken_refseq/ --output ../07_kraken/${d%%_R1*}_refseq --report ../07_kraken/${d%%_R1*}_refseq.report $d ${d%%_R1*}_R2.fastq.gz; \
    rm ../07_kraken/${d%%_R1*}_refseq; \
done
for d in *_R1.fastq.gz; \
    do kraken2 --paired --threads 24 --db ~/reference_databases/uhgg_v2.0.2/ --output ../07_kraken/${d%%_R1*}_uhgg --report ../07_kraken/${d%%_R1*}_uhgg.report $d ${d%%_R1*}_R2.fastq.gz; \
    rm ../07_kraken/${d%%_R1*}_uhgg; \
done
cd ../07_kraken

#Combine kraken results into one file

combine_kreports.py -r *_refseq.report -o kraken_refseq_summary.txt
combine_kreports.py -r *_uhgg.report -o kraken_uhgg_summary.txt

#Run bracken for microbial abundance estimation

for d in *_refseq.report; \
	do est_abundance.py -i $d -k ~/reference_databases/kraken_refseq/database150mers.kmer_distrib -o ../08_bracken/${d%%.report}.bracken; \
	   est_abundance.py -i $d -k ~/reference_databases/kraken_refseq/database150mers.kmer_distrib -l G -o ../08_bracken/${d%%.report}_genus.bracken; \
	   est_abundance.py -i $d -k ~/reference_databases/kraken_refseq/database150mers.kmer_distrib -l P -o ../08_bracken/${d%%.report}_phylum.bracken; \
done
for d in *_uhgg.report; \
	do est_abundance.py -i $d -k ~/reference_databases/uhgg_v2.0.2/database150mers.kmer_distrib -o ../08_bracken/${d%%.report}.bracken; \
	   est_abundance.py -i $d -k ~/reference_databases/uhgg_v2.0.2/database150mers.kmer_distrib -l G -o ../08_bracken/${d%%.report}_genus.bracken; \
	   est_abundance.py -i $d -k ~/reference_databases/uhgg_v2.0.2/database150mers.kmer_distrib -l P -o ../08_bracken/${d%%.report}_phylum.bracken; \
done

#Combine bracken output

cd ../08_bracken
combine_bracken_outputs.py --files *_refseq.bracken -o bracken_refseq_species.txt
combine_bracken_outputs.py --files *_refseq_genus.bracken -o bracken_refseq_genus.txt
combine_bracken_outputs.py --files *_refseq_phylum.bracken -o bracken_refseq_phylum.txt
combine_bracken_outputs.py --files *_uhgg.bracken -o bracken_uhgg_species.txt
combine_bracken_outputs.py --files *_uhgg_genus.bracken -o bracken_uhgg_genus.txt
combine_bracken_outputs.py --files *_uhgg_phylum.bracken -o bracken_uhgg_phylum.txt
conda deactivate

#Run ARGs-OAP for acquired resistance genes detection on the reads

cd ../05_fastq_unmapped
conda activate args_oap
args_oap stage_one -i ./ -o ../09_ARGs-OAP -f fastq.gz -t 24 
cd ../09_ARGs-OAP 
args_oap stage_two -i ./ -t 24

#Run mumame to detect mutations leading to AMR

cd ../05_fastq_unmapped
mumame -d ~/reference_databases/mumame/mutation_database -o ../11_mumame/mumame_output *.fastq.gz

#Assembly of unmapped reads with metaSPAdes

conda activate spades
mkdir ../12_metaSPAdes/scaffolds ../12_metaSPAdes/gfa
for d in *_R1.fastq.gz; \
    do spades.py -1 $d -2 ${d%%_R1*}_R2.fastq.gz -o ../12_metaSPAdes/${d%%_R1*} --meta -t 24; \
    cd ../12_metaSPAdes/${d%%_R1*}; \
    mv "scaffolds.fasta" ../scaffolds/"${d%%_R1*}.fasta"; \
    mv "assembly_graph_with_scaffolds.gfa" ../gfa/"${d%%_R1*}.gfa"; \
    cd ../; \
    rm -rd ${d%%_R1*}; \
    cd ../05_fastq_unmapped; \
done
conda deactivate

#metaBAT for binning

conda activate metabat
cd ../12_metaSPAdes/scaffolds
for d in *.fasta; \
	do metabat -i $d -o ../../13_metaBAT/${d%%.fasta} -t 24; \ 
	done
conda deactivate

#checkM for quality check of bins 

conda activate checkm
cd ../../13_metaBAT 
checkm2 predict --threads 24 --input ./ --output-directory ../14_checkM/ -x fa 
conda deactivate 

#Run HUMAnN for metabolic pathway prediction

cd ../05_fastq_unmapped
mkdir ../10_HUMAnN/output
conda activate humann
for d in *_R1.fastq.gz; \
	do gunzip $d ${d%%_R1*}_R2.fastq.gz; \
	cat ${d%%.gz} ${d%%_R1*}_R2.fastq > ../10_HUMAnN/${d%%_R1*}.fastq; \
	gzip *.fastq; \
	cd ../10_HUMAnN/; \
	humann -i ${d%%_R1*}.fastq -o ./output/${d%%_R1*} --threads 24; \
	rm *.fastq; \
	rm -rd ./output/${d%%_R1*}/*tem*; \
	cd ../05_fastq_unmapped/; \
	done

#Retrieve, normalize and merge all HUMAnN ouptut

cd ../10_HUMAnN/output
mkdir genefamilies pathcoverage pathabundance
for d in *; \
	do humann_renorm_table --input $d/*_genefamilies.tsv --output genefamilies/${d}_genefamilies_relab.tsv --units relab; \
	cp $d/*_pathcoverage.tsv pathcoverage/; \
	humann_renorm_table --input $d/*_pathabundance.tsv --output pathabundance/${d}_pathabundance_relab.tsv --units relab; \
done
humann_join_tables --input genefamilies --output humann_genefamilies.tsv --file_name genefamilies_relab
humann_join_tables --input pathcoverage --output humann_pathcoverage.tsv --file_name pathcoverage
humann_join_tables --input pathabundance --output humann_pathabundance.tsv --file_name pathabundance_relab
conda deactivate

