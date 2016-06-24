for fa in *1.sanger.fastq;

do
	mkdir ${fa%%.1.sanger.fastq}
	cd ${fa%%.1.sanger.fastq}/
	f=${fa%%.1.sanger.fastq}

	cutadapt -a CTGTCTCTTATA -A CTGTCTCTTATA -o ${f}.1.sanger_cutAdapt.fastq \
		-p ${f}.2.sanger_cutdapt.fastq ../$f*1*.fastq ../$f*2*.fastq

	java -Xmx16g -XX:ParallelGCThreads=4 -jar ${PICARD_DIR}/FastqToSam.jar \
		FASTQ=../${f}.1.sanger.fastq \
		FASTQ2=../${f}.2.sanger.fastq \
		OUTPUT=${f}_unmapped.sanger.fastq
		SAMPLE_NAME=${f}
		SORT_ORDER=queryname

	bowtie2 -x hg19 -1 ${f}.1.sanger.cutAdapt.fastq -2 ${f}.2.sanger.cutAdapt.fastq \
		-X 2000 -N 1 -p 4 --dovetail --un-conc ${f}_discordant.fastq -S ${f}_bowtie2_N1_hg19.sam

	java -Xmx16g -XX:ParallelGCThreads=4 -jar ${PICARD_DIR}/MergeBamAlignment.jar \
		ALIGNED=${f}_bowtie2_N1_hg19.sam
		UNMAPPED=${f}_unmapped.sanger.fastq
		OUTPUT=${f}_bowtie2_N1_hg19_merged.bam
		REFERENCE_SEQUENCE=~/hg19_reference_sequence/

	java -Xmx16g -XX:ParallelGCThreads=4 -jar ${PICARD_DIR}/SortSam.jar \
		I=${f}_bowtie2_N1_hg19_merged.bam
		O=${f}_bowtie2_N1_hg19_merged_sort.bam
		SORT_ORDER=coordinate 

	java -Xmx16g -XX:ParallelGCThreads=4 -jar ${PICARD_DIR}/MarkDuplicates.jar \
        		MAX_RECORDS_IN_RAM=5000000 \
        		TMP_DIR=/ddn/gs1/home/millimanej/tmp/ \ 
        		VERBOSITY=DEBUG \ 
				INPUT=${f}_bowtie2_N1_hg19_merged_sort.bam \ 
				OUTPUT=${f}_bowtie2_N1_hg19_merged_sort_dedupe.bam \ 
				METRICS_FILE=${f}_bowtie2_N1_hg19_merged_sort_dedupe.metrics.txt \ 
				REMOVE_DUPLICATES=true \
				ASSUME_SORTED=true \ 

	~/bin/jdk1.8.0_91/bin/java -jar ~/bin/GATK3.6/bin/GenomeAnalysisTK.jar -T RealignerTargetCreator \
		-R ~/hg19_reference_sequence/hg19_chromfa.fa
		-I ${f}_bowtie2_N1_hg19_merged_sort_dedupe.bam
		-o ${f}_forIndelRealigner.intervals
		-nt 4

	~/bin/jdk1.8.0_91/bin/java -jar ~/bin/GATK3.6/bin/GenomeAnalysisTK.jar -T IndelRealigner \
		-R ~/hg19_reference_sequence/hg19_chromfa.fa
		-I ${f}_bowtie2_N1_hg19_merged_sort_dedupe.bam
		-targetIntervals ${f}_forIndelRealigner.intervals
		-o ${f}_bowtie2_N1_hg19_merged_sort_dedupe_realigned.bam
	cd ..
done;
	