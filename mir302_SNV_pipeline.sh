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

if [ ! -f "miR302_all_merge.vcf" ]; then
	echo "~/bin/jdk1.8.0_91/bin/java -jar ~/bin/GATK_3.6/GenomeAnalysisTK.jar \\
		-T GenotypeGVCFs \\
		-R ~/hg19_reference_sequence/hg19_chromfa.fa \\
		-o miR302_all_merge.vcf \\" > mergeScript.sh
	for v in */*.vcf
		do
		    echo "-V ${v} \\"
		done;

	chmod 777 mergeScript.sh
	./mergeScript.sh
	rm mergeScript.sh	
	
	if [ ! -f "raw_miR302_SNPs.recal" ]; then
		~/bin/jdk1.8.0_91/bin/java -jar ~/bin/GATK_3.6/GenomeAnalysisTK.jar \
			-T VariantRecalibrator \
			-nt 4 \
			-R ~/hg19_reference_sequence/hg19_chromfa.fa \
			-input miR302_all_merge.vcf \
					##I am assuming that I use the merged variant 
					##calls and not the variant calls for the individual samples
			-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp_137_b37.vcf \
			-resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.b37.sites.vcf \
			-resource:omni,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.b37.sites.vcf \
			-resource:1000G,known=false,training=true,truth=false,prio=10.0 1000G_phase1.snps.high_confidence.vcf \
					##Double check tool documentation not sure about the vcf file. 
					##There tutorial says it is included in the distro...
			-an DP -an QD -an FS -an MQRankSum -an ReadPosRankSum -an InbreedingCoeff -an SOR -MQ \
			-mode SNP \
			-recalFile raw_miR302_SNPs.recal \
			-tranchesFile raw_miR302_SNPs.tranches \
			-rscriptFile snp.recal.plots.R
	fi

	if [ ! -f "raw_miR302_SNPs.recal" ]; then
		~/bin/jdk1.8.0_91/bin/java -jar ~/bin/GATK_3.6/GenomeAnalysisTK.jar \
			-T ApplyRecalibration \
			-R ~/hg19_reference_sequence/hg19_chromfa.fa \
			-input miR302_all_merge.vcf \
			-mode SNP \
			-recalFile raw_miR302_SNPs.recal \
			-tranchesFile raw_miR302_SNPs.tranches \
			-o miR302_recalSNPs.vcf \
			-ts_filter_level 95.0
	fi

	echo"
		############################################################################
		##################### Calibrate Indel call scorese #########################
		############################################################################
		"

	if [ ! -f "miR302_recalSNPs_Indel.recal" ]; then
		~/bin/jdk1.8.0_91/bin/java -jar ~/bin/GATK_3.6/GenomeAnalysisTK.jar \
			-T VariantRecalibrator \
			-nt 4 \
			-R ~/hg19_reference_sequence/hg19_chromfa.fa \
			-input miR302_recalSNPs.vcf \
					##I am assuming that I use the merged variant 
					##calls and not the variant calls for the individual samples
			-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp.b37.vcf \
			-resource:mills,known=false,training=true,truth=true,prior=12.0 Mills_and_1000G_gold_standard.indels.b37.sites.vcf \
					##Double check tool documentation not sure about the vcf file. 
					##There tutorial says it is included in the distro...
			-an DP -an QD -an FS -an MQRankSum -an ReadPosRankSum -an InbreedingCoeff \
			-mode INDEL \
			--maxGaussians 4 \
			-recalFile miR302_recalSNPs_Indel.recal \
			-tranchesFile miR302_recalSNPs_Indel.tranches \
			-rscriptFile indel.recal.plots.R
	fi

	if [ ! -f "raw_miR302_SNPs.recal" ]; then
		~/bin/jdk1.8.0_91/bin/java -jar ~/bin/GATK_3.6/GenomeAnalysisTK.jar \
			-T ApplyRecalibration \
			-R ~/hg19_reference_sequence/hg19_chromfa.fa \
			-input miR302_recalSNPs.vcf \
			-mode INDEL \
			-recalFile miR302_recalSNPs_Indel.recal \
			-tranchesFile miR302_recalSNPs_Indel.tranches \
			-o miR302_recalSNPs_recalIndels.vcf \
			-ts_filter_level 95.0
	fi