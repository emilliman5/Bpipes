##Script to perform variant calling analyses using MIT's GATK v3.6
##Written by: Eric Milliman
##Updated: 06/28/2016

##Note: Indel realignment is capped at a depth of 200000. This should be fine for the miR302 analysis but a down-sampling
##may need to be performed just in case.

for fa in *1.sanger.fastq;

do
	if [ ! -d "${fa%%.1.sanger.fastq}" ]; then
		mkdir ${fa%%.1.sanger.fastq}
	fi

	cd ${fa%%.1.sanger.fastq}/
	f=${fa%%.1.sanger.fastq}
	echo $f

	echo "
		######################################################		
		################Trim adaptor sequences################
		######################################################
		"
	if [ ! -f "${f}.1.sanger_cutAdapt.fastq" ] || [ ! -f "${f}.2.sanger_cutAdapt.fastq" ] \
	cutadapt -a CTGTCTCTTATA -A CTGTCTCTTATA -o ${f}.1.sanger_cutAdapt.fastq \
		-p ${f}.2.sanger_cutAdapt.fastq ../$f*1*.fastq ../$f*2*.fastq
	fi

	echo "
		#################################################################################
		################Merge unmapped Fastq files into unmapped Sam file################
		#################################################################################
		"
	if [ ! -f "*_unmapped.sanger.sanger_cutAdapt.fastq"] || [ ! -f "${f}.1.sanger_cutAdapt.fastq" ] || \
		[ ! -f "${f}.2.sanger_cutAdapt.fastq" ]; then
		java -Xmx16g -XX:ParallelGCThreads=4 -jar ${PICARD_DIR}/FastqToSam.jar \
			FASTQ=../${f}.1.sanger.fastq \
			FASTQ2=../${f}.2.sanger.fastq \
			OUTPUT=${f}_unmapped.sanger.fastq \
			SAMPLE_NAME=${f} \
			SORT_ORDER=queryname
	fi

	echo "
		###############################################
		################Align sequences################
		###############################################
		"
	if [ ! -f "*._bowtie2_N1_hg19.sam"] || [ ! -f "*_unmapped.sanger.sanger_cutAdapt.fastq"]; then
		bowtie2 -x hg19 -1 ${f}.1.sanger.cutAdapt.fastq -2 ${f}.2.sanger.cutAdapt.fastq \
			-X 2000 -N 1 -p 4 --dovetail --un-conc ${f}_discordant.fastq -S ${f}_bowtie2_N1_hg19.sam
	fi

	echo "
		###############################################################################################
		################Merge aligned reads with unaligned SAM file to rebuild metadata################
		###############################################################################################
		"
	if [ ! -f "*.merged.bam"] || [ ! -f "*._bowtie2_N1_hg19.sam"]; then
		java -Xmx16g -XX:ParallelGCThreads=4 -jar ${PICARD_DIR}/MergeBamAlignment.jar \
			ALIGNED=${f}_bowtie2_N1_hg19.sam \
			UNMAPPED=${f}_unmapped.sanger.fastq \
			OUTPUT=${f}_bowtie2_N1_hg19_merged.bam \
			REFERENCE_SEQUENCE=~/hg19_reference_sequence/
	fi

	echo "
		################################################################
		################Sort aligned reads by coordinate################
		################################################################
		"
	if [ ! -f "*sort.bam"] || [ ! -f "*.merged.bam"]; then
		java -Xmx16g -XX:ParallelGCThreads=4 -jar ${PICARD_DIR}/SortSam.jar \
			I=${f}_bowtie2_N1_hg19_merged.bam \
			O=${f}_bowtie2_N1_hg19_merged_sort.bam \
			SORT_ORDER=coordinate
	fi

	echo "
		###################################################################
		################Mark and Remove duplicate sequences################
		###################################################################
		"
	if [ ! -f "*dedupe.bam"] || [ ! -f "*sort.bam"]; then 
		java -Xmx16g -XX:ParallelGCThreads=4 -jar ${PICARD_DIR}/MarkDuplicates.jar \
		    		MAX_RECORDS_IN_RAM=5000000 \
		    		TMP_DIR=/ddn/gs1/home/millimanej/tmp/ \
		    		VERBOSITY=DEBUG \
					INPUT=${f}_bowtie2_N1_hg19_merged_sort.bam \
					OUTPUT=${f}_bowtie2_N1_hg19_merged_sort_dedupe.bam \
					METRICS_FILE=${f}_bowtie2_N1_hg19_merged_sort_dedupe.metrics.txt \
					REMOVE_DUPLICATES=true \
					ASSUME_SORTED=true 
	fi

		echo "
		###################################################################
		###################### Bam file indexing ##########################
		###################################################################
		"
	if [ ! -f "${f}_bowtie2_N1_hg19_merged_sort_dedupe.bam.bai" ]; then 
		samtools index ${f}_bowtie2_N1_hg19_merged_sort_dedupe.bam
	fi

	echo "###############################################################################################
		  ################Create intervals at which local realignment is to be performed ################
		  ###############################################################################################"
	if [ ! -f "*_forIndelRealigner.intervals" ] || [ ! -f "*dedupe.bam"]; then
		~/bin/jdk1.8.0_91/bin/java -jar ~/bin/GATK3.6/bin/GenomeAnalysisTK.jar -T RealignerTargetCreator \
			-R ~/hg19_reference_sequence/hg19_chromfa.fa \
			-I ${f}_bowtie2_N1_hg19_merged_sort_dedupe.bam \
			-o ${f}_forIndelRealigner.intervals \
			-nt 8 \
			-maxReadsForRealignment 200000 \
			-maxReadsInMemory 200000
	fi	

	echo "#########################################################
	      ################Perform local realignment################
	      #########################################################"
	if [ ! -f "*realigned.bam" ] || [ ! -f "*_forIndelRealigner.intervals" ]; then
		~/bin/jdk1.8.0_91/bin/java -Xmx16g -jar ~/bin/GATK3.6/bin/GenomeAnalysisTK.jar -T IndelRealigner \
			-R ~/hg19_reference_sequence/hg19_chromfa.fa
			-I ${f}_bowtie2_N1_hg19_merged_sort_dedupe.bam
			-targetIntervals ${f}_forIndelRealigner.intervals
			-o ${f}_bowtie2_N1_hg19_merged_sort_dedupe_realigned.bam
	fi

	echo "########################################################
		############## Index Realigned Bam File ################
		########################################################"
	if [ ! -f "${f}_bowtie2_N1_hg19_merged_sort_dedupe_realigned.bam.bai" ]; then 
		samtools index ${f}_bowtie2_N1_hg19_merged_sort_dedupe_realigned.bam
	fi

	echo "##############################################################
		############## Build Base Recalibration Table ################
		##############################################################"
	if [ ! -f "${f}_recal.table" ]; then
		~/bin/jdk1.8.0_91/bin/java -jar ~/bin/GATK_3.6/GenomeAnalysisTK.jar \
			-T BaseRecalibrator \
			-R ~/hg19_reference_sequence/hg19_chromfa.fa \
			-I ${f}_bowtie2_N1_hg19_merged_sort_dedupe_realigned.bam \
			-knownSites dbsnp147_GRCh37p13_common_all.vcf \
			-o ${f}_recal.table \
			-L ../miR302.intervals	## If not doing WGS the sequencing targets need to be specified to 
						## prevent the recalibration from going awry
	fi
	echo "##############################################################
		################# Recalibrate sequence data ##################
		##############################################################"
	if [ ! -f "${f}_bowtie2_N1_hg19_merged_sort_dedupe_realigned_recal.bam" ]; then
		~/bin/jdk1.8.0_91/bin/java -jar ~/bin/GATK_3.6/GenomeAnalysisTK.jar \
			-T PrintReads \
			-R ~/hg19_reference_sequence/hg19_chromfa.fa \
			-I ${f}_bowtie2_N1_hg19_merged_sort_dedupe_realigned.bam \
			-BQSR ${f}_recal.table \
			-o ${f}_bowtie2_N1_hg19_merged_sort_dedupe_realigned_recal.bam
	fi

	echo "##################################################################
		############## Build 2nd Base Recalibration Table ################
		##################################################################"
	if [ ! -f "${f}_after_recal.table" ]; then
		~/bin/jdk1.8.0_91/bin/java -jar ~/bin/GATK_3.6/GenomeAnalysisTK.jar \
			-T BaseRecalibrator \
			-R ~/hg19_reference_sequence/hg19_chromfa.fa \
			-I ${f}_bowtie2_N1_hg19_merged_sort_dedupe_realigned.bam \
			-knownSites dbsnp147_GRCh37p13_common_all.vcf \
			-BQSR ${f}_recal.table \
			-o ${f}_after_recal.table \
			-L ../miR302.intervals	## If not doing WGS the sequencing targets need to be specified to 
						## prevent the recalibration from going awry
	fi
	
	echo "##################################################################
		#################### Recalibration QC plots ######################
		##################################################################"

	if [ ! -f "${f}_recal_plots.pdf" ]; then
		~/bin/jdk1.8.0_91/bin/java -jar ~/bin/GATK_3.6/GenomeAnalysisTK.jar \
			-T AnalyzeCovariates \
			-R ~/hg19_reference_sequence/hg19_chromfa.fa \
			-before ${f}_recal.table \
			-after ${f}_after_recal.table \
			-plots ${f}_recal_plots.pdf
	fi

	if [ ! -f "${f}_recal_plots.pdf" ]; then
	~/bin/jdk1.8.0_91/bin/java -jar ~/bin/GATK_3.6/GenomeAnalysisTK.jar \
		-T HaplotypeCaller \
		-R ~/hg19_reference_sequence/hg19_chromfa.fa \
		-I ${f}_bowtie2_N1_hg19_merged_sort_dedupe_realigned_recal.bam
		-o ${f}_bowtie2_N1_hg19_merged_sort_dedupe_realigned_recal.g.vcf
		-L ../miR302.intervals
		-ERC GVCF
	fi

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
fi
	
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
				##Their tutorial says it is included in the distro...
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