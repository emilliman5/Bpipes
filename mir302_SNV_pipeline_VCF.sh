for fa in *.bam;

do
	f=${fa%_bowtie2*}
	echo "##############################################################
		############## Build Base Recalibration Table ################
		##############################################################"
	if [ ! -f "${f}_recal.table" ]; then
		~/bin/jdk1.8.0_91/bin/java -jar ~/bin/GATK_3.6/GenomeAnalysisTK.jar \
			-T BaseRecalibrator \
			-nct 8 \
			-R ~/hg19_reference_sequence/hg19_chromfa.fa \
			-I ${f}_bowtie2_N1_hg19_merged_sort_dedupe_realigned.bam \
			-knownSites ~/hg19_featureSets/dbsnp147_GRCh37p13_all_2.vcf \
			-o ${f}_recal.table \
			-L ../miR302.intervals
## If not doing WGS the sequencing targets need to be specified to 
## prevent the recalibration from going awry

	fi
	echo "##############################################################
		################# Recalibrate sequence data ##################
		##############################################################"
	if [ ! -f "${f}_bowtie2_N1_hg19_merged_sort_dedupe_realigned_recal.bam" ]; then
		~/bin/jdk1.8.0_91/bin/java -jar ~/bin/GATK_3.6/GenomeAnalysisTK.jar \
			-T PrintReads \
			-nct 8 \
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











