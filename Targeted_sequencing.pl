#!/usr/bin/perl
use Data::Dumper;

######################################################################################################
##                         Genome Variant Calling - From Fastq to VCF Pipeline
##                Applied to Target Sequencing - low sample size and low sequences size
##     Adapted from the Best Practices for Germline short variant discovery of the Broad Institute
##              https://software.broadinstitute.org/gatk/best-practices/workflow?id=11145
##
######################################################################################################

## Authors: Wagner Magalhaes, Rennan Moreira and Zsolt Balázs
## Date: 02-03-2017 - 02-04-2018
## Dependencies:  FastQC, Trimmomatic, Samtools, GATK, Picard,
## Requirements in the working directory:
## A FASTA file listing the adapters that were used in the experiments: haloplexadapters.fa
## A list of the paths to where the .vcfs of the joint genotype calling are expected: intervalmap.list
## A list of the samples and the path to where their .vcfs files are expected: samplemap_all.tsv
## FASTQ files (can be gzipped) containing the sequenced reads Sequences/$data.fq
## A list of the targeted intervals: Target3.interval.list
##
## The comments of each step that diverge from the Best Practices mention it.

system "mkdir /home/zsolt/Desktop/Target_pipeline/newqual/step01 | mkdir /home/zsolt/Desktop/Target_pipeline/newqual/step02 | mkdir /home/zsolt/Desktop/Target_pipeline/newqual/step03 | mkdir /home/zsolt/Desktop/Target_pipeline/newqual/step04 | mkdir /home/zsolt/Desktop/Target_pipeline/newqual/step05 | mkdir /home/zsolt/Desktop/Target_pipeline/newqual/step06 | mkdir /home/zsolt/Desktop/Target_pipeline/newqual/step07 | mkdir /home/zsolt/Desktop/Target_pipeline/newqual/step08 | mkdir /home/zsolt/Desktop/Target_pipeline/newqual/step09";
system "mkdir /home/zsolt/Desktop/Target_pipeline/newqual/step10 | mkdir /home/zsolt/Desktop/Target_pipeline/newqual/step11 | mkdir /home/zsolt/Desktop/Target_pipeline/newqual/step12 | mkdir /home/zsolt/Desktop/Target_pipeline/newqual/step13 | mkdir /home/zsolt/Desktop/Target_pipeline/newqual/step14 | mkdir /home/zsolt/Desktop/Target_pipeline/newqual/step15 | mkdir /home/zsolt/Desktop/Target_pipeline/newqual/step16 | mkdir /home/zsolt/Desktop/Target_pipeline/newqual/GenomeDB";


#Variables:
#$bwa=directory of bwa
$bwa='/programs/bwa-0.7.12/bwa';
#$trimmomatic= directory of trimmomatic
$trimmomatic='/programs/Trimmomatic-0.36/trimmomatic-0.36.jar';
#$picard= directory of picard.jar
$picard='/home/zsolt/programs/picard-2.17.10/picard.jar';
#$fastqc= directory of picard.jar
$fastqc='/home/zsolt/programs/FastQC/fastqc';
#$samtools= directory of samtools
$samtools='/programs/samtools-1.4/samtools';
#$oldGATK= directory of the old version of GATK (gatk jar file)
$oldGATK=' /programs/gatk-3.7/GenomeAnalysisTK.jar';
#$GATK= directory of newest version of GATK (gatk executable file)
$GATK='/home/zsolt/programs/gatk-4.0.1.2/gatk';
#$Ref_genome= directory of reference genome (human_g1k_v37.fasta)
$Ref_genome='/home/zsolt/Desktop/Genomes_pipeline/reference_files/human_g1k_v37.fasta';
#$Known_Indel= directory of known indels database (Mills_and_1000G_gold_standard.indels.b37.vcf)
$Known_Indel= '/home/zsolt/Desktop/Genomes_pipeline/reference_files/Mills_and_1000G_gold_standard.indels.b37.vcf';
#$Known_SNP= directory of known SNPs database (dbsnp_138.b37.vcf)
$Known_SNP= '/home/zsolt/Desktop/Genomes_pipeline/reference_files/dbsnp_138.b37.vcf';
#$Gen_dir= directory of target sequences input files (.fastq)
$Gen_dir= '/home/rennanmoreira/Desktop/Target_pipeline/Sequencias';
#$TargetList = directory of the file which contains coordinates of the targeted regions
$TargetList= '/home/zsolt/Desktop/Target_pipeline';
$out_step01= '/home/zsolt/Desktop/Target_pipeline/newqual/step01';
$out_step02= '/home/zsolt/Desktop/Target_pipeline/newqual/step02';
$out_step03= '/home/zsolt/Desktop/Target_pipeline/newqual/step03';
$out_step04= '/home/zsolt/Desktop/Target_pipeline/newqual/step04';
$out_step05= '/home/zsolt/Desktop/Target_pipeline/newqual/step05';
$out_step06= '/home/zsolt/Desktop/Target_pipeline/newqual/step06';
$out_step07= '/home/zsolt/Desktop/Target_pipeline/newqual/step07';
$out_step08= '/home/zsolt/Desktop/Target_pipeline/newqual/step08';
$out_step09= '/home/zsolt/Desktop/Target_pipeline/newqual/step09';
$out_step10= '/home/zsolt/Desktop/Target_pipeline/newqual/step10';
$out_step11= '/home/zsolt/Desktop/Target_pipeline/newqual/step11';
$out_step12= '/home/zsolt/Desktop/Target_pipeline/newqual/step12';
$out_step13= '/home/zsolt/Desktop/Target_pipeline/newqual/step13';
$out_step14= '/home/zsolt/Desktop/Target_pipeline/newqual/step14';
$out_step15= '/home/zsolt/Desktop/Target_pipeline/newqual/step15';
$out_step16= '/home/zsolt/Desktop/Target_pipeline/newqual/step16';


#############################################################################
## Reading the file names
##
#############################################################################

print "Reading the files...\n";

opendir(directory, $Gen_dir);
@list = readdir(directory);
closedir(directory);

#print Dumper @list;

for($i=0; $i<= $#list ;$i++)
{
	if(@list[$i] ne "." and @list[$i] ne "..")
	{
		@names=split("_L",$list[$i]);
		$names2{@names[0]}=0;
	}
}

@names2=keys(%names2);
print Dumper $names2[0];
#print "End of reading files...\n";

#print STDERR "Hit Enter > "; scalar <STDIN>;


#############################################################################
## Step 1 - Trimommatic
## DIFFERENCE FROM BEST PRACTICES: "HEADCROP:5" as a trimming option
#############################################################################

print "Executing Trimmomatic...\n";

for($i=0; $i<= $#names2 ;$i++)
{
	system "mkdir $out_step01/$names2[$i]_L001_R1_001/ | mkdir $out_step01/$names2[$i]_L001_R2_001 | mkdir $out_step01/$names2[$i]_L001_R1_trimmed | mkdir $out_step01/$names2[$i]_L001_R2_trimmed";
	system "$fastqc -o $out_step01/$names2[$i]_L001_R1_001 $Gen_dir/$names2[$i]_L001_R1_001.fastq.gz";
	system "$fastqc -o $out_step01/$names2[$i]_L001_R2_001 $Gen_dir/$names2[$i]_L001_R2_001.fastq.gz";
	system "java -jar $trimmomatic PE -threads 8 -phred33 $Gen_dir/$names2[$i]_L001_R1_001.fastq.gz $Gen_dir/$names2[$i]_L001_R2_001.fastq.gz $out_step01/$names2[$i]_file1.fq $out_step01/$names2[$i]_file2.fq $out_step01/$names2[$i]_file3.fq $out_step01/$names2[$i]_file4.fq ILLUMINACLIP:/home/zsolt/Desktop/Target_pipeline/newqual/haloplexadapters.fa:2:30:12 LEADING:12 TRAILING:12 SLIDINGWINDOW:4:15 MINLEN:36 HEADCROP:5";
	system "$fastqc -o $out_step01/$names2[$i]_L001_R1_trimmed $out_step01/$names2[$i]_file1.fq";
	system "$fastqc -o $out_step01/$names2[$i]_L001_R2_trimmed $out_step01/$names2[$i]_file3.fq";
}

print "End\n";

#print STDERR "Hit Enter > "; scalar <STDIN>;

#############################################################################
## Reading files names from trimmomatic
##
#############################################################################

#print "Reading the files...\n";

opendir(directory, $out_step1);
@lista = readdir(directory);
closedir(directory);

#print Dumper @lista;

for($i=0; $i<= $#lista ;$i++)
{
	if(@lista[$i] ne "." and @lista[$i] ne "..")
    {
		@names=split("_f",$lista[$i]);
		$names2{@names[0]}=0;
    }
}

@names2=keys(%names2);

#print Dumper $names2[0];
#print "End of reading files...\n";


#############################################################################
## Preparing to run BWA
##
#############################################################################

print "Executing BWA...\n";
system "$bwa index -a bwtsw $Ref_genome";
system "$samtools faidx $Ref_genome";
system "java -jar $picard CreateSequenceDictionary REFERENCE=$Ref_genome OUTPUT=/home/zsolt/Desktop/Target_pipeline/newqual/human_g1k_v37.dict";


#############################################################################
## Starting the analysis per sample
##
#############################################################################

for($i=0; $i<= $#names2 ;$i++)
{

#############################################################################
## Step 2 - Running BWA
##
#############################################################################

	print "Running BWA...\n";
	system "$bwa mem -M -t 16 $Ref_genome $out_step01/$names2[$i]_file1.fq $out_step01/$names2[$i]_file3.fq | $samtools import /home/zsolt/Desktop/Genomes_pipeline/reference_files/human_g1k_v37.fasta.fai - - | $samtools sort - -o $out_step02/$names2[$i]_out_temp.bam";
	system "$samtools index $out_step02/$names2[$i]_out_temp.bam";

	
#############################################################################
## Step 2.1 Picard - MarkDuplicates module 
## DIFFERENCE FROM BEST PRACTICES: Duplicates are marked, not removed
## "REMOVE_DUPLICATES=false", duplicates can be removed to save storage space
#############################################################################

	system "java -Xmx16g -jar $picard MarkDuplicates INPUT=$out_step02/$names2[$i]_out_temp.bam OUTPUT=$out_step02/$names2[$i]_out_temp_wd.bam METRICS_FILE=$out_step02/$names2[$i]_wd_metrics_file.txt REMOVE_DUPLICATES=false CREATE_INDEX=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT";
	system "$samtools index $out_step02/$names2[$i]_out_temp_wd.bam";

	
#############################################################################
## Step 3 - Running Picard - CollectInsertSizeMetrics module
##
#############################################################################

	print "Running Picard - CollectInsertSizeMetrics module: $i\n";
	system "java -Xmx16g -jar $picard CollectInsertSizeMetrics INPUT=$out_step02/$names2[$i]_out_temp_wd.bam OUTPUT=$out_step03/$names2[$i]_insertsize_metrics.txt HISTOGRAM_FILE=$out_step03/$names2[$i]_OUTPUT_InsertSizeHist.pdf DEVIATIONS=10.0 MINIMUM_PCT=0.05 VALIDATION_STRINGENCY=LENIENT";


#############################################################################
## Step 4 - Running SAM tools Filtering out improperly aligned pairs
##
#############################################################################

	print "Running SAM tools to filter out improperly aligned pairs: $i\n";
	system "$samtools view -b -f 3 $out_step02/$names2[$i]_out_temp_wd.bam > $out_step04/$names2[$i]_out_temp_PAIRED.bam";
	system "$samtools index $out_step04/$names2[$i]_out_temp_PAIRED.bam";


#############################################################################
## Step 5 - Running Picard AddOrReplaceReadGroups module
##
#############################################################################

	print "Picard AddOrReplaceReadGroups module\n";
	system "java -Xmx16g -jar $picard AddOrReplaceReadGroups INPUT=$out_step04/$names2[$i]_out_temp_PAIRED.bam  OUTPUT=$out_step05/$names2[$i]_OutputMOD_bam.bam RGID=$names2[$i] RGLB=$names2[$i] RGPL=ILLUMINA RGPU=N/A RGSM=$names2[$i] RGCN=$names2[$i]";
	print "Running SAMtools";
	system "$samtools index $out_step05/$names2[$i]_OutputMOD_bam.bam";


#############################################################################
## Step 6 - GATK - Analyze patterns of covariation in the sequence dataset
##
#############################################################################

	print "Analyzing patterns of covariation in the sequence dataset: $i+1 \n";
	system "$GATK BaseRecalibrator -R $Ref_genome --disable-read-filter NotDuplicateReadFilter -I $out_step05/$names2[$i]_OutputMOD_bam.bam --known-sites $Known_SNP --known-sites $Known_Indel -O $out_step06/$names2[$i]_recal_data.table";

#############################################################################
## Step 7 - GATK - 3.4)Apply the recalibration to your sequence data)
##
#############################################################################

	print "3.4)Applying the recalibration to sequence data: $i\n";
	system "$GATK ApplyBQSR -R $Ref_genome --disable-read-filter NotDuplicateReadFilter -I $out_step05/$names2[$i]_OutputMOD_bam.bam --bqsr-recal-file $out_step06/$names2[$i]_recal_data.table -O $out_step07/$names2[$i]_recal_reads.bam";


#############################################################################
## Step 8 - Picard - Validate Sam files
##
#############################################################################

	system "java -Xms6000m -jar $picard ValidateSamFile INPUT=$out_step07/$names2[$i]_recal_reads.bam OUTPUT=$out_step08/$names2[$i]_validationreport MODE=VERBOSE IS_BISULFITE_SEQUENCED=false";


#############################################################################
## Step 9 - GATK - SNP-Variant Calling
## DIFFERENCE FROM BEST PRACTICES: "NotDuplicateReadFilter"
#############################################################################

	print "SNP-Variant Calling: $i\n";

	system "$GATK --java-options '-Xmx8g' HaplotypeCaller -R $Ref_genome --disable-read-filter NotDuplicateReadFilter --ERC GVCF --max-reads-per-alignment-start 0 -new-qual -I $out_step07/$names2[$i]_recal_reads.bam -L $TargetList/Target3.interval.list -A FisherStrand -A QualByDepth -A ReadPosRankSumTest -A Coverage -A LikelihoodRankSumTest -O $out_step09/$names2[$i]_output.raw.snps.indels.g.vcf";

	$samplesbam[$i] = "-I $out_step07/$names2[$i]_recal_reads.bam ";

}

print "End of analysis per sample \n\n";


#############################################################################
## 
## Iterating through the intervals in Target3.interval.list
##
#############################################################################


open (Int_list, "$TargetList/Target3.interval.list") or die;
$count = 1;

while ( my $line = <Int_list> ) {
	chomp $line;
	$line=~ s/\r//gi;
	$dbfolder="$TargetList/GenomeDB/Int_$count";

	
#############################################################################
## 
## Step 10 GATK - Creating a database of sample variants
##
#############################################################################

	print "Creating a database of sample variants\n";
	system "$GATK --java-options '-Xmx8g' GenomicsDBImport --genomicsdb-workspace-path $dbfolder --batch-size 50 -L $line --sample-name-map $TargetList/samplemap_all.tsv --reader-threads 5";
	
	
#############################################################################
## Step 11 - GATK - Joint discovery
##
#############################################################################

	print "Joint discovery of short variants\n";
	system "$GATK --java-options '-Xmx8g' GenotypeGVCFs -R $Ref_genome -O $out_step11/interval$count.vcf -D $Known_SNP -G StandardAnnotation --only-output-calls-starting-in-intervals -V gendb://$dbfolder -new-qual -L $line"; 
	$count ++;
}

print "Ran correctly for $count-1 intervals\n";

#############################################################################
## Step 12 - GATK - Merge and annotate variants
##
#############################################################################
print "Merging and annotating variants\n";
system "java -Xmx16g -jar $picard MergeVcfs I=$TargetList/intervalmap.list O=$out_step12/1_150_output_raw_GenotypeGVCF.vcf";
system "java -Xmx16g -jar $oldGATK -T VariantAnnotator -R $Ref_genome @samplesbam -o $out_step12/1_150_output_raw_annotated.vcf -V $out_step12/1_150_output_raw_GenotypeGVCF.vcf -A MappingQualityRankSumTest -A ReadPosRankSumTest";

#############################################################################
## GATK - Validate VCF file
##
#############################################################################

print "Validate VCF file\n";
system "$GATK --java-options '-Xmx8g' ValidateVariants -V $out_step12/1_150_output_raw_annotated.vcf -R $Ref_genome -L $TargetList/Target3.interval.list --validation-type-to-exclude ALLELES --dbsnp $Known_SNP";


#############################################################################
## Step 13 - GATK - Extract SNPs from the call set
##
#############################################################################

print "Extract SNPs from the call set\n";
system "$GATK SelectVariants -R $Ref_genome -V $out_step12/1_150_output_raw_annotated.vcf -L $TargetList/Target3.interval.list --select-type-to-include SNP -O $out_step13/1_150_raw_snps_extracted.vcf";


#############################################################################
## Step 14 - GATK - Apply the filter to the SNP call set
## DIFFERENCE FROM BEST PRACTICES: Hard filtering is only advised for small datasets
#############################################################################

print "SNP Filtration\n";
system "java -jar $oldGATK -T VariantFiltration -R $Ref_genome -V $out_step13/1_150_raw_snps_extracted.vcf --filterExpression \"QD < 3.5 || QUAL < 100.0 || FS > 15.0 || MQ < 40.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\" --filterName 'SNPfilter' -o $out_step14/1_150_filtered_snps.vcf"; 
system "java -Xms16000m -jar $picard CollectVariantCallingMetrics INPUT=$out_step14/1_150_filtered_snps.vcf OUTPUT=$out_step14/1_150_filtered_snpsmetrics.txt DBSNP=/home/zsolt/Desktop/Genomes_pipeline/reference_files/dbsnp_138.b37.vcf SEQUENCE_DICTIONARY=$TargetList/human_g1k_v37.dict";

#############################################################################
## Step 15 - GATK - Extract the Indels from the call set
##
#############################################################################

print "Extract the Indels from the call set\n";
system "$GATK SelectVariants -R $Ref_genome -V $out_step12/1_150_output_raw_annotated.vcf -L $TargetList/Target3.interval.list  --select-type-to-include INDEL -O $out_step15/1_150_raw_indels_extracted.vcf";


#############################################################################
## Step 16 - GATK - Apply the filter to the Indel call set
## DIFFERENCE FROM BEST PRACTICES: Hard filtering is only advised for small datasets
#############################################################################

print "Apply the filter to the Indel call set\n";
system "java -jar $oldGATK -T VariantFiltration -R $Ref_genome -V $out_step15/1_150_raw_indels_extracted.vcf --filterExpression \"QD < 3.5 || QUAL < 100.0 || FS > 30.0 || SOR > 10.0 || ReadPosRankSum < -2.5\" --filterName 'indelfilter' -o $out_step16/1_150_filtered_indels.vcf"; 
system "java -Xms16000m -jar $picard CollectVariantCallingMetrics INPUT=$out_step16/1_150_filtered_indels.vcf OUTPUT=$out_step16/1_150_filtered_indelsmetrics.txt DBSNP=$Known_SNP SEQUENCE_DICTIONARY=$TargetList/human_g1k_v37.dict";

print "End of the whole pipeline :)   \n";