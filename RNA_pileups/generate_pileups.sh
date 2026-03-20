#!/usr/bin/env bash

## Check installation of required tools
for tool in bcftools gatk tabix bgzip; do
    if ! command -v $tool &> /dev/null; then
        echo "Error: Required tool '$tool' is not installed or not in PATH."
        exit 1
    fi
done

sample_name=$1
vcf_path=$2
bam_path=$3
output_dir=$4
delete_temp=$5

## Clear existing temporary dir and make a new one
if [[ -d $output_dir/${sample_name}_temp ]]; then
    rm -r $output_dir/${sample_name}_temp
fi
mkdir $output_dir/${sample_name}_temp

## Check if AF field in your vcf_file

if [ $(bcftools view -h $vcf_path | grep '##INFO=<ID=AF' | wc -l) -eq 0 ]; then
    echo "The VCF file does not contain the Population Allele Frequency field in the INFO column. 
    Will now add the AF field to the VCF file using bcftools annotate. "
    
    bcftools query -f '%CHROM\t%POS\t[ %AF]\n' $vcf_path | \
    awk -F"\t" '{ OFS="\t"} {split($3, arr, " "); print $1, $2, arr[2]}' |  bgzip -c > $output_dir/${sample_name}_temp/annot.txt.gz
    
    tabix -s1 -b2 -e2 $output_dir/${sample_name}_temp/annot.txt.gz

    echo -e '##INFO=<ID=AF,Number=1,Type=Float,Description="Population Allele Frequencies">' > $output_dir/${sample_name}_temp/hdr.txt

    gatk_vcf="$output_dir/${sample_name}_temp/${sample_name}.temp.vcf"

    #vcf_sample_id=$(bcftools view -h $vcf_path | tail -n 1 | awk -F"\t" '{print $11}') # for pulling out NYGC Sample name

    bcftools annotate -s $sample_name -a $output_dir/${sample_name}_temp/annot.txt.gz -h $output_dir/${sample_name}_temp/hdr.txt -c CHROM,POS,INFO/AF $vcf_path > $gatk_vcf

    echo -e "Successfully annotated file, wrote temporary vcf file here:\n$gatk_vcf"

else 
    gatk_vcf=$(echo $vcf_path)
fi 

## WARNING: This will throw an error if the bam file is not indexed

out_file="$output_dir/${sample_name}_pileup_gatk_output.txt"

if [[ ! -f $gatk_vcf ]]; then
    echo "Incorrect path for VCF_File: $gatk_vcf."
    exit 1
elif [[ ! -f $bam_path ]]; then
    echo "Incorrect path for BAM_File: $bam_path."
    exit 1
else
    
    gatk IndexFeatureFile \
    -F $gatk_vcf
    
    gatk GetPileupSummaries \
    --minimum-population-allele-frequency 0 \
    --maximum-population-allele-frequency 1 \
    --disable-read-filter MappingQualityAvailableReadFilter \
    --disable-read-filter MappingQualityNotZeroReadFilter \
    --disable-read-filter NotDuplicateReadFilter \
    --disable-read-filter WellformedReadFilter \
    -I $bam_path \
    -V $gatk_vcf \
    -L $gatk_vcf \
    -O $out_file
    
fi 

if [[ "$delete_temp" == "true" ]]; then
    rm -r "$output_dir/${sample_name}_temp"
fi


tail +2 $out_file > "$output_dir/${sample_name}_RNApileup_final.txt"

echo -e "Succesfully generated RNA pileup file:"
echo -e "$output_dir/${sample_name}_RNApileup_final.txt"
echo -e "Please use this file when generating predictions with FIFA with the \'--rna_annotations\' flag"