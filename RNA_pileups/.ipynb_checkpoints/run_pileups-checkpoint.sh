#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --time=03:59:59

#module purge
#module load bcftools
#module load sam tools

module load gatk

sample=$1
vcf=$2
bam=$3
dir=$4

#echo -e "Sample\tChrom\tPOS\tEBM\tTOTAL_COV\tREF\tCOV_REF\tALT\tCOV_ALT" > "$dir/${sample}_pileup.txt"

out_file="$dir/${sample}_pileup.txt"

if [[ ! -f $vcf ]]; then
    echo "Incorrect path for VCF_File: $vcf."
elif [[ ! -f $bam ]]; then
    echo "Incorrect path for BAM_File: $bam."
else
    indexJOB=$(sbatch --parsable --mem=8G \
    --output="/nfs/home/vgrether/test/ffpe_filtering/RNA_pileups/logs/${sample}.IndexFeatureFile.out" \
    --wrap="gatk IndexFeatureFile \
    -I $vcf")
    
    jobID=$(sbatch --parsable \
    -d afterok:$indexJOB \
    --mem=32G \
    --output="/nfs/home/vgrether/test/ffpe_filtering/RNA_pileups/logs/${sample}.GetPileupSummaries.out" \
    --wrap="gatk GetPileupSummaries \
    --minimum-population-allele-frequency 0 \
    --maximum-population-allele-frequency 1 \
    --disable-read-filter MappingQualityAvailableReadFilter \
    --disable-read-filter MappingQualityNotZeroReadFilter \
    --disable-read-filter NotDuplicateReadFilter \
    --disable-read-filter WellformedReadFilter \
    -I $bam \
    -V $vcf \
    -L $vcf \
    -O $out_file")
    
    sbatch -d afterok:$jobID --wrap="rm -r $dir/${sample}_temp"
fi 









