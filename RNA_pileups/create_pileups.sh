#!/bin/bash
#SBATCH --job-name=testingArrayJobs
#SBATCH --ntasks=1
#SBATCH --error=logs/%j_err.log
#SBATCH --output=logs/%j_out.log
#SBATCH --mem=8G
#SBATCH --time=07:59:59
#SBATCH --array=9,13,14,16

#module purge
#module load bcftools
#module load samtools

config="/nfs/home/vgrether/test/ffpe_filtering/temp/dlbcl_temp.txt"
# Extract the sample name, the vcf_file, and the bam_file for the current $SLURM_ARRAY_TASK_ID from the config file
sample=$(head -n $SLURM_ARRAY_TASK_ID $config | tail -n 1 | awk '{print $2}')
vcf_file=$(head -n $SLURM_ARRAY_TASK_ID $config | tail -n 1 | awk '{print $3}')
bam_file=$(head -n $SLURM_ARRAY_TASK_ID $config | tail -n 1 | awk '{print $4}')
dir="/gpfs/commons/groups/compbio/projects/FFPE_filtering/vgrether/ebm/DLBCL/pileup_summaries"

#Another way to do the same thing! 
#vcf_file=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)
#bam_file=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $config)
echo "This is array task ${SLURM_ARRAY_TASK_ID}, the sample name is ${sample}."

### Make temporary working directory
if [[ -d $dir/${sample}_temp ]]; then
    rm -r $dir/${sample}_temp
fi
mkdir $dir/${sample}_temp

cp $bam_file $dir/${sample}_temp
bam_file="$dir/${sample}_temp/$(basename $bam_file)"
samtools index $bam_file

bcftools query -f '%CHROM\t%POS\t[ %AF]\n' $vcf_file | \
awk -F"\t" '{ OFS="\t"} {split($3, arr, " "); print $1, $2, arr[2]}' |  bgzip -c > $dir/${sample}_temp/annot.txt.gz

tabix -s1 -b2 -e2 $dir/${sample}_temp/annot.txt.gz

echo -e '##INFO=<ID=AF,Number=1,Type=Float,Description="Population Allele Frequencies">' > $dir/${sample}_temp/hdr.txt

gatk_vcf="$dir/${sample}_temp/${sample}.gatk.vcf"

#sample_name=$(echo $sample | sed 's/--.*//')
sample_name=$(bcftools view -h $vcf_file | tail -n 1 | awk -F"\t" '{print $11}')

bcftools annotate -s $sample_name -a $dir/${sample}_temp/annot.txt.gz -h $dir/${sample}_temp/hdr.txt -c CHROM,POS,INFO/AF $vcf_file > $gatk_vcf

echo "Succesfully created GATK VCF: $gatk_vcf"

sbatch --output="/nfs/home/vgrether/test/ffpe_filtering/RNA_pileups/logs/${sample}.run_pileups.out" run_pileups.sh $sample $gatk_vcf $bam_file $dir

## EXTRA STUFF IGNORE!! 

#bam="/gpfs/commons/projects/p1000/data/Breast/Project_POLY_14658_B01_ESS_RNA_ER/Sample_25023660RNAt/analysis/25023660RNAt.final.bam"
#vcf_file="/gpfs/commons/groups/compbio/projects/FFPE_filtering/vgrether/ebm/SAMBAI/HighConfidence/vcfs/25023660DNAt--25023660-549.snv.indel.final.v6.annotated.ebm.high_confidence.vcf"

#while IFS=$'\t' read -r sample vcf_file bam; do 
#    bam=$(echo $bam | xargs)
#    echo $bam
#    echo $(samtools mpileup -r "chr1:1495725-1495725" $bam)
#done < $input