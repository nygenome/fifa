#!/bin/bash
#SBATCH --job-name=testingArrayJobs
#SBATCH --ntasks=1
#SBATCH --error=logs/%j_out.log
#SBATCH --output=logs/%j_out.log
#SBATCH --mem=24G
#SBATCH --time=07:59:59
#SBATCH --array=3-51

#Should be --array=1-51

#module purge
#module load singularity/4.2.1
#module load bcftools
#module load R

ref_seq='/gpfs/commons/groups/nygcfaculty/kancero/data/bwa-mem2/GRCh38.d1.vd1.fa'
dir='/gpfs/commons/groups/compbio/projects/FFPE_filtering/vgrether'
config="/nfs/home/vgrether/p1000_samples.csv"

tumor=$(head -n $SLURM_ARRAY_TASK_ID $config | tail -n 1 | awk -F, '{print $2}')
normal=$(head -n $SLURM_ARRAY_TASK_ID $config | tail -n 1 | awk -F, '{print $3}')
vcf_FFPE=$(head -n $SLURM_ARRAY_TASK_ID $config | tail -n 1 | awk -F, '{print $5}')
bam_FFPE=$(head -n $SLURM_ARRAY_TASK_ID $config | tail -n 1 | awk -F, '{print $6}')

sample="${tumor}--${normal}"
echo $sample

if [ ! -d "$dir/ffpolish/original_model/Sample_${sample}" ]; then
    mkdir -p "$dir/ffpolish/original_model/Sample_$sample"
fi 
if [ -f "$vcf_FFPE" ]; then
    if [[ "$vcf_FFPE" == *.vcf.gz ]]; then
        :
    else
        bgzip -c "$vcf_FFPE" > "$dir/ffpolish/original_model/Sample_$sample/$sample.vcf.gz"
        vcf_FFPE="$dir/ffpolish/original_model/Sample_$sample/$sample.vcf.gz"
    fi
    sbatch --mem=24G --output="/nfs/home/vgrether/test/ffpe_filtering/RNA_pileups/logs/$sample.out" --wrap="singularity exec -B '/gpfs' docker://matnguyen/ffpolish:0.1.0 ffpolish filter -o '$dir/ffpolish/original_model/Sample_$sample' -p '$tumor' '$ref_seq' '$vcf_FFPE' '$bam_FFPE'"
fi 
