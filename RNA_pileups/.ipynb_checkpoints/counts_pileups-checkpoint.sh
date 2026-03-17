#!/bin/bash

for sample in $(find $dir -name "*_pileup.txt" -exec basename {} \; | sed 's/_pileup.txt//'); do
    vcf="/gpfs/commons/groups/compbio/projects/FFPE_filtering/vgrether/ebm/SAMBAI/HighConfidence/vcfs/${sample}.snv.indel.final.v6.annotated.ebm.high_confidence.vcf"
    if [[ ! -f $vcf ]]; then 
        echo $sample 
        continue
    fi
    pileup=$(find $dir -name "${sample}_pileup.txt")
    cat $pileup| tail -n +2 > "${sample}_pileup_noHeader.txt"
    echo -e "SAMPLE\tCHROM\tPOS\tEBM" > $dir/${sample}_ebm_labels.txt
    bcftools query -f '%CHROM\t%POS\t%INFO/EBM\n' $vcf | awk -F"\t" '{OFS="\t"; print "'"${sample}"'", $1, $2, $3}'  >> $dir/${sample}_ebm_labels.txt
done


for file in $(find /gpfs/commons/groups/compbio/projects/FFPE_filtering/vgrether/ebm/DLBCL/pileup_summaries/ -name '*_pileup.txt'); do
    cat $file | awk -F'\t' '{print $6}' | sort | awk '{a[i++]=$1} END {print a[int(i/2)]}'
done

mutect2_dir="/gpfs/commons/groups/compbio/projects/FFPE_filtering/wliao/mutect2"
dir="/gpfs/commons/groups/compbio/projects/FFPE_filtering/vgrether/ebm/DLBCL/pileup_summaries"
for file in $(find $dir -name '*_pileup.txt' | xargs realpath); do
    sample=$(basename $file | sed 's/_pileup.txt//')
    echo $sample
    tail +2 $file > "/gpfs/commons/groups/compbio/projects/FFPE_filtering/vgrether/ebm/DLBCL/pileup_summaries/${sample}_pileup_noHeader.txt"
    #vcf_file=$(find $mutect2_dir -name "${sample}*.mutect2_forceCall.labeled.tumorSample.vcf.gz" | xargs realpath)
    #echo $vcf_file
    #label_file="$dir/${sample}.labels.txt"
    #echo -e "Chrom\tPOS\tLabel" > $label_file
    #bcftools query -f '%CHROM\t%POS\t%INFO/Label' $vcf_file >> $label_file
done

output="$dir/all_samples_ffpolish_truthset.txt"
echo -e 'SAMPLE\tCHROM\tPOS\tFFPOLISH' > $output

for vcf_file in $(find $dir --maxdepth 2 -name '*.vcf.gz'); do
    sample=$(basename "$vcf_file" | sed 's|\.vcf\.gz$||')
    sample_name=$(echo $sample | sed 's/--*//')
    echo $sample
    echo $sample_name
    bcftools query -f '%CHROM\t%POS' "$dir/Sample_${sample}/${sample_name}_filtered.vcf" | awk -v sample=$sample -F'\t' '{OFS=\t; print sample, $1, $2, 1}' >> $output
done