version 1.0

# ================== COPYRIGHT ================================================
# New York Genome Center
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2026) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.
#
#    Nico Robine (nrobine@nygenome.org)
#    Valentina Grether
#    Zoe R. Goldstein (zgoldstein@nygenome.org)
#    Jennifer M Shelton (jshelton@nygenome.org)
#    Timothy R. Chu (tchu@nygenome.org)
#    William F. Hooper (whooper@nygenome.org)
#    Heather Geiger (hgeigher @nygenome.org)
#    André Corvelo (acorvelo@nygenome.org)
#    Rachel Martini 
#    Melissa B. Davis
# 
#
# ================== /COPYRIGHT ===============================================

# Workflow from https://www.biorxiv.org/content/10.64898/2026.03.10.710815v1
# An explainable boosting machine model for identifying artifacts caused by formalin-fixed paraffin embedding
import "wdl/wdl_structs.wdl"
import "wdl/prediction_wkf.wdl" as extractionWkf


workflow ExtractionWkfs {
    input {
        Array[Bam] bams
        Array[String] sampleIds
        String projectId
        Array[IndexedVcf] vcfs
        IndexedReference referenceFa
        Array[File]? optionalRnaFiles
        Array[File] models
        # resources
        String qos = "compbio"
        String partition = "cpu"
        String cpuPlatform = "Intel Cascade Lake"
    }
    Array[File] rnaFiles = select_first([optionalRnaFiles, []])
    scatter (i in range(length(bams))) {
        if (defined(optionalRnaFiles)) {
            call extractionWkf.ExtractionWkf as rnaExtractionWkf {
                input:
                    bam = bams[i],
                    sampleId = sampleIds[i],
                    projectId = projectId,
                    vcf = vcfs[i],
                    optionalRnaFile = rnaFiles[i],
                    models = models,
                    referenceFa = referenceFa,
                    qos = qos,
                    partition = partition,
                    cpuPlatform = cpuPlatform
            }
        }
        if (!defined(optionalRnaFiles)) {
            call extractionWkf.ExtractionWkf {
                input:
                    bam = bams[i],
                    sampleId = sampleIds[i],
                    projectId = projectId,
                    vcf = vcfs[i],
                    models = models,
                    referenceFa = referenceFa,
                    qos = qos,
                    partition = partition,
                    cpuPlatform = cpuPlatform
            }
        }
        File extractedFeaturesRun = select_first([rnaExtractionWkf.extractedFeatures, ExtractionWkf.extractedFeatures])
        File fifaVcfRun = select_first([rnaExtractionWkf.fifaVcf, ExtractionWkf.fifaVcf])
    }
    output {
        Array[File] extractedFeatures = extractedFeaturesRun
        Array[File] fifaVcf = fifaVcfRun
    }
}
