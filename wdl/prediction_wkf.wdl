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
import "wdl_structs.wdl"
import "fifa.wdl" as fifa


workflow ExtractionWkf {
    input {
        Bam bam
        String sampleId
        String projectId
        IndexedVcf vcf
        IndexedReference referenceFa
        File? optionalRnaFile
        Array[File] models
        # resources
        String qos = "compbio"
        String partition = "cpu"
        String cpuPlatform = "Intel Cascade Lake"
    }
    Int diskSize = ceil(size(bam.bam, "GB")) * 3
    call fifa.Extraction {
        input:
            bam = bam,
            sampleId = sampleId,
            projectId = projectId,
            vcf = vcf,
            referenceFa = referenceFa,
            diskSize = diskSize,
            qos = qos,
            partition = partition,
            cpuPlatform = cpuPlatform
    }

    if (defined(optionalRnaFile)) {
        call fifa.PredictionWithRna {
            input:
                sampleId = sampleId,
                vcf = vcf,
                extractedFeatures = Extraction.extractedFeatures,
                models = models,
                optionalRnaFile = optionalRnaFile
        }
    }

    if (!defined(optionalRnaFile)) {
        call fifa.Prediction {
            input:
                sampleId = sampleId,
                vcf = vcf,
                extractedFeatures = Extraction.extractedFeatures,
                models = models
        }
    }
    File fifaVcfRun = select_first([PredictionWithRna.fifaVcf, Prediction.fifaVcf])
    output {
        File extractedFeatures = Extraction.extractedFeatures
        File fifaVcf = fifaVcfRun
    }
}
