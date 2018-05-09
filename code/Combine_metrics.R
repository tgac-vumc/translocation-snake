#!/usr/bin/env R
#########################################################

#Author: Tjitske de Vries
#date:   08-05-2018
#Name:   Combine_metrics.R

#this script is to merge the important output of multiple metrics files.
# Several metrics from gridss are used, so this script will not work if gridss is not included in the pipeline.
# The combined metrics file is not a complete overview of all metrics obtained.

 #########################################################

HSmetrics<-snakemake@input[["HSmetrics"]]
Alignment_met<-snakemake@input[["Alignment_met"]]
insert_size_metrics<-snakemake@input[["insert_size_metrics"]]
sv_metrics<-snakemake@input[["sv_metrics"]]
sample=snakemake@wildcards[["sample"]]
output=snakemake@output[["combined"]]

header_HSmetrics<-read.table(HSmetrics , skip=6, nrows=1, stringsAsFactors=F)
metrics<-read.table(HSmetrics, skip=7, nrows=1, col.names=header_HSmetrics[1,1:53])

header_alignmetrics<-read.table(Alignment_met , skip=6, nrows=1, stringsAsFactors=F)
alignmetrics<-read.table(Alignment_met, skip=9, nrows=1, col.names=header_alignmetrics[1,1:22])

header_insertmetrics<-read.table(insert_size_metrics , skip=6, nrows=1, stringsAsFactors=F)
insert_metrics<-read.table(insert_size_metrics, skip=7, nrows=1, col.names=header_insertmetrics[1,1:18])

header_svmetrics<-read.table(sv_metrics , skip=1, nrows=1, stringsAsFactors=F)
svmetrics<-read.table(sv_metrics, skip=3, nrows=1, col.names=header_svmetrics[1,1:14])

all<-cbind(metrics,alignmetrics, insert_metrics, svmetrics)
all$SAMPLE<-sample

all$PCT_STRUCTURAL_VAR_READS=all$STRUCTURAL_VARIANT_READS/all$TOTAL_READS
all$PCT_INDEL_READS=all$INDEL_READS/all$TOTAL_READS
all$PCT_SPLIT_READS=all$SPLIT_READS/all$TOTAL_READS
all$PCT_SOFT_CLIPPED_READS=all$SOFT_CLIPPED_READS/all$TOTAL_READS

summary=all[,c("SAMPLE","TOTAL_READS","PF_UNIQUE_READS","PCT_PF_UQ_READS","PF_UQ_READS_ALIGNED", "PCT_PF_UQ_READS_ALIGNED", "PCT_OFF_BAIT", "MEAN_BAIT_COVERAGE", "MEAN_TARGET_COVERAGE", "MEDIAN_TARGET_COVERAGE","PCT_USABLE_BASES_ON_BAIT", "PCT_USABLE_BASES_ON_TARGET", "FOLD_ENRICHMENT", "ZERO_CVG_TARGETS_PCT", "PCT_TARGET_BASES_30X", "PCT_TARGET_BASES_100X" ,"PCT_CHIMERAS","PCT_ADAPTER","MEDIAN_INSERT_SIZE","MEDIAN_ABSOLUTE_DEVIATION", "MEAN_INSERT_SIZE", "STANDARD_DEVIATION", "STRUCTURAL_VARIANT_READS", "PCT_STRUCTURAL_VAR_READS","STRUCTURAL_VARIANT_READ_PAIRS","INDEL_READS","PCT_INDEL_READS","SPLIT_READS","PCT_SPLIT_READS","SOFT_CLIPPED_READS","PCT_SOFT_CLIPPED_READS","UNMAPPED_READS","DISCORDANT_READ_PAIRS","UNMAPPED_MATE_READS","STRUCTURAL_VARIANT_READ_ALIGNMENTS","INDEL_READ_ALIGNMENTS","SPLIT_READ_ALIGNMENTS","SOFT_CLIPPED_READ_ALIGNMENTS","DISCORDANT_READ_PAIR_ALIGNMENTS","UNMAPPED_MATE_READ_ALIGNMENTS","GC_DROPOUT","AT_DROPOUT","PCT_EXC_DUPE", "PCT_EXC_MAPQ", "PCT_EXC_BASEQ", "PCT_EXC_OVERLAP", "PCT_EXC_OFF_TARGET")]

write.table(summary, file=output, col.names=TRUE, row.names=FALSE, quote=FALSE)

#HSmetrics<-"../CovMetrics/DLBCL-27_S38_L003_HSmetrics.txt"
#Alignment_met<-"../gridss/DLBCL-27_S38_L003/DLBCL-27_S38_L003_coordsorted.bam.gridss.working/DLBCL-27_S38_L003_coordsorted.bam.alignment_summary_metrics"
#insert_size_metrics<-"../gridss/DLBCL-27_S38_L003/DLBCL-27_S38_L003_coordsorted.bam.gridss.working/DLBCL-27_S38_L003_coordsorted.bam.insert_size_metrics"
#sv_metrics<-"../gridss/DLBCL-27_S38_L003/DLBCL-27_S38_L003_coordsorted.bam.gridss.working/DLBCL-27_S38_L003_coordsorted.bam.sv_metrics"
#sample="DLBCL-27_S38_L003"
#output="../CovMetrics/DLBCL-27_S38_L003_combinedmetrics.txt"

#HSmetrics<-"../CovMetrics/DLBCL-44_S28_L003_HSmetrics.txt"

#Alignment_met<-"../gridss/DLBCL-44_S28_L003/DLBCL-44_S28_L003_coordsorted.bam.gridss.working/DLBCL-44_S28_L003_coordsorted.bam.alignment_summary_metrics"
#insert_size_metrics<-"../gridss/DLBCL-44_S28_L003/DLBCL-44_S28_L003_coordsorted.bam.gridss.working/DLBCL-44_S28_L003_coordsorted.bam.insert_size_metrics"
#sv_metrics<-"../gridss/DLBCL-44_S28_L003/DLBCL-44_S28_L003_coordsorted.bam.gridss.working/DLBCL-44_S28_L003_coordsorted.bam.sv_metrics"

#HSmetrics<-"CovMetrics/Hyper_Official_D5_All_Prep_S30_HSmetrics.txt"
#Alignment_met<-"gridss/Hyper_Official_D5_All_Prep_S30/Hyper_Official_D5_All_Prep_S30_coordsorted.bam.gridss.working/Hyper_Official_D5_All_Prep_S30_coordsorted.bam.alignment_summary_metrics"
#insert_size_metrics<-"gridss/Hyper_Official_D5_All_Prep_S30/Hyper_Official_D5_All_Prep_S30_coordsorted.bam.gridss.working/Hyper_Official_D5_All_Prep_S30_coordsorted.bam.insert_size_metrics"
#sv_metrics<-"gridss/Hyper_Official_D5_All_Prep_S30/Hyper_Official_D5_All_Prep_S30_coordsorted.bam.gridss.working/Hyper_Official_D5_All_Prep_S30_coordsorted.bam.sv_metrics"
#sample="Hyper_Official_D5_All_Prep_S30"
