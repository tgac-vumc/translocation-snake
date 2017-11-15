
######## Snakemake header ########
library(methods)
Snakemake <- setClass(
    "Snakemake",
    slots = c(
        input = "list",
        output = "list",
        params = "list",
        wildcards = "list",
        threads = "numeric",
        log = "list",
        resources = "list",
        config = "list",
        rule = "character"
    )
)
snakemake <- Snakemake(
    input = list('../merged/TVU12-11384_S27_L004-trl_summary.csv', "summary" = '../merged/TVU12-11384_S27_L004-trl_summary.csv'),
    output = list('../reports/TVU12-11384_S27_L004-circlize.png', "circlize" = '../reports/TVU12-11384_S27_L004-circlize.png'),
    params = list(),
    wildcards = list('TVU12-11384_S27_L004', "sample" = 'TVU12-11384_S27_L004'),
    threads = 1,
    log = list(),
    resources = list(),
    config = list("gridss" = list("GRIDSS_JAR" = '/net/nfs/PAT/home/tjitske/GRIDSS/gridss-1.4.2-jar-with-dependencies.jar'), "all" = list("REF_CHR" = '/net/nfs/PAT/home/tjitske/files/ref/incl_chr/hg19.fa', "REF_NOCHR" = '/net/nfs/PAT/home/tjitske/files/ref/excl_chr/hg19_refgen_exclchr.fa', "Annotationfile" = 'code/annotation_new_BCNHL-v2.bed', "THREADS" = 4, "targets" = 'code/BCNHLv2_primary_coord.bed'), "novobreak" = list("NOVOBREAK" = '/net/nfs/PAT/home/tjitske/novobreak/nb_distribution/run_novoBreak.sh', "EXE_DIR" = '/net/nfs/PAT/home/tjitske/novobreak/nb_distribution', "NORMAL" = '/net/nfs/PAT/analysis/MPS-310/Blood_combined_bam/combined_healthy_blood_nochr.bam', "HEADER" = '/net/nfs/PAT/home/tjitske/novobreak/nb_distribution/header.txt'), "wham" = list("TRAIN" = '/net/nfs/PAT/home/tjitske/WHAM/wham/data/WHAM_training_data.txt', "WHAM" = '/net/nfs/PAT/home/tjitske/WHAM/wham/bin/wham', "BASEQUAL" = 5, "MAPQUAL" = 10, "CLASSIFY" = '/net/nfs/PAT/home/tjitske/WHAM/wham/utils/classify_WHAM_vcf.py'), "filters" = list("blacklist" = 'code/blacklist.bed', "novoscore" = 4, "gridssscore" = 450, "discread" = 3, "highSVevidence" = 15, "splitread" = 4)),
    rule = 'circlize'
)
######## Original script #########
#!/usr/bin/env Rscript
#########################################################
#Author: Tjitske de Vries
#date:   08-11-2017
#Name:   circlize.R

#########################################################
#library(QDNAseq)
suppressMessages(library(Biobase))
suppressMessages(library(circlize))

png(snakemake@output[["circlize"]], width=7, height=7, unit="in", res=300)

circos.par("start.degree" = 90)
par(cex=0.9)
circos.initializeWithIdeogram(plotType = c("labels", "ideogram"), track.height=convert_height(2, "mm"))

translocations<-read.delim(snakemake@input[["summary"]], stringsAsFactors = F, header = T)

side_one<-translocations[translocations$TPorFP=="TP",c("CHROM","POS")]
side_one[,3]<-side_one[,2]+1
side_one[,4]<-translocations[translocations$TPorFP=="TP","GENE"]
side_two<-translocations[translocations$TPorFP=="TP",c("CHROM2","POS2")]
side_two[,3]<-side_two[,2]+1
side_two[,4]<-translocations[translocations$TPorFP=="TP","GENE2"]
names(side_two)<-names(side_one)
genes<-rbind(side_one,side_two)
genes<-genes[!duplicated(genes[,4]),]

circos.genomicLabels(genes, labels.column = 4, cex=0.7, side="outside",  connection_height =convert_height(1.5, "mm") , labels_height = min(c(convert_height(2.5, "cm"))))

circos.genomicLink(side_one, side_two)

dev.off()
