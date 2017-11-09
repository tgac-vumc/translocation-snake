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
