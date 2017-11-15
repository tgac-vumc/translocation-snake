#!/usr/bin/env Rscript
#########################################################

#Author: Tjitske de Vries
#date:   07-11-2017
#Name:   reorder-novobreak-snake.R

# This script is written to reorder Novobreak output data.

# 171107-reorder-novobreak-snake.R

#########################################################
#               load required libraries                 #
#########################################################
source("code/reorder_functions.R")
suppressMessages(library(splitstackshape))
suppressMessages(library(bedr))
suppressMessages(library(VariantAnnotation))
suppressMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(plyr))

#########################################################
#               Locations annotation files              #
#########################################################

gns <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)

Annotationfile<-snakemake@params[["Annotationfile"]] # contain annotations of areas in genome that are in capturepanel BCNHLv2.
Annotations<-read.delim(Annotationfile, stringsAsFactors = FALSE, sep = "\t" , header=FALSE ) # V1=chr V2=start V3=stop V4=name
Annotations[,1]<-paste("chr",Annotations[,1],sep="")

capture_targets<-get_ext_capture_targets()
# capture_targets<-read.table("code/, sep="\t", stringsAsFactors =F, header= F)
# capture_targets[,2]<-capture_targets[,2]-300
# capture_targets[,3]<-capture_targets[,3]+300

#########################################################
#            functions	to annotate gene                #
#########################################################
#
# makeBed<-function(chr, start){
# 	bed<-data.frame(chr, start, start+1)
# 	colnames(bed)<-c("chr","start","stop")
# 	return(bed)
# }
#
# #this funtion find genes in UCSC hg19 database and returns genesymbol in vector.
# getGene<-function(chr, start){
# 	bed<-makeBed(chr,start)
# 	gr<-makeGRangesFromDataFrame(bed)
# 	hits <- as.data.frame(findOverlaps(gr, gns, ignore.strand=TRUE))
# 	hits$SYMBOL <- biomaRt::select(org.Hs.eg.db, gns[hits$subjectHits]$gene_id, "SYMBOL")$SYMBOL
# 	GENE<-rep("", length(chr))
# 	GENE[hits$queryHits] <- hits$SYMBOL
# 	return(GENE)
# }
#
# annotate_specific<-function(row){
# 	gene<-which(row[1]==Annotations$V1 & as.numeric(row[2]) > Annotations$V2 & as.numeric(row[2]) <= Annotations$V3)
# 	ifelse(length(gene)!=0,Annotations$V4[gene],row[3])
# }

#return true if one of the breakpoints fall within the capture region or an area 300 around it.
# captured<- function(row){
# 	position1<-which(row['CHROM'] == capture_targets$V1 & as.numeric(row['POS']) >=  capture_targets$V2 & as.numeric(row['POS']) <=  capture_targets$V3)
# 	position2<-which(row['CHROM2'] == capture_targets$V1 & as.numeric(row['POS2']) >=  capture_targets$V2 & as.numeric(row['POS2']) <=  capture_targets$V3)
# 	return( length(position1) != 0 | length(position2) != 0 )
# }

#########################################################
#            Function to sort lexographically           #
#########################################################

# orderSvLexoNovo<-function(df){
# 	#sort CHROM 1 and 2 lexographically, turn columns around if order need to be changed
# 		df[df$CHROM > df$CHROM2, c("CHROM", 'POS', 'CHROM2','POS2', 'GENE','GENE2', 'SR','SR2',"DR","DR2","STRAND","STRAND2", "BRKPT_COV","BRKPT_COV2", "TUM_BRKPT1_QUAL","TUM_BRKPT1_SR_HIGH_QUAL" , "TUM_BRKPT1_SR_QUAL", "TUM_BRKPT2_QUAL",  "TUM_BRKPT2_SR_HIGH_QUAL" , "TUM_BRKPT2_SR_QUAL","NORM_BRKPT1_DEP", "NORM_SR","NORM_BRKPT1_QUAL" , "NORM_BRKPT1_SR_HIGH_QUAL", "NORM_BRKPT1_SR_QUAL", "NORM_BRKPT2_DEP" ,"NORM_SR2" ,"NORM_BRKPT2_QUAL" ,"NORM_BRKPT2_SR_HIGH_QUAL", "NORM_BRKPT2_SR_QUAL" ,"NORM_DR", "NORM_DR2" )]<-df[df$CHROM > df$CHROM2, c('CHROM2','POS2', "CHROM",'POS','GENE2','GENE','SR2','SR',"DR2","DR","STRAND2","STRAND","BRKPT_COV2", "BRKPT_COV","TUM_BRKPT2_QUAL","TUM_BRKPT2_SR_HIGH_QUAL" , "TUM_BRKPT2_SR_QUAL", "TUM_BRKPT1_QUAL",  "TUM_BRKPT1_SR_HIGH_QUAL" , "TUM_BRKPT1_SR_QUAL","NORM_BRKPT2_DEP", "NORM_SR2","NORM_BRKPT2_QUAL" , "NORM_BRKPT2_SR_HIGH_QUAL", "NORM_BRKPT2_SR_QUAL", "NORM_BRKPT1_DEP" ,"NORM_SR" ,"NORM_BRKPT1_QUAL" ,"NORM_BRKPT1_SR_HIGH_QUAL", "NORM_BRKPT1_SR_QUAL" ,"NORM_DR2", "NORM_DR")]
#
# 	return(df)
# }

#########################################################
#            reordering novobreak	                #
#########################################################

orderNovobreak<-function(inputfile, output1, output2){

	bedr<-read.vcf(inputfile, split.info = T)	#the extra columns of the vcf from novobreak are not read by this function, so excols is required to read those. however this is easy to remove infofield prefixes.
	vcf<-bedr$vcf

	excols<-read.delim(inputfile, header = F, sep="\t", stringsAsFactors = F, skip=23 )
	colnames(excols)[11:39]<-c("CLUSTER_ID","CONTIG_NUM","SIZE","READS-ASSEMBLY", "COVIDENCE", "BRKPT_COV", "SR", "TUM_BRKPT1_QUAL", "TUM_BRKPT1_SR_HIGH_QUAL", "TUM_BRKPT1_SR_QUAL", "NORM_BRKPT1_DEP", "NORM_SR", "NORM_BRKPT1_QUAL", "NORM_BRKPT1_SR_HIGH_QUAL", "NORM_BRKPT1_SR_QUAL", "BRKPT_COV2", "SR2", "TUM_BRKPT2_QUAL", "TUM_BRKPT2_SR_HIGH_QUAL", "TUM_BRKPT2_SR_QUAL", "NORM_BRKPT2_DEP", "NORM_SR2", "NORM_BRKPT2_QUAL", "NORM_BRKPT2_SR_HIGH_QUAL", "NORM_BRKPT2_SR_QUAL", "DR" , "NORM_DR","DR2", "NORM_DR2")

	vcf<-cbind(vcf[,c(1,2,3,6,7,11:13,16,18,19)], excols[,11:39])  #only keep the informative fields of the vcf and add the extra fields

	#Create unique ID per hit
	vcf$ID<-paste("novobreak",seq(1:length(vcf$ID)), sep="-")

	#remove chrM - not expected to be real events.
	vcf<-vcf[vcf$CHROM !="M" & vcf$CHR2 != "M",]

	#remove weird events with higher DR counts than coverages.
	vcf<-vcf[!(vcf$DR > vcf$BRKPT_COV & vcf$DR2 > vcf$BRKPT_COV2),]

	#add chromosome prefix
	vcf$CHROM<-paste("chr",vcf[,"CHROM"], sep="")
	vcf$CHR2<-paste("chr",vcf[,"CHR2"], sep="")

	vcf<-plyr::rename(vcf, replace= c("CHR2"="CHROM2", "END"="POS2", "CONSENSUS"="CONTIG", "QUAL"="MAPQ"))
	vcf<-plyr::rename(vcf, replace= c("COVIDENCE"="QUAL"))
	vcf$QUAL<-as.numeric(gsub("cov", "", vcf$QUAL))

	#remove hits further than 300 bp outside capture areas.
	keeplist<-apply(vcf ,1, captured)
	vcf<-vcf[keeplist,]

	vcf<-as.data.table(vcf)

	#Annotate genes
	vcf$GENE<-getGeneBed(vcf[,"CHROM"], vcf[,"POS"])
	vcf$GENE2<-getGeneBed(vcf[,"CHROM2"],vcf[, "POS2"])

	#Annotate genes that are in the capture panel as translocation targets more specific (exons, upstr gene, IGH etc.)
	vcf$GENE<-apply(vcf[,c("CHROM","POS","GENE")],1,annotate_specific)
	vcf$GENE2<-apply(vcf[,c("CHROM2","POS2","GENE2")],1,annotate_specific)

	vcf$TOOL<-"novobreak"
	vcf$STRAND<-ifelse(substring(vcf$CT, first= 1, last=1)=="3" ,yes="+",no="-")
	vcf$STRAND2<-ifelse(substring(vcf$CT, first=4, last=4)=="3" ,yes="+",no="-")
	vcf$SVTYPE[vcf[,"SVTYPE"] == "TRA"]<-"TRL"
	vcf$POS2<-as.integer(vcf$POS2)

	vcf<-vcf[,c('CHROM', 'POS', 'CHROM2','POS2','GENE','GENE2','SVTYPE','SR','SR2',"DR","DR2","SVLEN", "BRKPT_COV", "BRKPT_COV2","STRAND","STRAND2", 'QUAL',"ID","TOOL", "CT", "MAPQ" ,"CLUSTER_ID", "CONTIG_NUM" ,"SIZE",  "READS-ASSEMBLY", "TUM_BRKPT1_QUAL","TUM_BRKPT1_SR_HIGH_QUAL" , "TUM_BRKPT1_SR_QUAL", "TUM_BRKPT2_QUAL",  "TUM_BRKPT2_SR_HIGH_QUAL" , "TUM_BRKPT2_SR_QUAL","NORM_BRKPT1_DEP", "NORM_SR","NORM_BRKPT1_QUAL" , "NORM_BRKPT1_SR_HIGH_QUAL", "NORM_BRKPT1_SR_QUAL", "NORM_BRKPT2_DEP" ,"NORM_SR2" ,"NORM_BRKPT2_QUAL" ,"NORM_BRKPT2_SR_HIGH_QUAL", "NORM_BRKPT2_SR_QUAL" ,"NORM_DR", "NORM_DR2", "FILTER", "CONTIG")]

	vcf<-orderSvLexoNovo(vcf)	#sort each SV lexographically,

	vcf<-arrange(vcf, CHROM,POS,CHROM2,POS2,plyr::desc(DR),plyr::desc(DR2),plyr::desc(SR),plyr::desc(SR2))   #sort vcf lexographically

	#remove all duplicate entries
	vcf2<-vcf[!duplicated(vcf[,1:7]),]

	write.table(vcf, file=output1,row.names=FALSE, sep="\t")
	write.table(vcf2, file=output2,row.names=FALSE, sep="\t")
}

orderNovobreak(inputfile=snakemake@input[[1]], output1=snakemake@output[["dups"]],output2=snakemake@output[["ordered"]])
