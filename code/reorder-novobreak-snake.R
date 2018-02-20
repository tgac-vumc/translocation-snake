#!/usr/bin/env Rscript
#########################################################

#Author: Tjitske de Vries
#date:   07-11-2017
#Name:   reorder-novobreak-snake.R

# This script is written to reorder Novobreak output data as part of the translocation-snake pipeline.

# reorder-novobreak-snake.R

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

capture_targets<-get_ext_capture_targets(snakemake@params[["targets"]])


#########################################################
#            reordering novobreak	                #
#########################################################

orderNovobreak<-function(inputfile, output1, output2){

	bedr<-read.vcf(inputfile, split.info = T)	#the extra columns of the vcf from novobreak are not read by this function, so excols is required to read those. however this is easy to remove infofield prefixes.
	vcf<-bedr$vcf

	excols<-tryCatch(read.delim(inputfile, header = F, sep = '\t', stringsAsFactors = F, skip=23), error=function(e) NULL)

	if(!is.null(excols)){
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

}else{vcf<-data.frame(CHROM=character(), POS=integer(), CHROM2=character(),POS2=integer(),GENE=character(),GENE2=character(),
	SVTYPE=character(),SR=integer(),SR2=integer(),DR=integer(),DR2=integer(),SVLEN=integer(), BRKPT_COV=integer(), BRKPT_COV2=integer(),
	STRAND=character(),STRAND2=character(), QUAL=numeric(),ID=character(),TOOL=character(), CT=integer(), MAPQ=numeric() ,
	CLUSTER_ID=character(), CONTIG_NUM=integer() ,SIZE=integer(),  READS_ASSEMBLY=integer(), TUM_BRKPT1_QUAL=numeric(),TUM_BRKPT1_SR_HIGH_QUAL =numeric(),
	TUM_BRKPT1_SR_QUAL=numeric(), TUM_BRKPT2_QUAL=numeric(),  TUM_BRKPT2_SR_HIGH_QUAL=integer() , TUM_BRKPT2_SR_QUAL=numeric(),
	NORM_BRKPT1_DEP=integer(), NORM_SR=integer(),NORM_BRKPT1_QUAL=numeric() , NORM_BRKPT1_SR_HIGH_QUAL=numeric(), NORM_BRKPT1_SR_QUAL=numeric(),
	NORM_BRKPT2_DEP=integer() ,NORM_SR2=integer() ,NORM_BRKPT2_QUAL=numeric() ,NORM_BRKPT2_SR_HIGH_QUAL=numeric(), NORM_BRKPT2_SR_QUAL=numeric() ,
	NORM_DR=integer(), NORM_DR2=integer(), FILTER=character(), CONTIG=character())

	vcf2<-vcf

}
write.table(vcf, file=output1,row.names=FALSE, sep="\t")
write.table(vcf2, file=output2,row.names=FALSE, sep="\t")
}

orderNovobreak(inputfile=snakemake@input[["vcf"]], output1=snakemake@output[["dups"]],output2=snakemake@output[["ordered"]])
