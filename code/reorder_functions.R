#!/usr/bin/env R
#########################################################

#Author: Tjitske de Vries
#date:   08-11-2017
#Name:   reorder_functions.R

#this script contain the functions used for reordering the output of Gridss, Breakmer
#Wham and novobreak.

 #########################################################

get_ext_capture_targets<-function(){
  capture_targets<-read.table(snakemake@params[["targets"]], sep="\t", stringsAsFactors =F, header= F) #"code/BCNHLv2_primary_coord.bed"
  capture_targets[,2]<-capture_targets[,2]-300
  capture_targets[,3]<-capture_targets[,3]+300
  return(capture_targets)
}


#####################################################
#      Annotate genes	                              #
#####################################################
#this funtion find genes in UCSC hg19 database and returns genesymbol in vector. input is a generange object
getGeneGr<-function(gr){
	hits <- as.data.frame(findOverlaps(gr, gns, ignore.strand=TRUE))
	hits$SYMBOL <- biomaRt::select(org.Hs.eg.db, gns[hits$subjectHits]$gene_id, "SYMBOL")$SYMBOL
	hits$gene_strand <- as.character(strand(gns[hits$subjectHits]))
	gr$SYMBOL <- ""
	gr$SYMBOL[hits$queryHits] <- hits$SYMBOL
	return(gr)
}

#this funtion find genes in UCSC hg19 database and returns genesymbol in vector. input is a genomic position chromosome and start
getGeneBed<-function(chr, start){
	bed<-makeBed(chr,start)
	gr<-makeGRangesFromDataFrame(bed)
	hits <- as.data.frame(findOverlaps(gr, gns, ignore.strand=TRUE))
	hits$SYMBOL <- biomaRt::select(org.Hs.eg.db, gns[hits$subjectHits]$gene_id, "SYMBOL")$SYMBOL
	GENE<-rep("", nrow(chr))
	GENE[hits$queryHits] <- hits$SYMBOL
	return(GENE)
}

#annotate specific breakpoint location of targetgenes (exons, introns, breakpointregion, IG regions)
annotate_specific<-function(row){
	gene<-which(row[1]==Annotations$V1 & as.numeric(row[2]) > Annotations$V2 & as.numeric(row[2]) <= Annotations$V3)
	ifelse(length(gene)!=0,Annotations$V4[gene],row[3])
}


#####################################################
#      remove events not in capture region          #
#####################################################

# returns false if both breakpoints lay outside capture targets.
captured<- function(row){
	position1<-which(row['CHROM'] == capture_targets$V1 & as.numeric(row['POS']) >=  capture_targets$V2 & as.numeric(row['POS']) <=  capture_targets$V3)
	position2<-which(row['CHROM2'] == capture_targets$V1 & as.numeric(row['POS2']) >=  capture_targets$V2 & as.numeric(row['POS2']) <=  capture_targets$V3)
	return( length(position1) != 0 | length(position2) != 0 )
}

#########################################################
#            Function to sort lexographically           #
#########################################################

orderSvLexo<-function(df){
	#sort CHROM 1 and 2 lexographically, turn columns around if order need to be changed
		df[df$CHROM > df$CHROM2, c("CHROM", 'POS', 'CHROM2','POS2', 'GENE','GENE2', 'SR','SR2',"DR","DR2","STRAND","STRAND2", "BRKPT_COV","BRKPT_COV2" )]<-df[df$CHROM > df$CHROM2, c('CHROM2','POS2', "CHROM",'POS','GENE2','GENE','SR2','SR',"DR2","DR","STRAND2","STRAND","BRKPT_COV2", "BRKPT_COV")]

	#sort on position if chrom 1 and 2 are equal, turn columns around if order need to be changed.
		df[df$CHROM == df$CHROM2 & df$POS > df$POS2, c("CHROM", 'POS', 'CHROM2','POS2', 'GENE','GENE2', 'SR','SR2',"DR","DR2", "BRKPT_COV","BRKPT_COV2","STRAND","STRAND2")]<-df[df$CHROM == df$CHROM2 & df$POS > df$POS2, c('CHROM2','POS2', "CHROM",'POS','GENE2','GENE','SR2','SR',"DR2","DR","BRKPT_COV2", "BRKPT_COV","STRAND2","STRAND")]
	return(df)
}

orderSvLexoNovo<-function(df){
	#sort CHROM 1 and 2 lexographically, turn columns around if order need to be changed
		df[df$CHROM > df$CHROM2, c("CHROM", 'POS', 'CHROM2','POS2', 'GENE','GENE2', 'SR','SR2',"DR","DR2","STRAND","STRAND2", "BRKPT_COV","BRKPT_COV2", "TUM_BRKPT1_QUAL","TUM_BRKPT1_SR_HIGH_QUAL" , "TUM_BRKPT1_SR_QUAL", "TUM_BRKPT2_QUAL",  "TUM_BRKPT2_SR_HIGH_QUAL" , "TUM_BRKPT2_SR_QUAL","NORM_BRKPT1_DEP", "NORM_SR","NORM_BRKPT1_QUAL" , "NORM_BRKPT1_SR_HIGH_QUAL", "NORM_BRKPT1_SR_QUAL", "NORM_BRKPT2_DEP" ,"NORM_SR2" ,"NORM_BRKPT2_QUAL" ,"NORM_BRKPT2_SR_HIGH_QUAL", "NORM_BRKPT2_SR_QUAL" ,"NORM_DR", "NORM_DR2" )]<-df[df$CHROM > df$CHROM2, c('CHROM2','POS2', "CHROM",'POS','GENE2','GENE','SR2','SR',"DR2","DR","STRAND2","STRAND","BRKPT_COV2", "BRKPT_COV","TUM_BRKPT2_QUAL","TUM_BRKPT2_SR_HIGH_QUAL" , "TUM_BRKPT2_SR_QUAL", "TUM_BRKPT1_QUAL",  "TUM_BRKPT1_SR_HIGH_QUAL" , "TUM_BRKPT1_SR_QUAL","NORM_BRKPT2_DEP", "NORM_SR2","NORM_BRKPT2_QUAL" , "NORM_BRKPT2_SR_HIGH_QUAL", "NORM_BRKPT2_SR_QUAL", "NORM_BRKPT1_DEP" ,"NORM_SR" ,"NORM_BRKPT1_QUAL" ,"NORM_BRKPT1_SR_HIGH_QUAL", "NORM_BRKPT1_SR_QUAL" ,"NORM_DR2", "NORM_DR")]
	return(df)
 }

#########################################################
#       Split colums (breakmer)                         #
#########################################################
split_columns<-function(df){
	#Split the split read count column into separate columns and place them at the end of the datatable
	df<-cSplit_f(df, c("split_read_count", "strands"), sep=",")

	#Split the breakpoint positions and genes in seperate colums
	df<-cSplit_f(df, "target_breakpoints", sep="," )
	df<-cSplit_f(df, c("target_breakpoints_1", "target_breakpoints_2"), sep=":" )
	df<-cSplit(df, c("genes","breakpoint_coverages"), sep=",")
	return(df)
}

#########################################################
#      Annotate SV-types	(gridss)                      #
#########################################################
#copied from gridss github, function to add rearrangement type information.
simpleEventType <- function(gr) {
  return(ifelse(seqnames(gr) != seqnames(partner(gr)), "TRL", # inter-chromosomosal
          ifelse(gr$insLen >= abs(gr$svLen) * 0.7, "INS",
           ifelse(strand(gr) == strand(partner(gr)), "INV",
            ifelse(xor(start(gr) < start(partner(gr)), strand(gr) == "-"), "DEL",
             "DUP")))))
}

#########################################################
#   Find paired event and combine in one row (gridss)   #
#########################################################
#find position of breakend in the paired event (other row) and add in a column of first event
getPair<-function(df){
	df$CHROM2<-NA_character_
	df$POS2<-NA_integer_
	df$GENE2<-NA_character_
	df$STRAND2<-NA_character_
	df[,c("CHROM2","POS2","GENE2","SR2", "RP2","STRAND2")]<-df[df$PARID,c("CHROM","POS","GENE","SR", "RP","STRAND")]
	return(df)
}


#########################################################
#      Calculate the length of the structural variant   #
#########################################################

svlen<-function(df){
	svlen<-rep(NA_character_, nrow(df))
	len<-abs(as.numeric(df$POS)-as.numeric(df$POS2))
	svlen[df$CHROM == df$CHROM2]<-len[df$CHROM == df$CHROM2]
	return(svlen)
}

#########################################################
#  make a bed file from chromosome and position         #
#########################################################

makeBed<-function(chr, start){
	bed<-data.frame(chr, start, start+1)
	colnames(bed)<-c("chr","start","stop")
	return(bed)
}
