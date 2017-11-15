#!/usr/bin/env Rscript
#########################################################

#Author: Tjitske de Vries
#date:   07-11-2017
#Name:   reorder-breakmer004-snake.R

# This script is written to reorder breakmer output data.
# Based upon output of Breakmer 0.0.4

# 171107-reorder-breakmer004-snake.R

#########################################################
#               load required libraries                 #
#########################################################
source("code/reorder_functions.R")
suppressMessages(library(data.table))
suppressMessages(library(splitstackshape))
suppressMessages(library(plyr))
#########################################################
#               Locations annotation files              #
#########################################################

Annotationfile<-snakemake@params[["Annotationfile"]]  # contain annotations of areas in genome that are in capturepanel BCNHLv2.
Annotations<-read.delim(Annotationfile, stringsAsFactors = FALSE, sep = "\t" , header=FALSE ) # V1=chr V2=start V3=stop V4=name

#########################################################
#   Functions to annotate genes, calculate length and split columns		#
#########################################################
#
# split_columns<-function(df){
# 	#Split the split read count column into separate columns and place them at the end of the datatable
# 	df<-cSplit_f(df, c("split_read_count", "strands"), sep=",")
#
# 	#Split the breakpoint positions and genes in seperate colums
# 	df<-cSplit_f(df, "target_breakpoints", sep="," )
# 	df<-cSplit_f(df, c("target_breakpoints_1", "target_breakpoints_2"), sep=":" )
# 	df<-cSplit(df, c("genes","breakpoint_coverages"), sep=",")
# 	return(df)
# }

# #Calculate the length of the structural variant
# svlen<-function(df){
# 	svlen<-rep(NA_character_, nrow(df))
# 	len<-abs(as.numeric(df$POS)-as.numeric(df$POS2))
# 	svlen[df$CHROM == df$CHROM2]<-len[df$CHROM == df$CHROM2]
# 	return(svlen)
# }

# annotate_specific<-function(row){
# 	gene<-which(row[1]==Annotations$V1 & as.numeric(row[2]) > Annotations$V2 & as.numeric(row[2]) <= Annotations$V3)
# 	ifelse(length(gene)!=0,Annotations$V4[gene],row[3])
# }

#########################################################
#            Function to sort lexographically           #
#########################################################
#
# orderSvLexo<-function(df){
# 	#sort CHROM 1 and 2 lexographically, turn columns around if order need to be changed
# 		df[df$CHROM > df$CHROM2, c("CHROM", 'POS', 'CHROM2','POS2', 'GENE','GENE2', 'SR','SR2',"DR","DR2","STRAND","STRAND2", "BRKPT_COV","BRKPT_COV2" )]<-df[df$CHROM > df$CHROM2, c('CHROM2','POS2', "CHROM",'POS','GENE2','GENE','SR2','SR',"DR2","DR","STRAND2","STRAND","BRKPT_COV2", "BRKPT_COV")]
#
# 	#sort on position if chrom 1 and 2 are equal, turn columns around if order need to be changed.
# 		df[df$CHROM == df$CHROM2 & df$POS > df$POS2, c("CHROM", 'POS', 'CHROM2','POS2', 'GENE','GENE2', 'SR','SR2',"DR","DR2", "BRKPT_COV","BRKPT_COV2", "STRAND","STRAND2")]<-df[df$CHROM == df$CHROM2 & df$POS > df$POS2, c('CHROM2','POS2', "CHROM",'POS','GENE2','GENE','SR2','SR',"DR2","DR","BRKPT_COV2", "BRKPT_COV", "STRAND2","STRAND")]
# 	return(df)
# }

#####################################################
#   Main Order BreaKmer functions					#
#####################################################
orderBreakmer<-function(inputfile, output1, output2, output3){

	#Import the file given by the filepath
	input<-read.delim(inputfile, stringsAsFactors = FALSE)

	#Change the characters of the discordant read counts to numbers - None becomes NA
	input<-transform(input, disc_read_count=as.numeric(disc_read_count))

	#BreaKmer also calls complex rearrangments consisting of more than one breakpoint, these contain a dash in the breakpointfield.
	#At the moment it is not possible to compare those with output from the other tools, moreover most often these are false positive results.
	#Keep only dataframe with simple structural variants "ssv"
	ssv<-input[!grepl("-",input$target_breakpoints),]

	ssv$ID<-paste("BreaKmer",seq(1:length(ssv[,1])), sep="-")
	ssv<-split_columns(ssv)

	#Rename the columns so that output of different programs is equal
	ssv<-rename(ssv, replace= c( "split_read_count_1"="SR",  "split_read_count_2"="SR2", "target_breakpoints_1_1"="CHROM", "target_breakpoints_1_2"="POS", "target_breakpoints_2_1" = "CHROM2", "target_breakpoints_2_2"="POS2", "genes_1"="GENE", "genes_2"="GENE2", "sv_subtype"="SVTYPE", "disc_read_count"="DR", "breakpoint_coverages_1"="BRKPT_COV", 	"breakpoint_coverages_2"="BRKPT_COV2", "contig_seq"="CONTIG", "contig_id"="CONTIG_ID",  "mismatches"="MISMATCHES" ,"strands_1"="STRAND", "strands_2"="STRAND2","total_matching"="MATCHING"))

	ssv$SVLEN<-svlen(ssv)

	#Improve the annotation using the annotationfile
	ssv$GENE<-apply(ssv[,c("CHROM","POS","GENE")],1,annotate_specific)
	ssv$GENE2<-apply(ssv[,c("CHROM2","POS2","GENE2")],1,annotate_specific)

	ssv$CHROM<-paste("chr",ssv[,CHROM], sep="")
	ssv$CHROM2<-paste("chr",ssv[,CHROM2], sep="")

	#Change naming of svtypes
	ssv$SVTYPE<-revalue(ssv$SVTYPE, replace = c("trl"="TRL","tandem_dup"="DUP","None"="UNKNOWN","inversion"="INV"))

	#create columns that doesn't exist in BreaKmer output for comparison compatibility with other tools
	ssv$DR2<-as.numeric(rep(NA_character_, nrow(ssv)))
	ssv$QUAL<-rep(NA_character_, nrow(ssv))

	ssv<-orderSvLexo(ssv)
	ssv$TOOL<-"breakmer"
	#Change the order of the columns of the file
	ssv<-ssv[,c('CHROM', 'POS', 'CHROM2','POS2','GENE','GENE2','SVTYPE', 'SR', 'SR2', 'DR','DR2','SVLEN',"BRKPT_COV","BRKPT_COV2","STRAND","STRAND2",'QUAL', "ID","TOOL","MISMATCHES", "MATCHING","CONTIG_ID","CONTIG")]

	ssv<-arrange(ssv,CHROM,POS,CHROM2,POS2,desc(DR),desc(DR2),desc(SR),desc(SR2))   #sort vcf lexographically

	dups<-ssv[duplicated(ssv[,1:7]),]
	#remove all duplicate entries
	ssv2<-ssv[!duplicated(ssv[,1:7]),]

	write.table(ssv, file=output1, sep="\t",row.names=FALSE)
	write.table(ssv2, file=output2, sep="\t",row.names=FALSE)

#########################################################################
#      Keep complex rearrangements in separate file						#
#########################################################################
	complex_rearrangements<-input[grepl("-",input$target_breakpoints),]

	complex_rearrangements<-cSplit(complex_rearrangements, c("split_read_count"), sep=",")

#Reorder the results on structural variant subtype, discordant read count and split read count
	if(length(complex_rearrangements)>20) {
  	complex_rearrangements<-arrange(complex_rearrangements,desc(sv_subtype), desc(disc_read_count), desc(split_read_count_01),desc(split_read_count_02), desc(split_read_count_03), desc(split_read_count_04))
	} else {
		complex_rearrangements<-arrange(complex_rearrangements,desc(sv_subtype), desc(disc_read_count), desc(split_read_count_1), desc(split_read_count_2), desc(split_read_count_3), desc(split_read_count_4))
	}

	SRnames<-names(complex_rearrangements)
	SRnames<-SRnames[grep("split_read_count",names(complex_rearrangements))]

	complex_rearrangements<-complex_rearrangements[,c('target_breakpoints','genes','sv_subtype',SRnames,'disc_read_count', 'breakpoint_coverages', 'total_matching', 'mismatches', 'strands', 'contig_id', 'contig_seq'),with=FALSE]

	write.table(complex_rearrangements, file=output3, sep="\t",row.names=FALSE)
}

#TODO Improve complex rearrangements
orderBreakmer(inputfile=snakemake@input[[1]], output1=snakemake@output[["dups"]],output2=snakemake@output[["ordered"]], output3=snakemake@output[["complexe"]])
