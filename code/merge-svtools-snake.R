#!/usr/bin/env Rscript
#########################################################

#Author: Tjitske de Vries
#date:   07-11-2017
#Name:   merge-svtools.r

# This script is written to merge the ordered output of breakmer, gridss, novobreak and wham as part of the translocation snake.

#########################################################
#               load required libraries                 #
#########################################################
suppressMessages(library(plyr))
source("code/merge_functions.R")

blacklist<-read.delim(snakemake@config[["filters"]][["blacklist"]], stringsAsFactors=F, header=F)

whamfile=snakemake@input[['wham']]
gridssfile=snakemake@input[['gridss']]
novofile=snakemake@input[['novobreak']]
breakmerfile=snakemake@input[['breakmer']]
output=snakemake@output[['trl']]
output2=snakemake@output[['IG']]
output3=snakemake@output[['otherSVs']]
summary=snakemake@output[['summary']]

#########################################################
#               merge all SVs from the four tools       #
#########################################################

merge_SVs<-function(whamfile,gridssfile,novofile,breakmerfile, output, output2,output3,summary){

	wham<-read.delim(whamfile, stringsAsFactors = F)
	gridss<-read.delim(gridssfile, stringsAsFactors = F)
	novo<-read.delim(novofile, stringsAsFactors = F)
	breakmer<-read.delim(breakmerfile, stringsAsFactors = F)

	#convert columns of the different tools to be of equal characteristics.
	gridss$BRKPT_COV<-as.integer(gridss$BRKPT_COV)
	gridss$BRKPT_COV2<-as.integer(gridss$BRKPT_COV2)

	wham$STRAND<-as.character(wham$STRAND)
	wham$STRAND2<-as.character(wham$STRAND2)
	wham$QUAL<-as.numeric(wham$QUAL)
	breakmer$QUAL<-as.numeric(breakmer$QUAL)

	wham<-filter_wham(wham,snakemake@config[["filters"]][["splitread"]])
	gridss<-filter_gridss(gridss,snakemake@config[["filters"]][["gridssscore"]])
	breakmer<-filter_breakmer(breakmer,snakemake@config[["filters"]][["discread"]], snakemake@config[["filters"]][["splitread"]])
	novo<-filter_novobreak(novo,snakemake@config[["filters"]][["novoscore"]])

	df<-rbind(wham[1:19],gridss[1:19],breakmer[1:19],novo[1:19])
	df<-arrange(df,CHROM,POS,CHROM2,POS2,desc(SR),desc(SR2),desc(DR),desc(DR2))
	df$GENE[is.na(df$GENE)]<-""
	df$GENE2[is.na(df$GENE2)]<-""

	#remove blacklist
	df<-df[Blacklist(df),]

	#remove events not called by two tools.
	df<-df[!is.na(df$CHROM),]
	keeplist<-removeSingles(df)
	df2<-df[keeplist,]

	df2<-Evidence(df2)

	IG<-df2[(grepl("^IG",df2$GENE)&grepl("^IG",df2$GENE2)),]
	df4<-df2[!(grepl("^IG",df2$GENE)&grepl("^IG",df2$GENE2)),]

	#remove medium and low evidence insertions, deletions, inversions and duplications
	df5<-df4[otherSV(df4),]
	df4<-df4[!otherSV(df4),]

	if(nrow(df4) != 0){df4$EVENT<-Event(df4)}else
	{df4$EVENT<-integer()}
	write.table(df4, file=output ,sep="\t", row.names=FALSE)
	
	summarydf<-Summarize(df4)
	write.table(summarydf, file=summary, sep="\t", row.names=FALSE)

	if(nrow(IG) != 0){IG$EVENT<-Event(IG)}
	write.table(IG, file=output2 ,sep="\t", row.names=FALSE)
	if(nrow(df5) != 0){df5$EVENT<-Event(df5)}
	write.table(df5, file=output3 ,sep="\t", row.names=FALSE)

}

merge_SVs(whamfile,gridssfile,novofile,breakmerfile, output,output2,output3,summary)
