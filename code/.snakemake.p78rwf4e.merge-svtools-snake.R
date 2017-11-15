
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
    input = list('../wham/TVU12-11384FF_S28_L004/TVU12-11384FF_S28_L004-wham_ordered.csv', '../novobreak/TVU12-11384FF_S28_L004/TVU12-11384FF_S28_L004-novobreak_ordered.csv', '../gridss/TVU12-11384FF_S28_L004/TVU12-11384FF_S28_L004-gridss_ordered.csv', '../breakmer/TVU12-11384FF_S28_L004/TVU12-11384FF_S28_L004-breakmer_ordered.csv', "wham" = '../wham/TVU12-11384FF_S28_L004/TVU12-11384FF_S28_L004-wham_ordered.csv', "novobreak" = '../novobreak/TVU12-11384FF_S28_L004/TVU12-11384FF_S28_L004-novobreak_ordered.csv', "gridss" = '../gridss/TVU12-11384FF_S28_L004/TVU12-11384FF_S28_L004-gridss_ordered.csv', "breakmer" = '../breakmer/TVU12-11384FF_S28_L004/TVU12-11384FF_S28_L004-breakmer_ordered.csv'),
    output = list('../merged/TVU12-11384FF_S28_L004-other_merged.csv', '../merged/TVU12-11384FF_S28_L004-trl_merged.csv', '../merged/TVU12-11384FF_S28_L004-IG_merged.csv', '../merged/TVU12-11384FF_S28_L004-trl_summary.csv', "IG" = '../merged/TVU12-11384FF_S28_L004-IG_merged.csv', "trl" = '../merged/TVU12-11384FF_S28_L004-trl_merged.csv', "otherSVs" = '../merged/TVU12-11384FF_S28_L004-other_merged.csv', "summary" = '../merged/TVU12-11384FF_S28_L004-trl_summary.csv'),
    params = list(),
    wildcards = list('TVU12-11384FF_S28_L004', "sample" = 'TVU12-11384FF_S28_L004'),
    threads = 1,
    log = list(),
    resources = list(),
    config = list("filters" = list("gridssscore" = 450, "discread" = 3, "highSVevidence" = 15, "splitread" = 4, "novoscore" = 4, "blacklist" = 'code/blacklist.bed'), "all" = list("Annotationfile" = 'code/annotation_new_BCNHL-v2.bed', "targets" = 'code/BCNHLv2_primary_coord.bed', "THREADS" = 4, "REF_NOCHR" = '/net/nfs/PAT/home/tjitske/files/ref/excl_chr/hg19_refgen_exclchr.fa', "REF_CHR" = '/net/nfs/PAT/home/tjitske/files/ref/incl_chr/hg19.fa'), "gridss" = list("GRIDSS_JAR" = '/net/nfs/PAT/home/tjitske/GRIDSS/gridss-1.4.2-jar-with-dependencies.jar'), "novobreak" = list("NORMAL" = '/net/nfs/PAT/analysis/MPS-310/Blood_combined_bam/combined_healthy_blood_nochr.bam', "HEADER" = '/net/nfs/PAT/home/tjitske/novobreak/nb_distribution/header.txt', "NOVOBREAK" = '/net/nfs/PAT/home/tjitske/novobreak/nb_distribution/run_novoBreak.sh', "EXE_DIR" = '/net/nfs/PAT/home/tjitske/novobreak/nb_distribution'), "wham" = list("MAPQUAL" = 10, "WHAM" = '/net/nfs/PAT/home/tjitske/WHAM/wham/bin/wham', "CLASSIFY" = '/net/nfs/PAT/home/tjitske/WHAM/wham/utils/classify_WHAM_vcf.py', "TRAIN" = '/net/nfs/PAT/home/tjitske/WHAM/wham/data/WHAM_training_data.txt', "BASEQUAL" = 5)),
    rule = 'merge_svtools'
)
######## Original script #########
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

	if(nrow(df4) != 0){df4$EVENT<-Event(df4)}
	write.table(df4, file=output ,sep="\t", row.names=FALSE)
	
	summarydf<-Summarize(df4)
	write.table(summarydf, file=summary, sep="\t", row.names=FALSE)

	if(nrow(IG) != 0){IG$EVENT<-Event(IG)}
	write.table(IG, file=output2 ,sep="\t", row.names=FALSE)
	if(nrow(df5) != 0){df5$EVENT<-Event(df5)}
	write.table(df5, file=output3 ,sep="\t", row.names=FALSE)

}

merge_SVs(whamfile,gridssfile,novofile,breakmerfile, output,output2,output3,summary)
