#!/usr/bin/env Rscript
#########################################################
#Author: Tjitske de Vries
#date:   20-10-2017
#Name:   reorder-gridss.R

# This script is written to reorder raw gridss output vcf data, add gene information and svtype and to remove chrM events and deduplicate the events.
# Part of the script is based on scripts provided by gridss (github examples)
# this script is part of the translocation snakemake

#call this script by:
# reorder-gridss.R

#########################################################
#               load required libraries                 #
#########################################################
suppressMessages(library(VariantAnnotation))
suppressMessages(library(StructuralVariantAnnotation))
suppressMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(bedr))
suppressMessages(library(plyr))

source("code/reorder_functions.R")

#########################################################
#               Locations annotation files              #
#########################################################

gns <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)

Annotationfile<-snakemake@params[["Annotationfile"]] # contain annotations of areas in genome that are in capturepanel BCNHLv2.
Annotations<-read.delim(Annotationfile, stringsAsFactors = FALSE, sep = "\t" , header=FALSE ) # V1=chr V2=start V3=stop V4=name
Annotations[,1]<-paste("chr",Annotations[,1],sep="")

capture_targets<-get_ext_capture_targets()
# capture_targets<-read.table(snakemake@params[["targets"]], sep="\t", stringsAsFactors =F, header= F) #"code/BCNHLv2_primary_coord.bed"
# capture_targets[,2]<-capture_targets[,2]-300
# capture_targets[,3]<-capture_targets[,3]+300

#########################################################
#      function to annotate SV-types	                #
#########################################################
#copied from gridss github, function to add rearrangement type information.
# simpleEventType <- function(gr) {
#   return(ifelse(seqnames(gr) != seqnames(partner(gr)), "TRL", # inter-chromosomosal
#           ifelse(gr$insLen >= abs(gr$svLen) * 0.7, "INS",
#            ifelse(strand(gr) == strand(partner(gr)), "INV",
#             ifelse(xor(start(gr) < start(partner(gr)), strand(gr) == "-"), "DEL",
#              "DUP")))))
# }
#########################################################
#      functions to annotate genes	                #
#########################################################
#this funtion find genes in UCSC hg19 database and returns genesymbol in vector.
# getGene<-function(gr){
# 	hits <- as.data.frame(findOverlaps(gr, gns, ignore.strand=TRUE))
# 	hits$SYMBOL <- biomaRt::select(org.Hs.eg.db, gns[hits$subjectHits]$gene_id, "SYMBOL")$SYMBOL
# 	hits$gene_strand <- as.character(strand(gns[hits$subjectHits]))
# 	gr$SYMBOL <- ""
# 	gr$SYMBOL[hits$queryHits] <- hits$SYMBOL
# 	return(gr)
# }

#annotate specific breakpoint location of targetgenes (exons, introns, breakpointregion, IG regions)
# annotate_specific<-function(row){
# 	gene<-which(row[1]==Annotations$V1 & as.numeric(row[2]) > Annotations$V2 & as.numeric(row[2]) <= Annotations$V3)
# 	ifelse(length(gene)!=0,Annotations$V4[gene],row[3])
# }

#find position of breakend in the paired event (other row) and add in a column of first event
# getPair<-function(df){
# 	df$CHROM2<-NA_character_
# 	df$POS2<-NA_integer_
# 	df$GENE2<-NA_character_
# 	df$STRAND2<-NA_character_
# 	df[,c("CHROM2","POS2","GENE2","SR2", "RP2","STRAND2")]<-df[df$PARID,c("CHROM","POS","GENE","SR", "RP","STRAND")]
# 	return(df)
# }

# returns false if both breakpoints lay outside capture targets.
# captured<- function(row){
# 	position1<-which(row['CHROM'] == capture_targets$V1 & as.numeric(row['POS']) >=  capture_targets$V2 & as.numeric(row['POS']) <=  capture_targets$V3)
# 	position2<-which(row['CHROM2'] == capture_targets$V1 & as.numeric(row['POS2']) >=  capture_targets$V2 & as.numeric(row['POS2']) <=  capture_targets$V3)
# 	return( length(position1) != 0 | length(position2) != 0 )
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
# 		df[df$CHROM == df$CHROM2 & df$POS > df$POS2, c("CHROM", 'POS', 'CHROM2','POS2', 'GENE','GENE2', 'SR','SR2',"DR","DR2", "BRKPT_COV","BRKPT_COV2","STRAND","STRAND2")]<-df[df$CHROM == df$CHROM2 & df$POS > df$POS2, c('CHROM2','POS2', "CHROM",'POS','GENE2','GENE','SR2','SR',"DR2","DR","BRKPT_COV2", "BRKPT_COV","STRAND2","STRAND")]
# 	return(df)
# }

#########################################################
#      order and filter inputfile	                #
#########################################################

orderGridss <- function(inputfile, output1, output2) {
	vcf <- readVcf(inputfile, "hg19")

	gr <- breakpointRanges(vcf)
	svtype <- simpleEventType(gr)

	info(vcf)$SIMPLE_TYPE <- NA_character_
	info(vcf[gr$vcfId])$SIMPLE_TYPE <- svtype

	info(vcf)$svlen <- NA_integer_
	info(vcf[gr$vcfId])$svlen <- gr$svLen

	# annotate breakends with gene names
	gr<-getGeneGr(gr)

	info(vcf)$GENE <- NA_character_
	info(vcf[gr$vcfId])$GENE <- gr$SYMBOL
	info(vcf)$STRAND  <- NA_character_
	info(vcf[gr$vcfId])$STRAND <- strand(gr)

	#Read VCF file using bedr to create a simple dataframe,
	bedr<-read.vcf(inputfile, split.info = T)
	#remove header information.
	df<-bedr[[2]]

	df<-data.frame(df, info(vcf)[,c("SIMPLE_TYPE", "svlen", "GENE","STRAND")])

	df<-getPair(df)

	#remove chrM - not expected to be real events.
	df<-df[df$CHROM != "chrM" & df$CHROM2 != "chrM",]

	#sort on QUAL
	df<-arrange(df, desc(QUAL))

	# remove the duplicate rows, remove the second of each based on EVENT name (without the o and h in the ID field), because it is filtered on QUAL, the highest QUAL remains.
	df<-df[!duplicated(df$EVENT),]
	df$FILTER[which(is.na(df$FILTER))]<-"PASS"

	# rename some column to get consistent names
	df<- rename(df, replace=c("SVLEN"="oldSVLEN" ,  "SVTYPE"="bnd", "SIMPLE_TYPE"="SVTYPE", "svlen"="SVLEN", "RP"="DR", "RP2"="DR2"))

	df$BRKPT_COV2<-rep(NA_character_, nrow(df))
	df$BRKPT_COV<-rep(NA_character_, nrow(df))

	#remove hits further than 300 bp outside capture areas.
	keeplist<-apply(df ,1, captured)
	df<-df[keeplist,]

	#Annotate genes that are in the capture panel as translocation targets more specific (exons, upstr gene, IGH etc.)
	df$GENE<-apply(df[,c("CHROM","POS","GENE")],1,annotate_specific)
	df$GENE2<-apply(df[,c("CHROM2","POS2","GENE2")],1,annotate_specific)

	#change the order of the columns of the file and remove some columns without (useful) information
	df<-df[!is.na(df$CHROM),]
	df<-orderSvLexo(df)
	df$TOOL<-"gridss"
	df<-df[,c('CHROM', 'POS', 'CHROM2','POS2','GENE','GENE2','SVTYPE', 'SR', 'SR2', 'DR', 'DR2','SVLEN', "BRKPT_COV" ,"BRKPT_COV2" ,"STRAND","STRAND2",'QUAL',"ID","TOOL",'FILTER', 'SRQ','RPQ', "REF", "ALT", "AS",   "ASQ",  "ASRP", "ASSR" ,"BA",   "BAQ" , "BEID", "BQ" ,  "BSC",  "BSCQ", "BUM"  ,"BUMQ" ,"CAS",  "CASQ","CIPOS" , "CIRPOS", "CQ", "HOMLEN"  ,  "HOMSEQ",    "IC",        "IHOMPOS",   "IMPRECISE", "IQ",  "RAS"   ,  "RASQ" ,   "REF.1"   ,"REFPAIR","RSI",'SC', 'SELF')]

	df<-arrange(df,CHROM,POS,CHROM2,POS2,desc(DR),desc(DR2),desc(SR),desc(SR2))

	#remove all duplicate entries
	df2<-df[!duplicated(df[,1:7]),]

	#write dataframe to new file
	write.table(df, file=output1 ,sep="\t", row.names=FALSE)
	write.table(df2, file=output2 ,sep="\t", row.names=FALSE)

}
orderGridss(inputfile=snakemake@input[[1]], output1=snakemake@output[["dups"]],output2=snakemake@output[["ordered"]])


##fileformat=VCFv4.2
##ALT=<ID=INV,Description="Inversion">
##FILTER=<ID=ASSEMBLY_ONLY,Description="Variant is supported only by assembly evidence.">
##FILTER=<ID=ASSEMBLY_TOO_FEW_READ,Description="Not enough reads contribute to this assembly as specified by 'assembly.minReads'">
##FILTER=<ID=ASSEMBLY_TOO_SHORT,Description="This assembly is shorter than a read length">
##FILTER=<ID=LOW_BREAKPOINT_SUPPORT,Description="Does not reach the required threshold quality for calling as specified by 'variantcalling.minScore'">
##FILTER=<ID=LOW_QUAL,Description="Low quality call as specified by 'variantcalling.lowQuality'">
##FILTER=<ID=NO_ASSEMBLY,Description="No assembly supporting this variant could be found.">
##FILTER=<ID=REF,Description="Breakpoint corresponds to reference allele">
##FILTER=<ID=SINGLE_ASSEMBLY,Description="Only one side of the breakpoint could be assembled.">
##FILTER=<ID=SINGLE_SUPPORT,Description="Supported by a single read or read pair only.">
##FILTER=<ID=SMALL_EVENT,Description="Event size is smaller than the minimum reportable size specified by 'variantcalling.minSize'">
##FORMAT=<ID=ASQ,Number=1,Type=Float,Description="Pro-rata quality score contribution of assemblies supporting breakpoint">
##FORMAT=<ID=ASRP,Number=1,Type=Integer,Description="Count of read pairs incorporated into any breakpoint assembly">
##FORMAT=<ID=ASSR,Number=1,Type=Integer,Description="Count of split, soft clipped or indel-containing reads incorporated into any breakpoint assemblies">
##FORMAT=<ID=BAQ,Number=1,Type=Float,Description="Pro-rata quality score contribution of assemblies supporting just local breakend">
##FORMAT=<ID=BQ,Number=1,Type=Float,Description="Quality score of breakend evidence after evidence reallocation">
##FORMAT=<ID=BSC,Number=1,Type=Integer,Description="Count of soft clips supporting just local breakend per category">
##FORMAT=<ID=BSCQ,Number=1,Type=Float,Description="Quality score of soft clips supporting just local breakend per category">
##FORMAT=<ID=BUM,Number=1,Type=Integer,Description="Count of read pairs (with one read unmapped) supporting just local breakend per category">
##FORMAT=<ID=BUMQ,Number=1,Type=Float,Description="Quality score of read pairs (with one read unmapped) supporting just local breakend per category">
##FORMAT=<ID=CASQ,Number=1,Type=Float,Description="Pro-rata quality score of complex compound breakpoint assemblies supporting breakpoint from elsewhere">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=IC,Number=1,Type=Integer,Description="Count of read indels supporting breakpoint per category">
##FORMAT=<ID=IQ,Number=1,Type=Float,Description="Quality score of read indels supporting breakpoint per category">
##FORMAT=<ID=QUAL,Number=1,Type=Float,Description="Quality score of breakend evidence after evidence reallocation">
##FORMAT=<ID=RASQ,Number=1,Type=Float,Description="Pro-rata quality score contribution of assemblies supporting breakpoint from remote breakend">
##FORMAT=<ID=REF,Number=1,Type=Integer,Description="Count of reads mapping across this breakend">
##FORMAT=<ID=REFPAIR,Number=1,Type=Integer,Description="Count of reference read pairs spanning this breakpoint supporting the reference allele">
##FORMAT=<ID=RP,Number=1,Type=Integer,Description="Count of read pairs supporting breakpoint per category">
##FORMAT=<ID=RPQ,Number=1,Type=Float,Description="Quality score of read pairs supporting breakpoint per category">
##FORMAT=<ID=SR,Number=1,Type=Integer,Description="Count of split reads supporting breakpoint per category">
##FORMAT=<ID=SRQ,Number=1,Type=Float,Description="Quality score of split reads supporting breakpoint per category">
##INFO=<ID=AS,Number=1,Type=Integer,Description="Count of assemblies supporting breakpoint">
##INFO=<ID=ASQ,Number=1,Type=Float,Description="Quality score of assemblies supporting breakpoint">
##INFO=<ID=ASRP,Number=1,Type=Integer,Description="Count of read pairs incorporated into any breakpoint assembly">
##INFO=<ID=ASSR,Number=1,Type=Integer,Description="Count of split, soft clipped or indel-containing reads incorporated into any breakpoint assemblies">
##INFO=<ID=BA,Number=1,Type=Integer,Description="Count of assemblies supporting just local breakend">
##INFO=<ID=BAQ,Number=1,Type=Float,Description="Quality score of assemblies supporting just local breakend">
##INFO=<ID=BEID,Number=.,Type=String,Description="Breakend assemblies contributing support to the breakpoint.">
##INFO=<ID=BQ,Number=1,Type=Float,Description="Quality score of breakend evidence">
##INFO=<ID=BSC,Number=1,Type=Integer,Description="Count of soft clips supporting just local breakend">
##INFO=<ID=BSCQ,Number=1,Type=Float,Description="Quality score of soft clips supporting just local breakend">
##INFO=<ID=BUM,Number=1,Type=Integer,Description="Count of read pairs (with one read unmapped) supporting just local breakend">
##INFO=<ID=BUMQ,Number=1,Type=Float,Description="Quality score of read pairs (with one read unmapped) supporting just local breakend">
##INFO=<ID=CAS,Number=1,Type=Integer,Description="Count of complex compound breakpoint assemblies supporting breakpoint from elsewhere">
##INFO=<ID=CASQ,Number=1,Type=Float,Description="Quality score of complex compound breakpoint assemblies supporting breakpoint from elsewhere">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=CIRPOS,Number=2,Type=Integer,Description="Confidence interval around remote breakend POS for imprecise variants">
##INFO=<ID=CQ,Number=1,Type=Float,Description="Breakpoint quality score before evidence reallocation">
##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
##INFO=<ID=EVENT,Number=1,Type=String,Description="ID of event associated to breakend">
##INFO=<ID=HOMLEN,Number=.,Type=Integer,Description="Length of base pair identical micro-homology at event breakpoints">
##INFO=<ID=HOMSEQ,Number=.,Type=String,Description="Sequence of base pair identical micro-homology at event breakpoints">
##INFO=<ID=IC,Number=1,Type=Integer,Description="Count of read indels supporting breakpoint">
##INFO=<ID=IHOMPOS,Number=2,Type=Integer,Description="Position of inexact homology">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=IQ,Number=1,Type=Float,Description="Quality score of read indels supporting breakpoint">
##INFO=<ID=PARID,Number=1,Type=String,Description="ID of partner breakend">
##INFO=<ID=RAS,Number=1,Type=Integer,Description="Count of assemblies supporting breakpoint from remote breakend">
##INFO=<ID=RASQ,Number=1,Type=Float,Description="Quality score of assemblies supporting breakpoint from remote breakend">
##INFO=<ID=REF,Number=1,Type=Integer,Description="Count of reads mapping across this breakend">
##INFO=<ID=REFPAIR,Number=1,Type=Integer,Description="Count of reference read pairs spanning this breakpoint supporting the reference allele">
##INFO=<ID=RP,Number=1,Type=Integer,Description="Count of read pairs supporting breakpoint">
##INFO=<ID=RPQ,Number=1,Type=Float,Description="Quality score of read pairs supporting breakpoint">
##INFO=<ID=RSI,Number=.,Type=Integer,Description="Support interval offsets of partner breakend.">
##INFO=<ID=SC,Number=1,Type=String,Description="CIGAR for displaying anchoring alignment of any contributing evidence and microhomologies.">
##INFO=<ID=SELF,Number=0,Type=Flag,Description="Indicates a breakpoint is self-intersecting">
##INFO=<ID=SI,Number=.,Type=Integer,Description="Support interval offsets from breakend position in which at least one supporting read/read pair/assembly is mapped.">
##INFO=<ID=SR,Number=1,Type=Integer,Description="Count of split reads supporting breakpoint">
##INFO=<ID=SRQ,Number=1,Type=Float,Description="Quality score of split reads supporting breakpoint">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
