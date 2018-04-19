#!/usr/bin/env Rscript
#########################################################

#Author: Tjitske de Vries
#date:   07-11-2017
#Name:   reorder-wham-snake.R

# This script is written to reorder wham output data with the translocation-snake pipeline.

# reorder-wham-snake.R

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
#            reordering wham		                #
#########################################################
orderWham<-function(inputfile, output1, output2){

	bedr<-read.vcf(inputfile, split.info = T)
	FORMAT<-unlist(strsplit(bedr$vcf[1,"FORMAT"],":"))

	vcf<-bedr$vcf
	vcf<-cSplit_f(vcf, sep=":", splitCols = 28)
	colnames(vcf)[28:33]<-FORMAT

	#Create unique ID per hit
	vcf$ID<-paste("wham",seq(1:nrow(vcf)), sep="-")

	#remove chrM - not expected to be real events.
	vcf<-vcf[vcf$CHROM !="M" & vcf$CHR2 != "M",]

	vcf$CHROM<-paste("chr",vcf[,CHROM], sep="")
	vcf$CHR2<-paste("chr",vcf[,CHR2], sep="")

	vcf<-cSplit(vcf, sep=",", splitCols ="SP", direction = "wide") # the SP column contain "Number of reads supporting endpoint: mate-position,split-read,alternative-mapping"

	vcf$SVTYPE<-apply(vcf,1,function(x){ifelse(x["CHROM"] != x["CHR2"],"TRL","other")})
	vcf$QUAL<-vcf$MQ

	vcf<-plyr::rename(vcf, replace=c("CHR2"="CHROM2", "END"="POS2", "NS"="SR","SP_1"="DR", "SP_2"="SR2", "RD"="BRKPT_COV" ))

	#remove hits further than 300 bp outside capture areas.
	keeplist<-apply(vcf ,1, captured)
	vcf<-vcf[keeplist,]

	#Annotate genes
	if(nrow(vcf)!=0){
	  vcf$GENE<-getGeneBed(vcf[,"CHROM"], vcf[,"POS"])
	  vcf$GENE2<-getGeneBed(vcf[,"CHROM2"],vcf[, "POS2"])

	  #Annotate genes that are in the capture panel as translocation targets more specific (exons, upstr gene, IGH etc.)
	  vcf$GENE<-apply(vcf[,c("CHROM","POS","GENE")],1,annotate_specific)
	  vcf$GENE2<-apply(vcf[,c("CHROM2","POS2","GENE2")],1,annotate_specific)
	}else{vcf$GENE<-rep(NA_character_, nrow(vcf))
	vcf$GENE2<-rep(NA_character_, nrow(vcf))}
	  #create required columns in correct format.
	vcf$DR2<-as.integer(rep(NA_character_, nrow(vcf)))
	vcf$STRAND<-rep(NA_character_, nrow(vcf))
	vcf$STRAND2<-rep(NA_character_, nrow(vcf))
	vcf$BRKPT_COV2<-as.integer(rep(NA_character_, nrow(vcf)))
	vcf$POS2<-as.integer(vcf$POS2)
	vcf$DR<-as.integer(as.character(vcf$DR))
	vcf$TOOL<-"wham"
	vcf<-orderSvLexo(vcf)

	vcf<-vcf[,c('CHROM', 'POS', 'CHROM2','POS2','GENE','GENE2','SVTYPE','SR', 'SR2',"DR", "DR2" ,"SVLEN", "BRKPT_COV" ,"BRKPT_COV2" ,"STRAND","STRAND2",'QUAL',"ID","TOOL", "PU" , "SU" ,"CU" ,"NC" ,"MQF", "DI", "NR", "NA", "SP_3" ,"GT","GL", "CISTART" ,"CIEND"  ,"CF" ,"ALT","LRT","WAF","GC","AT")]
	#excluded: "REF", "FILTER" , "MQ" ,"FORMAT","DP" , because these are non infomative

	vcf<-arrange(vcf,CHROM,POS,CHROM2,POS2,plyr::desc(DR),plyr::desc(DR2),plyr::desc(SR),plyr::desc(SR2))

	#remove all duplicate entries
	vcf2<-vcf[!duplicated(vcf[,1:7]),]

	write.table(vcf, file=output1,row.names=FALSE, sep="\t")
	write.table(vcf2, file=output2,row.names=FALSE, sep="\t")
}

orderWham(inputfile=snakemake@input[["vcf"]], output1=snakemake@output[["dups"]],output2=snakemake@output[["ordered"]])


#$header$INFO
 #     ID        Number Type        Description
# [1,] "LRT"     "1"    "Float"     "Likelihood Ratio Test statistic"
# [2,] "WAF"     "3"    "Float"     "Allele frequency of: background,target,combined"
# [3,] "GC"      "2"    "Integer"   "Number of called genotypes in: background,target"
# [4,] "AT"      "15"   "Float"     "Pileup attributes"
# [5,] "CF"      "1"    "Float"     "Fraction of reads with more than three cigar operations"
# [6,] "CISTART" "2"    "Integer"   "PE confidence interval around POS"
# [7,] "CIEND"   "2"    "Integer"   "PE confidence interval around END"
# [8,] "PU"      "1"    "Integer"   "Number of reads supporting position"
# [9,] "SU"      "1"    "Integer"   "Number of supplemental reads supporting position"
#[10,] "CU"      "1"    "Integer"   "Number of neighboring all soft clip clusters across all individuals at pileup position"
#[11,] "DP"      "1"    "Integer"   "Number of reads at pileup position across individuals passing filters"
#[12,] "NC"      "1"    "Float"     "Number of soft clipped sequences collapsed into consensus"
# [13,] "MQ"/"QUAL"      "1"    "Float"     "Average mapping quality"      -
#[14,] "MQF"     "1"    "Float"     "Fraction of reads with MQ less than 50"
#[15,] "SP"      "3"    "Integer"   "Number of reads supporting endpoint: mate-position,split-read,alternative-mapping"  #splitted up in RP, and SR2 and SP_3
#[16,] "CHR2"/"CHROM2    "3"    "String"    "Other seqid"
#[17,] "DI"      "1"    "Character" "Consensus is from front or back of pileup: f,b"
#[18,] "END"/"POS2"     "1"    "Integer"   "End position of the variant described in this record"
#[19,] "SVLEN"   "1"    "Integer"   "Difference in length between POS and END"
#[20,] "WC" /SVTYPE     "1"    "String"    "WHAM classifier variant type"
#[21,] "WP"      "4"    "Float"     "WHAM probability estimate for each structural variant classification from RandomForest model"

#$header$FORMAT
 #    ID   Number Type      Description
#[1,] "GT" "1"    "String"  "Genotype"
#[2,] "GL" "A"    "Float"   "Genotype Likelihood"
#[3,] "NR" "1"    "Integer" "Number of reads that do not support a SV"
#[4,] "NA" "1"    "Integer" "Number of reads supporting a SV"
#[5,] "NS"/"SR" "1"    "Integer" "Number of reads with a softclip at POS for individual"
#[6,] "RD"/BRKPT_COV "1"    "Integer" "Number of reads passing filters"
