#!/usr/bin/env Rscript
#########################################################
#Author: Tjitske de Vries
#date:   5-12-2017
#Name:   Calculate_brkpt_freq.R

#########################################################

coverage_file<-snakemake@input[["cov"]]
summary<-snakemake@input[["summary"]]

output<-snakemake@output[["brkpt_freq"]]

Calculate_brkpt_freq<-function(coverage_file, summary, output){

coverage<-tryCatch(read.delim(coverage_file, header = F, stringsAsFactors = F), error=function(e) NULL)
summarydf<-read.delim(summary, stringsAsFactors = F)
if(!is.null(coverage)){
  summarydf[,c("COV_LOW", "COV_HIGH","CALCULATED_BRKPT_FREQ", "COV_LOW2", "COV_HIGH2","CALCULATED_BRKPT_FREQ2")]<-NA
    for( i in seq(1,nrow(coverage),by=4)){
    COV_HIGH<-max(coverage[i:(i+3),"V5"])
    COV_LOW<-min(coverage[i:(i+3),"V5"])
    difference<-COV_HIGH-COV_LOW
    CALCULATED_BRKPT_FREQ<-round(difference/COV_HIGH*100, digits=1)
    if(i <= nrow(coverage)/2 ){
      summarydf[(((i-1)/4)+1),c("COV_LOW", "COV_HIGH","CALCULATED_BRKPT_FREQ")]<-c(COV_LOW, COV_HIGH,CALCULATED_BRKPT_FREQ)
    }else{ summarydf[(((i-nrow(coverage)/2-1)/4)+1),c("COV_LOW2", "COV_HIGH2","CALCULATED_BRKPT_FREQ2")]<-c(COV_LOW, COV_HIGH,CALCULATED_BRKPT_FREQ)
    }
  }
}else{summarydf<-data.frame(summarydf,"COV_LOW"=integer(), "COV_HIGH"=integer(),"CALCULATED_BRKPT_FREQ"=numeric(), "COV_LOW2"=integer(), "COV_HIGH2"=integer(),"CALCULATED_BRKPT_FREQ2"=numeric())}
write.table(summarydf,file=output, row.names=FALSE, sep="\t")
}
Calculate_brkpt_freq(coverage_file, summary, output)
