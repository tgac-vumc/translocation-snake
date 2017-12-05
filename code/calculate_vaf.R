#!/usr/bin/env Rscript
#########################################################
#Author: Tjitske de Vries
#date:   5-12-2017
#Name:   calculate_vaf.R

#########################################################

coverage_file<-snakemake@input[["cov"]]
summary<-snakemake@input[["summary"]]

output<-snakemake@output[["vaf"]]

calculate_vaf<-function(coverage_file, summary, output){

coverage<-read.delim(coverage_file, header = F, stringsAsFactors = F)
summarydf<-read.delim(summary, stringsAsFactors = F)
summarydf[,c("COV_LOW", "COV_HIGH","CALCULATED_VAF", "COV_LOW2", "COV_HIGH2","CALCULATED_VAF2")]<-NA
  for( i in seq(1,nrow(coverage),by=4)){
    COV_HIGH<-max(coverage[i:(i+3),"V5"])
    COV_LOW<-min(coverage[i:(i+3),"V5"])
    difference<-COV_HIGH-COV_LOW
    CALCULATED_VAF<-round(difference/COV_HIGH*100, digits=1)
    if(i <= nrow(coverage)/2 ){
      summarydf[(((i-1)/4)+1),c("COV_LOW", "COV_HIGH","CALCULATED_VAF")]<-c(COV_LOW, COV_HIGH,CALCULATED_VAF)
    }else{ summarydf[(((i-1)/4)-3),c("COV_LOW2", "COV_HIGH2","CALCULATED_VAF2")]<-c(COV_LOW, COV_HIGH,CALCULATED_VAF)
    }
  }
write.table(summarydf,file=output, row.names=FALSE, sep="\t")
}
calculate_vaf(coverage_file, summary, output)
