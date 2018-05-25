#!/usr/bin/env R
#########################################################

#Author: Tjitske de Vries
#date:   08-11-2017
#Name:   merge_functions.R

#this script contain the functions used for filtering and merging
# the output of Gridss, Breakmer, Wham and novobreak.

 #########################################################

 #########################################################
 #   function to remove lines without multiple matches   #
 #########################################################
 #Remove hits which are found with one tool only.
 removeSingles<-function(df){
 	keeplist<-c()
 	for (i in 1:nrow(df)){
 		range=max(0,i-50):min(i+50,nrow(df))

 		TF<-df[i,"CHROM"]== df[range,"CHROM"] &  abs(df[i,"POS"] - df[range,"POS"])< 10 &  df[i,"CHROM2"]== df[range,"CHROM2"] &  abs(df[i,"POS2"] - df[range,"POS2"])< 10 & df[i,"TOOL"] != df[range,"TOOL"]

 		keeplist<-c(keeplist,(sum(TF)>0))
 	}
 	return(keeplist)
 }

 #count number of events, an event is a single structural variant, a reciprocal translocation can consist of one or two events depending on the distance of the two breakpoints from each other. within 10 bp is called 1 event.
 Event<-function(df){
  df$EVENT<-NA
 	event<-0
 	if(nrow(df)>1){
 		for(i in 1:(nrow(df))){
      same_event<-which(df[i,"CHROM"]== df["CHROM"] &  abs(df[i,"POS"] - df["POS"])< 10 &  df[i,"CHROM2"]== df["CHROM2"] &  abs(df[i,"POS2"] - df["POS2"])< 10)
        if(sum(!is.na(df$EVENT[same_event]))==0){
          event<-event+1
          df$EVENT[same_event]<-event
        }else{df$EVENT[same_event]<-event}
      }
    }else{df$EVENT<-1}
  return(df$EVENT)
  }

RemoveSingleToolEvents<-function(df){
allevents<-unique(df$EVENT)
newdf<-df
  for(i in allevents){
    if(length(unique(df[df$EVENT==i,"TOOL"]))==1){
      newdf<-newdf[newdf$EVENT!=i,]
    }
  }
  return(newdf)
}

CalculateVAF<-function(df){
  df$VAF<-round(df$SR/df$BRKPT_COV*100, digits=1)
  df$VAF2<-round(df$SR2/df$BRKPT_COV2*100, digits=1)
return(df)
}

 # Annotate if there is a hign level of evidence for a translocation based on discordant and split reads.
 Evidence<-function(df, highSVevidence){
 	df$EVIDENCE_LEVEL[df$SR >= highSVevidence & df$SR2 >= highSVevidence & (df$DR >=highSVevidence | is.na(df$DR)) &  (df$DR2 >=highSVevidence | is.na(df$DR2))]<-"HIGH"
 	df$EVIDENCE_LEVEL[!(df$SR >= highSVevidence & df$SR2 >= highSVevidence & (df$DR >=highSVevidence | is.na(df$DR)) &  (df$DR2 >=highSVevidence | is.na(df$DR2)))]<-"LOW"

 	return(df)
 }

 Create_bed<-function(df){
   bed<-data.frame(chr=df$CHROM, start=df$POS-2, stop=df$POS+2)
   bed2<-data.frame(chr=df$CHROM2, start=df$POS2-2, stop=df$POS2+2)
   combinedbed<-rbind(bed, bed2)
   return(combinedbed)
 }

 #########################################################
 #               Filter functions for the 4 tools        #
 #########################################################
 # Remove breakmer hits not meeting criteria
 filter_breakmer<-function(df, dr , sr ){

 	df_filt<-df[ df$SR >= sr & df$SR2 >= sr & (df$DR >= dr | is.na(df$DR)) & (df$DR2 >= dr | is.na(df$DR2)),]
 	return(df_filt)
 }

 # Remove wham hits not meeting criteria
 filter_wham<-function(df, sr ){

 	sumDR_SR<-rowSums(df[c("SR","SR2","DR","DR2")], na.rm = T)
 	df_filt<-df[sumDR_SR>= 2*sr,]
 	return(df_filt)
 }

 # Remove gridss hits not meeting criteria
 filter_gridss<-function(df, gridssscore ){

 	df_filt<-df[df$QUAL >=gridssscore ,]
 	return(df_filt)
 }

 #Remove novobreak hits not meeting criteria
 filter_novobreak<-function(df,novoscore){

 	df_filt<-df[df$QUAL >= novoscore ,]
 	return(df_filt)
 }

 #Remove hits of which one of the breakpoints is located in blacklisted region.
Blacklist<- function(row){
        position1<-which(row['CHROM'] == blacklist$V1 & as.numeric(row['POS']) >=  blacklist$V2 & as.numeric(row['POS']) <=  blacklist$V3)
        position2<-which(row['CHROM2'] == blacklist$V1 & as.numeric(row['POS2']) >=  blacklist$V2 & as.numeric(row['POS2']) <=  blacklist$V3)
        return( length(position1) != 0 | length(position2) != 0 )
}

 #function to find out if event is translocation or High evidence other event
 otherSV<-function(df){
 	!(df$SVTYPE == "TRL")
 }

 tra_otherSV_HIGH<-function(df){
   (df$SVTYPE == "TRL" | df$EVIDENCE_LEVEL == "HIGH")
 }

 #function to keep only one line per event with highest level of evidence
 Summarize<-function(df){
 	df<-arrange(df,EVENT, EVIDENCE_LEVEL,plyr::desc(DR),plyr::desc(DR2),plyr::desc(SR),plyr::desc(SR2))
  allevents<-unique(df$EVENT)
  for(i in allevents){
    df$ALL_TOOLS[df$EVENT==i]<-paste(unique(df[df$EVENT==i,"TOOL"]),collapse=",")
    df$OCCURENCE[df$EVENT==i]<-sum(df$EVENT==i)
  }
 	if(nrow(df)>0){df$TPorFP<-"TP"}else{df$TPorFP<-character()} #this column is used to exclude rows from circlize plot if they are marked as FP
 	summarydf<-df[!duplicated(df$EVENT),]
 }
