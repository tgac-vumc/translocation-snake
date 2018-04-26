#!/bin/bash

export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8
export LANGUAGE=en_US.UTF-8
export LC_NUMERIC=en_US.UTF-8

title=$1
samples=$2


echo -e '<?xml version="1.0" encoding="UTF-8"?>'
echo -e '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.1//EN"'
echo -e '\t"http://www.w3.org/TR/xhtml11/DTD/xhtml11.dtd">'
echo -e '<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en">'
echo -e '\t<head>'
echo -e "\t\t<title>$title</title>"
echo -e '\t\t<style type="text/css">'
echo -e '\t\t\ttable {border-collapse: collapse; font-size: smaller;}'
echo -e '\t\t\tth {border: 1px solid gray; padding: 5px;}'
echo -e '\t\t\ttd {border: 1px solid gray; padding: 5px; text-align: right;}'
echo -e '\t\t</style>'
#echo -e '\t\t<link rel="stylesheet" href="http://ccagc-gen01.vumc.nl/js/lightbox2-master/dist/css/lightbox.css">'
echo -e '\t\t<link rel="stylesheet" href="../translocation-snake/lb2/css/lightbox.css">'
echo -e '\t</head>'
echo -e '\t<body>'
echo -e "\t\t<h1>$title</h1>"
#echo -e '\t\t<p><a href="../index.html">Back to overview</a></p>'
echo -e '\t\t<table>'
echo -e '\t\t\t<tr>'
echo -e '\t\t\t\t<th>sample</th>'
echo -e '\t\t\t\t<th>total</th>'
echo -e '\t\t\t\t<th colspan="2">unique</th>'
echo -e '\t\t\t\t<th colspan="2">aligned</th>'
echo -e '\t\t\t\t<th>mean target cov</th>'
echo -e '\t\t\t\t<th>%30x</th>'
echo -e '\t\t\t\t<th>%100x</th>'
echo -e '\t\t\t\t<th>%chimeric reads</th>'
echo -e '\t\t\t\t<th>qc</th>'
echo -e '\t\t\t\t<th><a href="circlize/index.html">circlize</a></th>'
echo -e '\t\t\t\t<th><a href="../merged/">translocations</a></th>'
echo -e '\t\t\t\t<th><a href="../merged/">other SVs</a></th>'
echo -e '\t\t\t\t<th><a href="../merged/">IG events</a></th>'
#echo -e '\t\t\t\t<th><a href="profiles/reCalled/index.html">reCalled</a></th>'
echo -e '\t\t\t</tr>'

sumtotal=0
sumaligned=0
sumunique=0
summtc=0
sum30x=0
sum100x=0
sumchimeric=0
sumtrl=0
sumoth=0
sumig=0
n=0

for sample in $samples
do
  n=`echo $n + 1 | bc -l`
  sample=`basename $sample _coordsorted.bam`
  sample3=`echo ../fastqc/${sample}*_fastqc.html`
  sample2=`echo $sample | sed -r 's/^[0-9]{6}_[A-Z0-9]{9}_L[1-8]{1}_//'`   #TOTO - remove line and remove sample lines
  total=`awk 'NR==8 { print $6; }' ../CovMetrics/${sample}_HSmetrics.txt`
  unique=`awk 'NR==8 { print $8 }' ../CovMetrics/${sample}_HSmetrics.txt`
  percunique=`awk 'NR==8 { print $10 }' ../CovMetrics/${sample}_HSmetrics.txt`
  aligned=`awk 'NR==8 { print $11; }' ../CovMetrics/${sample}_HSmetrics.txt`
  percaligned=`awk 'NR==8 { print $12; }' ../CovMetrics/${sample}_HSmetrics.txt`
  mtc=`awk 'NR==8 { print $23; }' ../CovMetrics/${sample}_HSmetrics.txt`
  x30=`awk 'NR==8 { print $39; }' ../CovMetrics/${sample}_HSmetrics.txt`
  x100=`awk 'NR==8 { print $42; }' ../CovMetrics/${sample}_HSmetrics.txt`
  chimeric=`awk 'NR==10 { print $21; }' ../gridss/${sample}/${sample}_coordsorted.bam.gridss.working/${sample}_coordsorted.bam.alignment_summary_metrics`

  numtrl=`wc -l < "../merged/$sample-trl_summary_brkpt_freq.csv"`
  numoth=`wc -l < "../merged/$sample-other_merged.csv"`
  numig=`wc -l < "../merged/$sample-IG_merged.csv"`

  sumtotal=`echo $sumtotal + $total | bc -l`
  sumaligned=`echo $sumaligned + $aligned | bc -l`
  sumunique=`echo $sumunique + $unique | bc -l`
  summtc=`echo $summtc + $mtc | bc -l`
  sum30x=`echo $sum30x + $x30 | bc -l`
  sum100x=`echo $sum100x + $x100 | bc -l`
  sumchimeric=`echo $sumchimeric + $chimeric | bc -l`
if [ -f "../merged/$sample-trl_summary_brkpt_freq.csv" ]
then
  sumtrl=`echo $sumtrl + $numtrl | bc -l`
  sumoth=`echo $sumoth + $numoth | bc -l`
  sumig=`echo $sumig + $numig | bc -l`
fi
  echo -e '\t\t\t<tr>'
  echo -e '\t\t\t\t<td style="text-align: left;">'$sample2'</td>'
  printf "\t\t\t\t<td>%'i</td>\n" $total
  printf "\t\t\t\t<td>%'i</td>\n" $unique
  printf "\t\t\t\t<td>%.2f%%</td>\n" `echo $percunique*100 | bc -l`
  printf "\t\t\t\t<td>%'i</td>\n" $aligned
  printf "\t\t\t\t<td>%.2f%%</td>\n" `echo $percaligned*100 | bc -l`
  printf "\t\t\t\t<td>%.0f</td>\n" $mtc
  printf "\t\t\t\t<td>%.2f%%</td>\n" `echo $x30*100 | bc -l`
  printf "\t\t\t\t<td>%.2f%%</td>\n" `echo $x100*100 | bc -l`
  printf "\t\t\t\t<td>%.2f%%</td>\n" `echo $chimeric*100 | bc -l`
  echo -en '\t\t\t\t<td><a href='$sample3'>fastq</a>, '
  echo -e '\t\t\t\t<td><a href="circlize/'$sample'-circlize.png">circlize</a></td>'
  echo -e '\t\t\t\t<td><a href="../merged/'$sample'-trl_summary_brkpt_freq.csv">'`expr $numtrl - 1`'</a></td>'
  echo -e '\t\t\t\t<td><a href="../merged/'$sample'-other_merged.csv">'`expr $numoth - 1`'</a></td>'
  echo -e '\t\t\t\t<td><a href="../merged/'$sample'-IG_merged.csv"> '`expr $numig - 1`'</a></td>'
  echo -e '\t\t\t</tr>'
done
echo -e '\t\t\t<tr style="border-top: double;">'
echo -e '\t\t\t\t<td style="text-align: left;">Average</td>'
printf "\t\t\t\t<td>%.0f</td>\n" `echo $sumtotal/$n | bc -l`
printf "\t\t\t\t<td>%.0f</td>\n" `echo $sumunique/$n| bc -l`
printf "\t\t\t\t<td>%.2f%%</td>\n" `echo $sumunique/$sumtotal*100 | bc -l`
printf "\t\t\t\t<td>%.0f</td>\n" `echo $sumaligned/$n | bc -l`
printf "\t\t\t\t<td>%.2f%%</td>\n" `echo $sumaligned/$sumunique*100 | bc -l`
printf "\t\t\t\t<td>%.0f</td>\n" `echo $summtc/$n | bc -l`
printf "\t\t\t\t<td>%.2f%%</td>\n" `echo $sum30x/$n*100 | bc -l`
printf "\t\t\t\t<td>%.2f%%</td>\n" `echo $sum100x/$n*100 | bc -l`
printf "\t\t\t\t<td>%.2f%%</td>\n" `echo $sumchimeric/$n*100 | bc -l`
echo -e '\t\t\t\t<td>&nbsp;</td>'
echo -e '\t\t\t\t<td>&nbsp;</td>'
printf "\t\t\t\t<td>%.2f</td>\n" `echo "($sumtrl-$n)/$n" | bc -l`
printf "\t\t\t\t<td>%.2f</td>\n" `echo "($sumoth-$n)/$n"  | bc -l`
printf "\t\t\t\t<td>%.2f</td>\n" `echo "($sumig-$n)/$n" | bc -l`
echo -e '\t\t\t</tr>'
#echo -e '\t\t\t\t<td>&nbsp;</td>'
#echo -e '\t\t\t</tr>'

echo -e '\t\t</table>'
#echo -e '\t<script src="http://ccagc-gen01.vumc.nl/js/lightbox2-master/dist/js/lightbox-plus-jquery.min.js"></script>'
echo -e '\t<script src=""../translocation-snake/lb2/js/lightbox-plus-jquery.min.js"></script>'
echo -e '\t</body>'
echo -e '</html>'
