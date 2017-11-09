# Translocation detection pipeline by next generation sequencing
## a combined approach of BreaKmer, Gridss, Wham and novoBreak.

The translocation detection pipeline is designed to increase the sensitivity and specificity of translocation detection by using a combined approach of multiple publically available bioinformatic tools. This is especially important for the analysis of DNA derived from formalin fixed paraffin embedded (FFPE) material, because this material is in general of less quality compared to fresh material and contain more noise. The analysis is

**requirements**
- Gridss v1.4.2  https://github.com/PapenfussLab/gridss/releases
- WHAM https://github.com/zeeev/wham
- BreaKmer v0.0.4 https://github.com/ryanabo/BreaKmer
- Novobreak v1.1.3   https://github.com/czc/nb_distribution
    For novobreak replace infer_bp_v4.pl in the novobreak installation folder with the infer_bp_v4.pl file which you can find in the bin folder.
- Snakemake http://snakemake.readthedocs.io/en/stable/  (follow tuturial setup steps)
- The other requirements can be installed using the environment.yaml file
    [conda env create --name translocation-snake --file environment.yaml]


**input file creation**

Pipeline is tested on paired-end ilumina data and processed as follow:
- Adapters are removed with SeqPurge v0.1-104 (-min_len 20)
- DNA is aligned with bwa_mem v0.7.12-r1039 against reference hg19. (-M -R <readgroupinfo> )
- Realignment around indels is performed with ABRA v0.96
- Reads are sorted on queryname by sambamba v0.5.6
- Duplicates are marked by Picardtools (ASSUME_SORT_ORDER=queryname)
- Reads are sorted on coordinate by sambamba
- a copy of the bamfile is made in which chromosome prefixes are removed for breakmer compatibility
- BreaKmer is performed using target regions of 5000 bp each

bamfiles need to be of the following make up:

stored in the folder bam
  {name}_coordsorted_nochr.bam
          and
      {name}_coordsorted.bam

Gridss, Novobreak and wham are performed by snakemake

[source activate translocation-snake]
[snakemake -s translocation.snakefile]


 
