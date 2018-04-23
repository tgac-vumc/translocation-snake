# Translocation detection pipeline by next generation sequencing
## A combined approach of BreaKmer, Gridss, Wham and novoBreak.

##### Work in progress :)


The translocation detection pipeline is designed to increase the sensitivity and specificity of translocation detection by using a combined approach of multiple publically available bioinformatic tools. This is especially important for the analysis of DNA derived from formalin fixed paraffin embedded (FFPE) material, because this material is in general of less quality compared to fresh material and contain more noise. The analysis is

**Requirements**

- Snakemake http://snakemake.readthedocs.io/en/stable/

    for easy installation you need (mini)conda.

    install Miniconda :
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
    follow the instructions of the installation process, give the location where you want miniconda to be installed and answer YES to add miniconda to your path.

    Snakemake and the other programs required are installed in an conda environment using conda for this you need the environment file of the translocation-snake:

    go to the directory where the analysis need to be performed and download the translocation-snake
    ```
        cd </path/to/analysis/dir>
        git clone https://github.com/tgac-vumc/translocation-snake.git
        cd translocation-snake

        conda env create --name translocation-snake --file environment.yaml
    ```

unfortunately not all programs are yet available through conda the following programs need to be installed manually, descriptions can be found at there own github.

- Gridss v1.4.2  https://github.com/PapenfussLab/gridss/releases
- WHAM https://github.com/zeeev/wham
- BreaKmer v0.0.4 https://github.com/ryanabo/BreaKmer/archive/v0.0.4-beta.zip
- Novobreak v1.1.3   https://github.com/czc/nb_distribution
    For novobreak replace infer_bp_v4.pl in the novobreak installation folder with the infer_bp_v4.pl file which you can find in the bin folder.
- R package StructuralVariantAnnotation and bedr aren't available in conda and need to be installed within R:

 ```
source activate translocation-snake
R

library(devtools)
options(unzip = 'internal')
install_github("PapenfussLab/StructuralVariantAnnotation"
install.packages("bedr")

q()
 ```
**Input file creation**

Pipeline is tested on paired-end ilumina data processed as follow:
- Adapters are removed with SeqPurge v0.1-104 (-min_len 20)
- DNA is aligned with bwa_mem v0.7.12-r1039 against reference hg19. (-M -R <readgroupinfo> )
- Realignment around indels is performed with ABRA v0.96
- Reads are sorted on queryname by sambamba v0.5.6
- Duplicates are marked by Picardtools (ASSUME_SORT_ORDER=queryname)
- Reads are sorted on coordinate by sambamba
- a copy of the bamfile is made in which chromosome prefixes are removed for breakmer compatibility
- BreaKmer is performed using target regions of 5000 bp each

bamfiles need to be stored in the folder bam located at the same level as the folder containing the snakefile (translocation-snake) and have names of the followin make-up:

  {name}_coordsorted_nochr.bam and {name}_coordsorted.bam

Gridss, Novobreak and wham are performed by snakemake

[source activate translocation-snake]

[snakemake ]




** installation tips/errors **

- install BreaKmer
	[wget https://github.com/ryanabo/BreaKmer/archive/v0.0.4-beta.zip]  
	[python setup.py install --user]
    - if error: Setup script exited with Missing required dependency NumPy (Numerical Python).   
    	[python2 -m pip install --user NumPy]
    - Breakmer requires Jellyfish we use https://github.com/gmarcais/Jellyfish/releases/download/v2.2.6/jellyfish-2.2.6.tar.gz
    - error door pysam: 0.9.0

    the pysam.sort() commands does not work with new versions of samtools, therefor two breakmer scripts needed to be changed:

    breakmer/processor/target.py
    	oorspronkelijk: 	pysam.sort(self.files['sv_bam'], self.files['sv_bam_sorted'].replace('.bam', ''))

   	New: pysam.sort("-o", self.files['sv_bam_sorted'], self.files['sv_bam'])

    breakmer/assembly/contig.py

    oorspronkelijk:  pysam.sort(bam_out_fn,bam_out_sorted_fn.replace('.bam',''))

    new:	pysam.sort("-o", bam_out_sorted_fn, bamOutFn)
- installation  novobreak
    - SSAKE need to be added to the PATH
