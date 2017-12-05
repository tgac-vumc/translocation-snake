configfile: "config.yaml"
from snakemake.utils import report
SAMPLES, = glob_wildcards("../bam/{sample}_coordsorted.bam")
#SAMPLES = ["",""]

rule all:
    input:
        #expand("../merged/{sample}-trl_merged.csv",sample=SAMPLES),
        expand("../reports/{sample}-report.html",sample=SAMPLES),
        #expand("../reports/{sample}-circlize.png",sample=SAMPLES),
	    #expand("../breakmer/{sample}/{sample}.cfg", sample=SAMPLES)

rule run_gridss:
    input:
        bam="../bam/{sample}_coordsorted.bam"
    output:
        vcf="../gridss/{sample}/{sample}-gridss.vcf",
    log:
        "../gridss/log/{sample}.log"
    params:
        GRIDSS_JAR=config['gridss']['GRIDSS_JAR'],
        DIR="../gridss/{sample}/",
        ref=config['all']['REF_CHR'],
        assembly="../gridss/{sample}/{sample}-gridss.assembly.bam"
    threads:config["all"]["THREADS"]
    shell:
        "(java -ea -Xmx31g "
        "-Dsamjdk.create_index=true "
        "-Dsamjdk.use_async_io_read_samtools=true "
        "-Dsamjdk.use_async_io_write_samtools=true "
        "-Dsamjdk.use_async_io_write_tribble=true "
        "-cp {params.GRIDSS_JAR} gridss.CallVariants "
        "TMP_DIR= {params.DIR} "
        "WORKING_DIR= {params.DIR} "
        "REFERENCE_SEQUENCE={params.ref} "
        "INPUT={input.bam} "
        "OUTPUT={output.vcf} "
        "ASSEMBLY={params.assembly} "
        "THREADS={threads} "
        ") 2> {log} "

rule reorder_gridss:
    input:
        vcf="../gridss/{sample}/{sample}-gridss.vcf",
        script="code/reorder-gridss-snake.R"
    output:
        dups="../gridss/{sample}/{sample}-gridss_dups.csv",
        ordered="../gridss/{sample}/{sample}-gridss_ordered.csv"
    params:
        Annotationfile=config['all']["Annotationfile"],
        targets=config['all']['targets']
    script:
        'code/reorder-gridss-snake.R'

rule run_wham:
    input:
        bam="../bam/{sample}_coordsorted_nochr.bam"
    output:
        vcf="../wham/{sample}/{sample}-wham.vcf"
    log:
        "../wham/log/{sample}.log"
    params:
        WHAM=config['wham']['WHAM'],
        ref=config['all']['REF_NOCHR'],
        mapqual=config['wham']['MAPQUAL'],
        basequal=config['wham']['BASEQUAL']
    threads:
        config["all"]["THREADS"]
    shell:
        " {params.WHAM} -f {params.ref} "
        " -p {params.mapqual} -q {params.basequal} "
        " -x {threads} -t {input.bam} > {output.vcf} 2> {log} "

rule classify_wham:
    input:
        vcf="../wham/{sample}/{sample}-wham.vcf"
    output:
        classified="../wham/{sample}/{sample}-wham-class.vcf"
    log:
        "../wham/log/{sample}_class.log"
    params:
        TRAIN=config['wham']['TRAIN'],
        CLASSIFY=config['wham']['CLASSIFY'],
    threads:
        config["all"]["THREADS"]
    shell:
        "python2 {params.CLASSIFY} --proc {threads} {input.vcf} {params.TRAIN} "
        "> {output.classified} 2> {log} "

rule reorder_wham:
    input:
        classified="../wham/{sample}/{sample}-wham-class.vcf",
        script="code/reorder-wham-snake.R"
    output:
        dups="../wham/{sample}/{sample}-wham_dups.csv",
        ordered="../wham/{sample}/{sample}-wham_ordered.csv"
    params:
        Annotationfile=config['all']["Annotationfile"],
        targets=config['all']['targets']
    script:
        'code/reorder-wham-snake.R'

rule run_novobreak:
    input:
        bam="../bam/{sample}_coordsorted_nochr.bam"
    output:
        vcf="../novobreak/{sample}/{sample}-novobreak.vcf"
    log:
        "../novobreak/log/{sample}.log"
    params:
        NOVOBREAK=config['novobreak']['NOVOBREAK'],
        NORMAL=config['novobreak']['NORMAL'],
        EXE_DIR=config['novobreak']['EXE_DIR'],
        REF=config["all"]["REF_NOCHR"],
        HEADER=config["novobreak"]["HEADER"]
    threads:
        config["all"]["THREADS"]
    shell:
        "{params.NOVOBREAK} {params.EXE_DIR} {params.REF} {input.bam} "
        "{params.NORMAL} {threads} ../novobreak/{wildcards.sample} 2> {log} ; "
        " cat {params.HEADER} ../novobreak/{wildcards.sample}/ssake/split/*.sp.vcf > {output.vcf} "

rule reorder_novobreak:
    input:
        vcf="../novobreak/{sample}/{sample}-novobreak.vcf",
        script="code/reorder-novobreak-snake.R"
    output:
        dups="../novobreak/{sample}/{sample}-novobreak_dups.csv",
        ordered="../novobreak/{sample}/{sample}-novobreak_ordered.csv"
    params:
        Annotationfile=config['all']["Annotationfile"],
        targets=config['all']['targets']
    script:
        'code/reorder-novobreak-snake.R'

rule create_config_breakmer:
   input:
        bam="../bam/{sample}_coordsorted_nochr.bam"
   output:
        config="../breakmer/{sample}/{sample}.cfg"
        #"../breakmer/Breakmer_output/{sample}_BCNHL_Seq_V2_allTRL_svs.out"
   params:
        targets_bed=config["breakmer"]["targets_bed_file"],
        #analysis_name="{sample}"+config["breakmer"]["analysis_name"],
        analysis_name="{sample}_BCNHL_Seq_V2_allTRL",
	    refdir=config["breakmer"]["reference_data_dir"],
        cutadapt_config=config["breakmer"]["cutadapt_config_file"],
        cutadapt=config["breakmer"]["cutadapt"],
	    jellyfish=config["breakmer"]["jellyfish"],
	    ref=config["all"]["REF_NOCHR"],
        annotation=config["breakmer"]["gene_annotation_file"],
        kmer=config["breakmer"]["kmer_size"]
   shell:
      'echo "analysis_name={params.analysis_name} \n\
targets_bed_file={params.targets_bed} \n\
sample_bam_file={input.bam} \n\
analysis_dir=../breakmer/{wildcards.sample} \n\
reference_data_dir={params.refdir} \n\
cutadapt_config_file={params.cutadapt_config} \n\
cutadapt={params.cutadapt} \n\
jellyfish={params.jellyfish} \n\
blat=blat \n\
gfclient=gfClient \n\
gfserver=gfServer \n\
fatotwobit=faToTwoBit \n\
reference_fasta={params.ref} \n\
gene_annotation_file={params.annotation} \n\
kmer_size={params.kmer}\
" > {output.config}'

rule run_breakmer:
    input:
        bam="../bam/{sample}_coordsorted_nochr.bam",
        config="../breakmer/{sample}/{sample}.cfg",
    output:
         out="../breakmer/Breakmer_output/{sample}_BCNHL_Seq_V2_allTRL_svs.out"
    params:
        breakmer=config["breakmer"]["BREAKMER"],
         #analysis_name="{sample}"+config["breakmer"]["analysis_name"]
        analysis_name="{sample}_BCNHL_Seq_V2_allTRL"
    shell:
         "python2 {params.breakmer} run -c {input.config} ;"
         "cp ../breakmer/{wildcards.sample}/output/{params.analysis_name}_svs.out {output.out}"

rule reorder_breakmer:
    input:
        "../breakmer/Breakmer_output/{sample}_BCNHL_Seq_V2_allTRL_svs.out",
        script="code/reorder-breakmer004-snake.R"
    output:
        dups="../breakmer/{sample}/{sample}-breakmer_dups.csv",
        ordered="../breakmer/{sample}/{sample}-breakmer_ordered.csv",
        complexe="../breakmer/{sample}/{sample}-breakmer_complex.csv"
    params:
        Annotationfile=config['all']["Annotationfile"]
    script:
        'code/reorder-breakmer004-snake.R'

rule merge_svtools:
    input:
        novobreak="../novobreak/{sample}/{sample}-novobreak_ordered.csv",
        breakmer="../breakmer/{sample}/{sample}-breakmer_ordered.csv",
        wham="../wham/{sample}/{sample}-wham_ordered.csv",
        gridss="../gridss/{sample}/{sample}-gridss_ordered.csv",
        script="code/merge-svtools-snake.R"
    output:
        trl="../merged/{sample}-trl_merged.csv",
        IG="../merged/{sample}-IG_merged.csv",
        otherSVs="../merged/{sample}-other_merged.csv",
        summary="../merged/{sample}-trl_summary.csv"
        bedfile="./merged/{sample}-trl.bed"
    script:
        'code/merge-svtools-snake.R'

rule Calculate_VAF:
    input:
        summary="../merged/{sample}-trl_summary.csv",
        bedfile="./merged/{sample}-trl.bed",
        bam="../bam/{sample}_coordsorted.bam"
    output:
        summary="../merged/{sample}-trl_summary-vaf.csv"
    run:
        'samtools view -b -F 0x400 {input.bam} | bedtools coverage -d -abam stdin -b {input.bed}'


rule report:
    input:
        summary="../merged/{sample}-trl_summary.csv",
        circlize="../reports/{sample}-circlize.png",
        translocations="../merged/{sample}-trl_merged.csv",
        IG_rearrangments="../merged/{sample}-IG_merged.csv",
        otherSVs_low_evidence="../merged/{sample}-other_merged.csv"

    output:
        "../reports/{sample}-report.html"
    run:                                     #run instead of shell, so plain code will run directly.
        report("""

        ===================================================
        Translocations in {wildcards.sample}
        ===================================================

        Reads were mapped to hg19 with BWA mem

        translocations are called with novoBreak, BreaKmer, wham and gridss.

        Used filters were:

        - BreaKmer splitread: {config[filters][splitread]} and discordant reads: {config[filters][discread]}
        - novoBreak a minimal coverage score above: {config[filters][novoscore]}
        - Gridss a minimal quality score of: {config[filters][gridssscore]}
        - wham: the sum of reads is above two times {config[filters][splitread]}

        Events are devided in high and low evidence, an event is called high evidence in case the amount of discordant
        reads and split reads are both above: {config[filters][highSVevidence]}

        This resulted in  variants.

        """, output[0], **input)


rule circlize:
    input:
        summary="../merged/{sample}-trl_summary.csv",
        script="code/circlize.R"
    output:
        circlize="../reports/{sample}-circlize.png"
    script:
        'code/circlize.R'
