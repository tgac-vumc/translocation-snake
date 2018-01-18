configfile: "config.yaml"
from snakemake.utils import report
SAMPLES, = glob_wildcards("../bam/{sample}_coordsorted.bam")
'''
SAMPLES = ['DLBCL-27_S38_L003', 'DLBCL-90_S20_L001', 'DLBCL-20_S14_L002', 'DLBCL-122_S3', 'DLBCL-110_S34_L002', 'DLBCL-218_S76_L006', 'DLBCL-8_S8_L002', 'DLBCL-205_S79_L006', 'DLBCL-140_S35_L007', '10-028416_S66_L004', 'DLBCL-155_S41_L007', 'DLBCL-15_S22_L002', 'DLBCL-198_S52_L005', 'DLBCL-175_S32_L004', 'DLBCL-183_S64_L005', 'DLBCL-144_S27_L007', 'DLBCL-181_S44_L005', 'DLBCL-37_S32_L003', 'HE11-033749_S49_L003', 'DLBCL-139_S44_L007', 'DLBCL-171_S38_L004', 'DLBCL-215_S73_L006', 'T02-8834_S67_L004', 'DLBCL-182_S51_L005', 'DLBCL-159_S33_L007', 'DLBCL-209_S66_L006', 'DLBCL-56_S67', 'DLBCL-82_S9_L001', 'DLBCL-229_S85_L007', 'DLBCL-72A_S77', 'DLBCL-199_S62_L005', 'DLBCL-94_S12_L001', 'DLBCL-85_S10_L001', 'DLBCL-75_S76', 'DLBCL-97_S18_L001', 'DLBCL-123_S2', 'T10-10171_S33_L002', 'DLBCL-147_S28_L007', 'DLBCL-156_S39_L007', 'DLBCL-104_S27_L002', 'DLBCL-126_S20', 'DLBCL-190_S42_L005', 'DLBCL-233_S96_L007', 'DLBCL-84_S5_L001', 'DLBCL-112_S33_L002', '10-032737_S69_L004', 'DLBCL-131_S18', 'DLBCL-99_S13_L001', 'DLBCL-193_S55_L005', 'DLBCL-212_S70_L006', 'DLBCL-125_S21', 'DLBCL-107_S25_L002', 'DLBCL-92_S3_L001', 'DLBCL-18_S24_L002', 'DLBCL-116_S5', 'DLBCL-150_S38_L007', 'DLBCL-226_S92_L007', 'DLBCL-214_S77_L006', 'T09-23317_S64_L004', 'DLBCL-211_S69_L006', 'DLBCL-188_S59_L005', 'DLBCL-3_S3_L002', 'DLBCL-21_S5_L002', 'DLBCL-225_S95_L007', 'DLBCL-240_S98_L008', 'DLBCL-77_S75', 'DLBCL-22_S1_L002', '11-103151_S63_L004', 'DLBCL-141_S31_L007', 'DLBCL-138_S14', '11-104529_S48_L003', 'DLBCL-164_S41_L007', 'DLBCL-184_S63_L005', 'DLBCL-145_S30_L007', 'RH10-015545_S31_L002', 'DLBCL-113_S32_L002', 'DLBCL-87_S15_L001', 'DLBCL-133_S15', 'DLBCL-180_S33_L004', 'T10-14373_S52_L003', '11-101255_S44_L003', 'DLBCL-231_S86_L007', 'DLBCL-47_S36_L003', 'DLBCL-236_S97_L007', 'DLBCL-239_S102_L008', 'DLBCL-222_S90_L007', 'RH09-035008_S24_L002', 'DLBCL-120_S9', 'DLBCL-70_S72', 'DLBCL-42_S31_L003', 'DLBCL-9_S23_L002', 'DLBCL-13_S17_L002', 'DLBCL-143_S29_L007', 'DLBCL-223_S89_L007']
'''

rule all:
    input:
        #expand("../merged/{sample}-trl_merged.csv",sample=SAMPLES),
        expand("../reports/{sample}-report.html",sample=SAMPLES),
        #expand("../reports/{sample}-circlize.png",sample=SAMPLES),
	#expand("../breakmer/{sample}/{sample}.cfg", sample=SAMPLES)
        #expand("../merged/{sample}-trl_summary_brkpt_freq.csv",sample=SAMPLES),
        "../reports/circlize/index.html"


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
        ref=config['ALL']['REF_CHR'],
        assembly="../gridss/{sample}/{sample}-gridss.assembly.bam"
    threads:config["ALL"]["THREADS"]
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
        Annotationfile=config['ALL']["Annotationfile"],
        targets=config['ALL']['targets']
    script:
        '{input.script}'

rule run_wham:
    input:
        bam="../bam/{sample}_coordsorted_nochr.bam"
    output:
        vcf="../wham/{sample}/{sample}-wham.vcf"
    log:
        "../wham/log/{sample}.log"
    params:
        WHAM=config['wham']['WHAM'],
        ref=config['ALL']['REF_NOCHR'],
        mapqual=config['wham']['MAPQUAL'],
        basequal=config['wham']['BASEQUAL']
    threads:
        config["ALL"]["THREADS"]
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
        config["ALL"]["THREADS"]
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
        Annotationfile=config['ALL']["Annotationfile"],
        targets=config['ALL']['targets']
    script:
        '{input.script}'

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
        REF=config["ALL"]["REF_NOCHR"],
        HEADER=config["novobreak"]["HEADER"]
    threads:
        config["ALL"]["THREADS"]
    shell:
        "{params.NOVOBREAK} {params.EXE_DIR} {params.REF} {input.bam} "
        "{params.NORMAL} {threads} ../novobreak/{wildcards.sample} 2> {log} && "
        " cat {params.HEADER} ../novobreak/{wildcards.sample}/ssake/split/*.sp.vcf > {output.vcf} "

rule reorder_novobreak:
    input:
        vcf="../novobreak/{sample}/{sample}-novobreak.vcf",
        script="code/reorder-novobreak-snake.R"
    output:
        dups="../novobreak/{sample}/{sample}-novobreak_dups.csv",
        ordered="../novobreak/{sample}/{sample}-novobreak_ordered.csv"
    params:
        Annotationfile=config['ALL']["Annotationfile"],
        targets=config['ALL']['targets']
    script:
        '{input.script}'

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
	    ref=config["ALL"]["REF_NOCHR"],
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
        Annotationfile=config['ALL']["Annotationfile"]
    script:
        '{input.script}'

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
        summary="../merged/{sample}-trl_summary.csv",
        bedfile="../merged/{sample}-trl.bed"
    script:
        '{input.script}'

rule obtain_breakpoint_cov:
    input:
        bedfile="../merged/{sample}-trl.bed",
        bam="../bam/{sample}_coordsorted.bam"
    output:
        cov="../merged/{sample}-trl_cov.txt"
    shell:
        'samtools view -b -L {input.bedfile} -F 0x400 {input.bam} | bedtools coverage -d -b stdin -a {input.bedfile} > {output.cov}'

rule Calculate_brkpt_freq:
    input:
    	summary=temp("../merged/{sample}-trl_summary.csv"),
        cov="../merged/{sample}-trl_cov.txt",
        script='code/calculate_brkpt_freq.R'
    output:
        brkpt_freq="../merged/{sample}-trl_summary_brkpt_freq.csv"
    script:
        '{input.script}'

rule report:
    input:
        summary="../merged/{sample}-trl_summary_brkpt_freq.csv",
        circlize="../reports/circlize/{sample}-circlize.png",
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
        summary="../merged/{sample}-trl_summary_brkpt_freq.csv",
        script="code/circlize.R"
    output:
        circlize="../reports/circlize/{sample}-circlize.png"
    script:
        '{input.script}'

rule lightBox_circlize:
    input:
        circlize=expand("../reports/circlize/{sample}-circlize.png", sample=SAMPLES),
        script='code/createLightBox.sh',
    output:
        index="../reports/circlize/index.html",
    params:
        profiles="../reports/circlize/",
        lb2dir="lb2/",
    shell:
        "{input.script} {params.profiles} {params.lb2dir} > {output.index}"
