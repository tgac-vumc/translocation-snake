configfile: "config.yaml"
from snakemake.utils import report
SAMPLES, = glob_wildcards("../bam/{sample}_coordsorted.bam")
targetregionnumber=['{:02d}'.format(item) for item in list(range(config['ALL']['THREADS']))]


rule all:
    input:
        expand("../reports/{sample}-report.html",sample=SAMPLES),
        #expand("../merged/{sample}-trl_summary_brkpt_freq.csv",sample=SAMPLES),
        "../reports/circlize/index.html",
        "../reports/summary.html",
        #expand("../CovMetrics/{sample}_Combined_metrics.txt" , sample=SAMPLES),
        "../CovMetrics/All_samples_Combined_metrics.txt",

rule run_gridss:
    input:
        bam="../bam/{sample}_coordsorted.bam"
    output:
        vcf="../gridss/{sample}/{sample}-gridss.vcf",
        Alignment_met="../gridss/{sample}/{sample}_coordsorted.bam.gridss.working/{sample}_coordsorted.bam.alignment_summary_metrics",
        insert_size_metrics="../gridss/{sample}/{sample}_coordsorted.bam.gridss.working/{sample}_coordsorted.bam.insert_size_metrics",
        sv_metrics="../gridss/{sample}/{sample}_coordsorted.bam.gridss.working/{sample}_coordsorted.bam.sv_metrics",
    log:
        "../gridss/log/{sample}.log"
    params:
        GRIDSS_JAR=config['gridss']['GRIDSS_JAR'],
        DIR="../gridss/{sample}/",
        ref=config['ALL']['REF'],
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
        bam="../bam/{sample}_coordsorted.bam"
    output:
        vcf="../wham/{sample}/{sample}-wham.vcf"
    log:
        "../wham/log/{sample}.log"
    params:
        WHAM=config['wham']['WHAM'],
        ref=config['ALL']['REF'],
        mapqual=config['wham']['MAPQUAL'],
        basequal=config['wham']['BASEQUAL']
    threads:
        config["ALL"]["THREADS"]
    shell:
        " {params.WHAM} -f {params.ref} "
        " -p {params.mapqual} -q {params.basequal} "
        " -x {threads} -t {input.bam} > {output.vcf} 2> {log} "

rule reorder_wham:
    input:
        vcf="../wham/{sample}/{sample}-wham.vcf",
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
        bam="../bam/{sample}_coordsorted.bam"
    output:
        vcf="../novobreak/{sample}/{sample}-novobreak.vcf"
    log:
        "../novobreak/log/{sample}.log"
    params:
        NOVOBREAK=config['novobreak']['NOVOBREAK'],
        NORMAL=config['novobreak']['NORMAL'],
        EXE_DIR=config['novobreak']['EXE_DIR'],
        REF=config["ALL"]["REF"],
        HEADER=config["novobreak"]["HEADER"]
    threads:
        config["ALL"]["THREADS"]
    shell:
        "{params.NOVOBREAK} {params.EXE_DIR} {params.REF} {input.bam} "
        "{params.NORMAL} {threads} ../novobreak/{wildcards.sample} &> {log} && "
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

rule create_targetfiles:
	input:
		config['breakmer']['targets_bed_file']
	output:
		temp(expand("../breakmer/{{sample}}/targetregions-{number}",number=targetregionnumber))
	params:
		threads=config["ALL"]["THREADS"],
		output="../breakmer/{sample}/targetregions-"
	shell:
		'split -d -n l/{params.threads} {input} {params.output}'


rule create_config_breakmer:
   input:
        bam="../bam/{sample}_coordsorted_nochr.bam",
        targetfile="../breakmer/{sample}/targetregions-{number}"
   output:
        config="../breakmer/{sample}/{number}/{sample}-{number}.cfg"
        #"../breakmer/Breakmer_output/{sample}_BCNHL_Seq_V2_allTRL_svs.out"
   params:
        targets_bed="../breakmer/{sample}/targetregions-{number}",
        #analysis_name="{sample}"+config["breakmer"]["analysis_name"],
        analysis_name="{sample}_BCNHL_Seq_V2_allTRL-{number}",
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
analysis_dir=../breakmer/{wildcards.sample}/{wildcards.number} \n\
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
        config="../breakmer/{sample}/{number}/{sample}-{number}.cfg",
        targetfile="../breakmer/{sample}/targetregions-{number}",
    output:
         out="../breakmer/Breakmer_output/{sample}_BCNHL_Seq_V2_allTRL-{number}_svs.out"
    params:
        breakmer=config["breakmer"]["BREAKMER"],
         #analysis_name="{sample}"+config["breakmer"]["analysis_name"]
        analysis_name="{sample}_BCNHL_Seq_V2_allTRL"
    shell:
         "python2 {params.breakmer} run -c {input.config} ;"
         "cp ../breakmer/{wildcards.sample}/{wildcards.number}/output/{params.analysis_name}-{wildcards.number}_svs.out {output.out}"

rule concat_breakmer:
	input:
		expand("../breakmer/Breakmer_output/{{sample}}_BCNHL_Seq_V2_allTRL-{number}_svs.out", number=targetregionnumber)
	output:
		"../breakmer/Breakmer_output/{sample}_BCNHL_Seq_V2_allTRL_svs.out",
	params:
		header="code/breakmer_header.txt"
	shell:
		'awk FNR-1 {params.header} {input} > {output} '


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
        breakmer="../breakmer/{sample}/{sample}-breakmer_ordered.csv" if config["options"]["breakmer"] else "headers/empty_breakmer.csv",
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

rule Metrics:
	input:
		HSmetrics="../CovMetrics/{sample}_HSmetrics.txt",
		Alignment_met="../gridss/{sample}/{sample}_coordsorted.bam.gridss.working/{sample}_coordsorted.bam.alignment_summary_metrics",
		insert_size_metrics="../gridss/{sample}/{sample}_coordsorted.bam.gridss.working/{sample}_coordsorted.bam.insert_size_metrics",
		sv_metrics="../gridss/{sample}/{sample}_coordsorted.bam.gridss.working/{sample}_coordsorted.bam.sv_metrics",
		script="code/Combine_metrics.R",
	output:
		combined="../CovMetrics/{sample}_Combined_metrics.txt"
	script:
		'{input.script}'

rule All_metrics:
	input:
		combined=sorted(expand("../CovMetrics/{sample}_Combined_metrics.txt", sample=SAMPLES))
	output:
		"../CovMetrics/All_samples_Combined_metrics.txt"
	shell:
		"head -1 {input.combined[0]} > {output} && tail -n +2 -q {input.combined} >> {output}"


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

rule summary:
    input:
        script='code/summary.sh',
        break_freq=expand("../merged/{sample}-trl_summary_brkpt_freq.csv", sample=SAMPLES),
    output:
        "../reports/summary.html"
    params:
        bamfolder="../bam/'*'_coordsorted.bam",
    shell:
        "{input.script} Translocations {params.bamfolder} > {output}"

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


# rule Covmetrics:
#     input:
#         script1='code/splitHSmetrics.sh',
#         script2='code/c'
#         hsmetrics="../Covmetrics/{sample}_HSmetrics.txt",
#         perTargetCov="../Covmetrics/{sample}_PerTargetCov.txt",
#     output:
#         HsMetrics="../Covmetrics/{sample}_HsMetrics.txt",
#         histogram="../Covmetrics/{sample}_CovHisto.txt"
#     shell:
#         ""

#rule classify_wham:
#    input:
#        vcf="../wham/{sample}/{sample}-wham.vcf"
#    output:
#        classified="../wham/{sample}/{sample}-wham-class.vcf"
#    log:
#        "../wham/log/{sample}_class.log"
#    params:
#        TRAIN=config['wham']['TRAIN'],
#        CLASSIFY=config['wham']['CLASSIFY'],
#    threads:
#        config["ALL"]["THREADS"]
#    shell:
#        "python2 {params.CLASSIFY} --proc {threads} {input.vcf} {params.TRAIN} "
#        "> {output.classified} 2> {log} "
