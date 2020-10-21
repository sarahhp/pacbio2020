"""Rules to collapse isoforms using tama (ported from the tofu cupcake script), 
And apply a custom filter to remove sense transcripts,
and single exon transcripts  (from sam or gff files).  
Also formats the final gff file to include the number 
of full length reads contributing to each transcript.  
"""

__author__ = "Sarah Hazell Pickering (s.h.pickering@medisin.uio.no)"
__date__ = "2020-10-14"

FASTQ_DIR = config["fastq_dir"]
INDIR = config["cpc_indir"]
OUTDIR = "tama"

#import pandas
#import get_fastq_sum as tools

rule tama_all:
    input:
        expand(OUTDIR + "/{sample}/{sample}.bed",
                sample = config["samples"]),
        expand(OUTDIR + "/{sample}/{sample}.cfilt.bed",
                sample =  config["samples"]),
        expand(OUTDIR + "/{sample}/{sample}.min_fl_50.bed",
                sample = config["samples"])

rule raw_abundance:
    input:
        expand(OUTDIR + "/{sample}/{sample}_read_support.txt",
                sample = config["samples"])

rule tama_collapse:
    input: 
        short_genome = config["genome_dir"] + "/chr12.Homo_sapiens.GRCh38.fa",
        sam = INDIR + "/hotair_selection/{sample}_hq.sam"
    params:
        out_prefix = lambda wildcards, output: output.bed[:-4],
        extra = config["params"]["tama"],
        tama_dir = config["tama_dir"]
    conda: "envs/tama.yml"
    output:
        report = OUTDIR + "/{sample}/{sample}_trans_report.txt",
        bed = OUTDIR + "/{sample}/{sample}.bed",
        transbed =  OUTDIR + "/{sample}/{sample}_trans_read.bed" 
    shell:
        "python {params.tama_dir}/tama_collapse.py "
            "-s {input.sam} -f {input.short_genome} "
            "{params.extra} -p {params.out_prefix} "


rule custom_tama_filter:
    input:
        bed = "{dir}/{sample}.bed",
    output:
        #ncluded = "{dir}/{sample}.included.txt",
        filt = "{dir}/{sample}.cfilt.bed"
    shell:
        "awk '$6 == \"-\"' {input.bed} | "   #select - strand
            "awk '$10 > 1' > "               #must have > 1 exon
            "{output.filt}; " 


rule generate_counts:
    input:
       report = "new_transcripts/reports/{sample}_polished.cluster_report.csv",
       transbed =  OUTDIR + "/{sample}/{sample}_trans_read.bed"
    params:
        script_dir = config["tama_dir"] + "/tama_go/read_support",
        inprefix = lambda wildcards, input: input.transbed[:-15]
    output:
        config = "tama/config_files/{sample}.list.txt",
        tab = OUTDIR + "/{sample}/{sample}_read_support.txt"
    shell:
        "echo 'cluspol	{input.report}	cluster' > {output.config} ; "
        "python {params.script_dir}/tama_read_support_levels.py "
        "-f {output.config} -m {input.transbed} -o {params.inprefix} "

rule filter_by_cov:
    input:
        bed = OUTDIR + "/{sample}/{sample}.cfilt.bed",
        counts = OUTDIR + "/{sample}/{sample}_read_support.txt",
    wildcard_constraints:
        threshold = "\d+"
    output:
        included = OUTDIR + "/config_files/{sample}.included_min_fl_{threshold}.txt",
        filtered = OUTDIR + "/{sample}/{sample}.min_fl_{threshold}.bed"
    shell:
       "awk '$4 > 50' {input.counts} | cut -f 2 > {output.included} ; " 
       "grep -wf {output.included} {input.bed} > {output.filtered} ; "

