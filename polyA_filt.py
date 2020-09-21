"""Filter out transcripts without polyA tails (a small percentage)
AND/OR also remove those without/ with incorrect 5' overhang sequence: "ATGGG"   """

__author__ = "Sarah Hazell Pickering (s.h.pickering@medisin.uio.no)"
__date__ = "2020-09-21"

DIR = config["raw_dir"]


rule all:
    input:
        expand("{dir}/{prfx}{sample}{sfx}.polyA.fastq",
                dir = DIR,
                prfx = config["datatype_prefix"], #if any
                sample=config["samples"],
                sfx = config["datatype_suffix"]) #if any

##For ccs reads

rule filter_non_polyA:
    """filter non_polyA transcripts and those missing a 5' ATGGG,
       adding a suffix to the filename. 
    """
    input:
        DIR + "/{sample}.fastq"
    output:
        DIR + "/{sample}.polyA.fastq"
    shell:
        "egrep '^ATGGG.+A{{20}}$' -B 1 -A 2 --no-group-separator {input} "
            "> {output} "


###For hq transcripts
rule all_polyA_only:
    input:
	    expand("{dir}/polyA_only/{prfx}{sample}{sfx}.fastq",
	            dir = DIR,
	            prfx = config["datatype_prefix"], #if any
	            sample=config["samples"],
	            sfx = config["datatype_suffix"]) #if any


rule polyA_only:
    """filter non_polyA transcripts and deposit in a new dir"""
    input:
         DIR + "/{sample}.fastq"
    output:
        DIR + "/polyA_only/{sample}.fastq"
    shell:
        "egrep 'A{{20}}$' -B 1 -A 2 --no-group-separator {input} "
            "> {output} " 
