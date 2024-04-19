from subprocess import call
import sys
import os
import re
from itertools import chain
from pathlib import Path
from itertools import combinations
shell.prefix("source ~/.bash_profile; ")

if not os.path.exists("00log"):
    os.makedirs("00log")
if not os.path.exists("01seq"):
    os.makedirs("01seq")
if not os.path.exists("02fqc"):
    os.makedirs("02fqc")
if not os.path.exists("03aln"):
    os.makedirs("03aln")
if not os.path.exists("04quant"):
    os.makedirs("04quant")
if not os.path.exists("05norm"):
    os.makedirs("05norm")
if not os.path.exists("06stats"):
    os.makedirs("06stats")
if not os.path.exists("07vis"):
    os.makedirs("07vis")

ALL_GROUPS = config["groups"]
GROUP_NAMES = list(ALL_GROUPS.keys())

ALL_SAMPLES = list(chain(*ALL_GROUPS.values()))
comb = list(combinations(GROUP_NAMES,2))
ALL_DAY_ZERO = ALL_GROUPS[config["baseline"]]
#print(ALL_DAY_ZERO)
ALL_DAY_ZERO =  ",".join(ALL_DAY_ZERO)

ALL_COMB = {}
ALL_COMB_FILTERED = {}

for c in comb:
    ALL_COMB[c[0] + "_vs_" + c[1]] = [c[0],c[1]]
    if c[0] == config["baseline"]:
        ALL_COMB_FILTERED[c[0] + "_vs_" + c[1]] = [c[0],c[1]]

ALL_COMB_FILTERED_CONTROL_VS_TREATMENT = {}

for k, v in ALL_COMB_FILTERED.items():
    if re.search(config["control_postfix"] + '$' , v[1]):
        for k2, v2 in ALL_COMB_FILTERED.items():
            if not(v[1] == v2[1]): 
                ALL_COMB_FILTERED_CONTROL_VS_TREATMENT[v[0] + "_vs_" + v[1] + "_COMPARED_" + v2[0] + "_vs_" + v2[1]] = [v[0],v[1],v2[0],v2[1]]


COMB_NAMES = list(ALL_COMB.keys())
if not("stranded" in config):
   config["stranded"] = False 

if not("featureCounts") in config:
    config["featureCounts"] = False  

if not("all_genes") in config:
    config["all_genes"] = True 

if not("genes_to_plot") in config:
    genes = "CD22,ATR,CHEK2,CHEK1,ATM"
else:
    genes  = config["genes_to_plot"]
    

ALL_FASTQC = expand("02fqc/{sample}_fastqc.zip".split(), sample = ALL_SAMPLES)
ALL_TRIMMED = expand("01seq/{sample}_trimmed.fastq".split(), sample = ALL_SAMPLES)
ALL_ALIGNED = expand("03aln/{sample}.bam".split(), sample = ALL_SAMPLES)
ALL_SAMPLE_NAMES = ",".join(ALL_SAMPLES)


QUANT = ["04quant/mageck.count_filtered.txt",
        "04quant/mageck_median_filtered.count_normalized.txt",
        "04quant/mageck_total_filtered.count_normalized.txt"]

ALL_VIS = ["07vis/not_normalized_combined.pdf",
            "07vis/not_normalized_filtered_combined.pdf", 
            "07vis/median_filtered_combined.pdf",
            "07vis/total_filtered_combined.pdf"]

if not("raw_counts" in config):
    QUANT.append("04quant/mageck.count.txt")
else:
    if not os.path.exists("04quant/mageck.count.txt"):
        os.system("cp " + config["raw_counts"] + " 04quant/mageck.count.txt")


control_sgRNA = workflow.basedir + "/published_data/dummy.txt"
essential_genes = workflow.basedir + "/published_data/dummy.txt"
if ("essential_genes" in config):
    essential_genes = config["essential_genes"]
#print(control_sgRNA)
#exit(0)
if ("mageck_control_sgRNA" in config["sgRNA"]):
    control_sgRNA = config["sgRNA"]["mageck_control_sgRNA"]
    QUANT.append("04quant/mageck_control_filtered.count_normalized.txt")
    ALL_VIS.append("07vis/control_filtered_combined.pdf")
    ALL_STARS_GENELIST = expand("06stats/STARS_control_NORM_{comb}_genelist.csv 06stats/STARS_total_NORM_{comb}_genelist.csv".split() ,comb = ALL_COMB_FILTERED_CONTROL_VS_TREATMENT)
    ALL_SCREEN_BEAM_GENELIST = expand("06stats/ScreenBEAM_default_NORM_{comb}_genelist.csv 06stats/ScreenBEAM_control_NORM_{comb}_genelist.csv 06stats/ScreenBEAM_total_NORM_{comb}_genelist.csv".split() ,comb = ALL_COMB_FILTERED_CONTROL_VS_TREATMENT)
    ALL_PBNPA_GENELIST = expand("06stats/pbnpa_total_NORM_{comb}_genelist.csv 06stats/pbnpa_not_normalized_NORM_{comb}_genelist.csv 06stats/pbnpa_control_NORM_{comb}_genelist.csv".split() ,comb = ALL_COMB_FILTERED_CONTROL_VS_TREATMENT)
    ALL_MAGECK_GENELIST = expand("06stats/mageck_total_NORM_{comb}_genelist.csv 06stats/mageck_default_NORM_{comb}_genelist.csv 06stats/mageck_control_NORM_{comb}_genelist.csv".split() ,comb = ALL_COMB_FILTERED_CONTROL_VS_TREATMENT)
    ALL_MAGECK = expand("06stats/mageck_default_{comb}.gene_summary.txt 06stats/mageck_control_{comb}.gene_summary.txt 06stats/mageck_total_{comb}.gene_summary.txt".split() ,comb = ALL_COMB)
else: 
    ALL_STARS_GENELIST = expand("06stats/STARS_total_NORM_{comb}_genelist.csv".split() ,comb = ALL_COMB_FILTERED_CONTROL_VS_TREATMENT)
    ALL_SCREEN_BEAM_GENELIST = expand("06stats/ScreenBEAM_default_NORM_{comb}_genelist.csv 06stats/ScreenBEAM_total_NORM_{comb}_genelist.csv".split() ,comb = ALL_COMB_FILTERED_CONTROL_VS_TREATMENT)
    ALL_PBNPA_GENELIST = expand("06stats/pbnpa_total_NORM_{comb}_genelist.csv 06stats/pbnpa_not_normalized_NORM_{comb}_genelist.csv".split() ,comb = ALL_COMB_FILTERED_CONTROL_VS_TREATMENT)
    ALL_MAGECK_GENELIST = expand("06stats/mageck_total_NORM_{comb}_genelist.csv 06stats/mageck_default_NORM_{comb}_genelist.csv".split() ,comb = ALL_COMB_FILTERED_CONTROL_VS_TREATMENT)
    ALL_MAGECK = expand("06stats/mageck_default_{comb}.gene_summary.txt 06stats/mageck_total_{comb}.gene_summary.txt".split() ,comb = ALL_COMB)

#print(QUANT)

ALL_BAM = expand("03aln/{sample}.bam".split() ,sample = ALL_SAMPLES)

ALL_PERT = expand("05norm/control_NORM_{comb}_pertubation.txt 05norm/total_NORM_{comb}_pertubation.txt".split() ,comb = ALL_COMB)

ALL_TOOLS = config["tools"]




ruleorder:  run_make_genelists_pbnpa > run_make_genelists_mageck > run_make_genelists_stars > run_make_genelists_screenbeam > run_mageck_total > run_mageck_control > run_mageck_default > run_PBNPA_default > run_PBNPA_total > run_PBNPA_control > run_stars_total > run_stars_control > run_screen_beam > run_stars_null_total > run_stars_null_control > get_pertubations_total > get_pertubations_control


if ("bowtie2_index" in config["sgRNA"]) :
    
    bowtie2_index = config["sgRNA"]["bowtie2_index"]

    rule all:
        input: ALL_FASTQC  + QUANT  + ALL_VIS +  ALL_PBNPA_GENELIST + ALL_MAGECK_GENELIST + ALL_SCREEN_BEAM_GENELIST + ALL_STARS_GENELIST 
else: 
    bowtie2_index = ""
    rule all:
        input: QUANT + ALL_MAGECK + ALL_VIS +  ALL_PBNPA_GENELIST + ALL_MAGECK_GENELIST + ALL_SCREEN_BEAM_GENELIST + ALL_STARS_GENELIST 
        


rule fastqc:
    input:  "01seq/{sample}.fastq.gz"
    output: "02fqc/{sample}_fastqc.zip"
    log:    "00log/{sample}.fastqc"
    threads: 6
    params: 
        mem = "2G"
    message: "run fastqc for {input[0]}"
    shell:
        """
        conda activate /home/hayerk/miniconda3/envs/snakemake_new2/envs/mageck-vispr
        fastqc -t 5 -o 02fqc {input[0]}  2> {log}
        """

rule cutadapt:
    input:
        "01seq/{sample}.fastq.gz"
    output:
        trimmed = "01seq/{sample}_trimmed.fastq",
        too_short = "01seq/{sample}_too_short.fastq",
        untrimmed = "01seq/{sample}_untrimmed.fastq"
    message: "run cutadapt for {input[0]}"
    threads: 2
    params: 
        mem = "2G"
    shell:
        """
        conda activate /home/hayerk/miniconda3/envs/snakemake_new2/envs/mageck-vispr
        cutadapt -O 5 -e 0.1 -m 20 -g {config[sgRNA][adapter]} -l {config[sgRNA][length]} \
        --too-short-output={output[too_short]} --untrimmed-output={output[untrimmed]} {input} > {output[trimmed]}
        """

rule align:
    input:  "01seq/{sample}_trimmed.fastq"
    output: "03aln/{sample}.bam"
    log:    "00log/{sample}.align"
    threads: 10
    params: 
        mem = "1G", 
        bowtie2 = " -x " + bowtie2_index
    message: "aligning {input}: {threads} threads / {params.mem}"
    shell:
        """
        conda activate bioinf
        which samtools
        {config[tools][bowtie2]} {params.bowtie2} -p {threads} -U {input[0]} -N 1 -L 15 --norc | samtools view -bS -q 5 -  2> {log}  >  {output}
        """

rule bam_to_counts:
    input: 
        lib = config["sgRNA"]["index_mageck"],
        bams = ALL_BAM,
        #control_sgRNA = config["sgRNA"]["mageck_control_sgRNA"]
    output: 
        original_count = "04quant/mageck.count.txt", 
        first_norm = "04quant/mageck.count_normalized.txt"
    log:    "00log/mageck_count"
    threads: 2
    params: 
        mem = "4G",
        filter = 1,
        day_zero = ALL_DAY_ZERO
    message: "bam_to_counts {input}: {threads} threads / {params.mem}"
    shell:
        """
        conda activate /home/hayerk/miniconda3/envs/snakemake_new2/envs/mageck-vispr
        mageck count --norm-method total -l {input.lib} -n 04quant/mageck --sample-label "{ALL_SAMPLE_NAMES}" --fastq {input.bams} --day0-label {params.day_zero}
        """

rule normalize_counts:
    input: 
        lib = config["sgRNA"]["index_mageck"],
        original_count = "04quant/mageck.count.txt"
    output:  
        original_filtered = "04quant/mageck.count_filtered.txt",
        median = "04quant/mageck_median_filtered.count_normalized.txt",
        total = "04quant/mageck_total_filtered.count_normalized.txt"
    log:    "00log/mageck_count"
    threads: 2
    params: 
        mem = "4G",
        filter = 1,
        day_zero = ALL_DAY_ZERO
    message: "normalizing {input}: {threads} threads / {params.mem}"
    shell:
        """
        conda activate /home/hayerk/miniconda3/envs/snakemake_new2/envs/mageck-vispr
        Rscript {workflow.basedir}/scripts/filter_rows.R {input.original_count} {output.original_filtered} {params.filter} {ALL_SAMPLE_NAMES}
        mageck count --norm-method total -l {input.lib} -n 04quant/mageck_total_filtered -k {output.original_filtered} --day0-label {params.day_zero}
        mageck count --norm-method median -l {input.lib} -n 04quant/mageck_median_filtered -k {output.original_filtered} --day0-label {params.day_zero}
        """


def get_sam_names(wildcards):

    k = expand("{sample}", sample = list(config["groups"][ALL_COMB[wildcards[0]][0]]))
    return ",".join(k)#expand("{sample}", sample = list(config["groups"][ALL_COMB[wildcards[0]][0]]))

def get_sam_names2(wildcards):
    k = expand("{sample}", sample = list(config["groups"][ALL_COMB[wildcards[0]][1]]))
    return ",".join(k)

def get_name1(wildcards):
    #print(wildcards[0])
    return ALL_COMB[wildcards[0]][0]

def get_name2(wildcards):
    return ALL_COMB[wildcards[0]][1]


rule normalize_all_counts_by_control:
    input: 
        lib = config["sgRNA"]["index_mageck"],
        counts = "04quant/mageck.count_filtered.txt",
        control_sgRNA = control_sgRNA
    output: 
        control = "04quant/mageck_control_filtered.count_normalized.txt"
    log: "00log/normalize_all_counts_by_control"
    threads: 2
    params:
        mem = "4G",
        filter = 1,
        day_zero = ALL_DAY_ZERO
    message: "normalize_all_counts_by_control {input}: {output} : {threads} threads"
    shell:
        """
        conda activate /home/hayerk/miniconda3/envs/snakemake_new2/envs/mageck-vispr
        mageck count --norm-method control --control-sgrna {input.control_sgRNA} -l {input.lib} -n 04quant/mageck_control_filtered -k {input.counts} --day0-label {params.day_zero}
        """  

rule normalize_counts_by_control:
    input: 
        counts = "04quant/mageck.count_filtered.txt",
        control_sgRNA = control_sgRNA
    output: 
        normed = "05norm/control_{comb}.normalized.txt"    
    log: "00log/{comb}_normalize_counts_by_control"
    threads: 2
    params:
        mem  = "5G",
        name1 = get_name1,
        name2 = get_name2,
        sam_names1 = get_sam_names,
        sam_names2 = get_sam_names2,
        control_sgRNA = control_sgRNA
    message: "normalize_counts_by_control {input}: {output} : {threads} threads"
    shell:
        """
        conda activate /home/hayerk/miniconda3/envs/snakemake_new2/envs/mageck-vispr
        mageck test --normcounts-to-file --norm-method control --control-sgrna {params.control_sgRNA} -k {input.counts} -t {params.sam_names2} -c {params.sam_names1} -n 05norm/control_{params.name1}_vs_{params.name2}  2> {log}
        """  



rule run_plot_counts:
    input:  
        not_normalized = '04quant/mageck.count.txt',
        not_normalized_filtered = '04quant/mageck.count_filtered.txt',
        total_filtered = '04quant/mageck_total_filtered.count_normalized.txt',
        median_filtered = '04quant/mageck_median_filtered.count_normalized.txt'
    output: "07vis/not_normalized_combined.pdf", "07vis/total_filtered_combined.pdf",
        "07vis/not_normalized_filtered_combined.pdf", "07vis/median_filtered_combined.pdf"
    log:    "00log/plot_counts"
    threads: 10
    params:
        mem  = "5G",
        essential_genes = essential_genes,
        prefix = "07vis/not_normalized",
        prefix_filtered = "07vis/not_normalized_filtered",
        prefix_median = "07vis/median_filtered",
        prefix_total = "07vis/total_filtered"
    message: "plot_counts {input}: {output} : {threads} threads" #"/ {params.mem}"
    shell: 
        """
        Rscript {workflow.basedir}/scripts/plot_counts.R {input.not_normalized} not_normalized {params.essential_genes} {params.prefix}
        Rscript {workflow.basedir}/scripts/plot_counts.R {input.not_normalized_filtered} not_normalized_filtered {params.essential_genes} {params.prefix_filtered}
        Rscript {workflow.basedir}/scripts/plot_counts.R {input.total_filtered} total {params.essential_genes} {params.prefix_total}
        Rscript {workflow.basedir}/scripts/plot_counts.R {input.median_filtered} median {params.essential_genes} {params.prefix_median} 
        """

rule run_plot_control:
    input:  
        control_normalized_filtered = "04quant/mageck_control_filtered.count_normalized.txt"
    output: "07vis/control_filtered_combined.pdf"
    log:    "00log/plot_counts"
    threads: 10
    params:
        mem  = "5G",
        essential_genes = essential_genes,
        prefix_control = "07vis/control_filtered"
    message: "plot_counts_control {input}: {output} : {threads} threads" #"/ {params.mem}"
    shell: 
        """
        Rscript {workflow.basedir}/scripts/plot_counts.R {input.control_normalized_filtered} control {params.essential_genes} {params.prefix_control}
        """


rule get_pertubations_control:
    input:  "05norm/control_{comb}.normalized.txt"
    output: "05norm/control_NORM_{comb}_pertubation.txt"
    log:    "00log/{comb}_pertubation"
    threads: 10
    params:
        mem  = "5G",
        name1 = get_name1,
        name2 = get_name2,
        sam_names1 = get_sam_names,
        sam_names2 = get_sam_names2
    message: "get_pertubations_control {input}: {output} : {threads} threads" #"/ {params.mem}"
    shell: 
        """
        grep -v NTCONTROL {input} > {input}_{params.name1}_{params.name2}_STARS_tmp
        ruby {workflow.basedir}/scripts/make_pertubations.rb {params.name1} {params.name2} {input}_{params.name1}_{params.name2}_STARS_tmp {params.sam_names1} {params.sam_names2} > {output}
        rm {input}_{params.name1}_{params.name2}_STARS_tmp
        """

rule get_pertubations_total:
    input:  "04quant/mageck_total_filtered.count_normalized.txt"
    output: "05norm/total_NORM_{comb}_pertubation.txt"
    log:    "00log/{comb}_pertubation_total"
    threads: 10
    params:
        mem  = "5G",
        name1 = get_name1,
        name2 = get_name2,
        sam_names1 = get_sam_names,
        sam_names2 = get_sam_names2
    message: "get_pertubations_total {input}: {output} : {threads} threads" #"/ {params.mem}"
    shell: 
        """
        grep -v NTCONTROL {input} > {input}_{params.name1}_{params.name2}_STARS_tmp
        ruby {workflow.basedir}/scripts/make_pertubations.rb {params.name1} {params.name2} {input}_{params.name1}_{params.name2}_STARS_tmp {params.sam_names1} {params.sam_names2} > {output}
        rm {input}_{params.name1}_{params.name2}_STARS_tmp
        """

rule run_stars_null_control:
    input:  
        pert = "05norm/control_NORM_{comb}_pertubation.txt",
        lib = config["sgRNA"]["star_chip"]
    output: 
        "06stats/control_NORM_{comb}_Null_STARSOutput8_10_N.txt",
        "06stats/control_NORM_{comb}_Null_STARSOutput8_10_P.txt"
    log:    "00log/control_{comb}_stars_null"
    threads: 10
    params:
        mem  = "5G",
        name1 = get_name1,
        name2 = get_name2
    message: "run_stars_null_control {input}: {output} : {threads} threads" #"/ {params.mem}"
    shell: 
        """
        conda activate stars_env
        python {workflow.basedir}/STARS_v1.3/stars_null_v1.3.py --input-file {input.pert} --chip-file {input.lib} --thr 20 --num-ite 1000 --dir P --out-file {output[0]} 2> {log}
        python {workflow.basedir}/STARS_v1.3/stars_null_v1.3.py --input-file {input.pert} --chip-file {input.lib} --thr 20 --num-ite 1000 --dir N --out-file {output[1]} 2> {log}
        """

rule run_stars_null_total:
    input:  
        pert = "05norm/total_NORM_{comb}_pertubation.txt",
        lib = config["sgRNA"]["star_chip"]
    output: 
        "06stats/total_NORM_{comb}_Null_STARSOutput8_10_N.txt",
        "06stats/total_NORM_{comb}_Null_STARSOutput8_10_P.txt"
    log:    "00log/total_{comb}_stars_null"
    threads: 10
    params:
        mem  = "5G",
        name1 = get_name1,
        name2 = get_name2
    message: "run_stars_null_total {input}: {output} : {threads} threads" #"/ {params.mem}"
    shell: 
        """
        conda activate stars_env
        python {workflow.basedir}/STARS_v1.3/stars_null_v1.3.py --input-file {input.pert} --chip-file {input.lib} --thr 20 --num-ite 1000 --dir P --out-file {output[0]} 2> {log}
        python {workflow.basedir}/STARS_v1.3/stars_null_v1.3.py --input-file {input.pert} --chip-file {input.lib} --thr 20 --num-ite 1000 --dir N --out-file {output[1]} 2> {log}
        """

rule run_stars_control:
    input:  
        pert = "05norm/control_NORM_{comb}_pertubation.txt",
        lib = config["sgRNA"]["star_chip"],
        neg = "06stats/control_NORM_{comb}_Null_STARSOutput8_10_N.txt",
        pos = "06stats/control_NORM_{comb}_Null_STARSOutput8_10_P.txt"
    output: 
        # args.out_prefix + "_" + c +'_STARSOutput_'+direction+'.txt
        neg = "06stats/control_NORM_{comb}_Score_STARSOutput_N.txt",
        pos = "06stats/control_NORM_{comb}_Score_STARSOutput_P.txt"
    log:    "00log/control_NORM_{comb}_stars"
    threads: 2
    params:
        mem  = "2G",
        name1 = get_name1,
        name2 = get_name2,
        prefix = "06stats/control_NORM_{comb}"
    message: "run_stars_control {input}: {output} : {threads} threads" #"/ {params.mem}"
    shell: 
        """
        conda activate stars_env
        python {workflow.basedir}/STARS_v1.3/stars_v1.3.py --input-file {input.pert} --chip-file {input.lib} --dir P --thr 20 --null {input.pos} --use-first-pert N --out-prefix {params.prefix} 2> {log}
        python {workflow.basedir}/STARS_v1.3/stars_v1.3.py --input-file {input.pert} --chip-file {input.lib} --dir N --thr 20 --null {input.neg} --use-first-pert N --out-prefix {params.prefix} 2> {log}
        """

rule run_stars_total:
    input:  
        pert = "05norm/total_NORM_{comb}_pertubation.txt",
        lib = config["sgRNA"]["star_chip"],
        neg = "06stats/total_NORM_{comb}_Null_STARSOutput8_10_N.txt",
        pos = "06stats/total_NORM_{comb}_Null_STARSOutput8_10_P.txt"
    output: 
        # args.out_prefix + "_" + c +'_STARSOutput_'+direction+'.txt
        neg = "06stats/total_NORM_{comb}_Score_STARSOutput_N.txt",
        pos = "06stats/total_NORM_{comb}_Score_STARSOutput_P.txt"
    log:    "00log/total_NORM_{comb}_stars"
    threads: 2
    params:
        mem  = "2G",
        name1 = get_name1,
        name2 = get_name2,
        prefix = "06stats/total_NORM_{comb}"
    message: "run_stars_total {input}: {output} : {threads} threads" #"/ {params.mem}"
    shell: 
        """
        conda activate stars_env
        python {workflow.basedir}/STARS_v1.3/stars_v1.3.py --input-file {input.pert} --chip-file {input.lib} --dir P --thr 20 --null {input.pos} --use-first-pert N --out-prefix {params.prefix} 2> {log}
        python {workflow.basedir}/STARS_v1.3/stars_v1.3.py --input-file {input.pert} --chip-file {input.lib} --dir N --thr 20 --null {input.neg} --use-first-pert N --out-prefix {params.prefix} 2> {log}
        """


rule run_screen_beam:
    input:  "04quant/mageck.count_filtered.txt"
    output: "06stats/{comb}_ScreenBEAM_default.csv"
    log: "00log/{comb}_screen_beam"
    threads: 1
    params:
        mem  = "5G",
        name1 = get_name1,
        name2 = get_name2,
        sam_names1 = get_sam_names,
        sam_names2 = get_sam_names2
    message: "screen_beam {input}: {output} : {threads} threads"
    shell:
        """
        grep -v NTCONTROL {input} > {input}_{params.name1}_{params.name2}_ScreenBEAM_tmp
        Rscript {workflow.basedir}/scripts/run_screen_beam.R {input}_{params.name1}_{params.name2}_ScreenBEAM_tmp {output} {params.name1} {params.name2} {params.sam_names1} {params.sam_names2} 2> {log}
        rm {input}_{params.name1}_{params.name2}_ScreenBEAM_tmp
        """

rule run_screen_beam_control:
    input:  "05norm/control_{comb}.normalized.txt"
    output: "06stats/{comb}_ScreenBEAM_control.csv"
    log: "00log/{comb}_screen_beam_control"
    threads: 1
    params:
        mem  = "5G",
        name1 = get_name1,
        name2 = get_name2,
        sam_names1 = get_sam_names,
        sam_names2 = get_sam_names2
    message: "screen_beam_control {input}: {output} : {threads} threads"
    shell:
        """
        grep -v NTCONTROL {input} > {input}_{params.name1}_{params.name2}_ScreenBEAM_tmp_control
        Rscript {workflow.basedir}/scripts/run_screen_beam.R {input}_{params.name1}_{params.name2}_ScreenBEAM_tmp_control {output} {params.name1} {params.name2} {params.sam_names1} {params.sam_names2} 2> {log}
        rm {input}_{params.name1}_{params.name2}_ScreenBEAM_tmp_control
        """
rule run_screen_beam_total:
    input:  "04quant/mageck_total_filtered.count_normalized.txt"
    output: "06stats/{comb}_ScreenBEAM_total.csv"
    log: "00log/{comb}_screen_beam_total"
    threads: 1
    params:
        mem  = "5G",
        name1 = get_name1,
        name2 = get_name2,
        sam_names1 = get_sam_names,
        sam_names2 = get_sam_names2
    message: "screen_beam_total {input}: {output} : {threads} threads"
    shell:
        """
        grep -v NTCONTROL {input} > {input}_{params.name1}_{params.name2}_ScreenBEAM_tmp_total
        Rscript {workflow.basedir}/scripts/run_screen_beam.R {input}_{params.name1}_{params.name2}_ScreenBEAM_tmp_total {output} {params.name1} {params.name2} {params.sam_names1} {params.sam_names2} 2> {log}
        rm {input}_{params.name1}_{params.name2}_ScreenBEAM_tmp_total
        """

# ALL_PBNPA = expand("06stats/pbnpa_total_{comb}.csv 06stats/pbna_not_normalized_{comb}.csv".split() ,comb = ALL_COMB)
rule run_PBNPA_default:
    input:  "04quant/mageck.count.txt"
    output: "06stats/pbnpa_not_normalized_{comb}.csv"
    log: "00log/{comb}_pbnpa"
    threads: 1
    params:
        mem  = "5G",
        name1 = get_name1,
        name2 = get_name2,
        sam_names1 = get_sam_names,
        sam_names2 = get_sam_names2
    message: "PBNPA_default {input}: {output} : {threads} threads"
    shell:
        """
        Rscript {workflow.basedir}/scripts/run_pbnpa.R {input} {output} {params.name1} {params.name2} {params.sam_names1} {params.sam_names2} 2> {log}
        """

rule run_PBNPA_total:
    input:  "04quant/mageck_total_filtered.count_normalized.txt"
    output: "06stats/pbnpa_total_{comb}.csv"
    log: "00log/{comb}_pbnpa_total"
    threads: 1
    params:
        mem  = "5G",
        name1 = get_name1,
        name2 = get_name2,
        sam_names1 = get_sam_names,
        sam_names2 = get_sam_names2
    message: "PBNPA_total {input}: {output} : {threads} threads"
    shell:
        """
        Rscript {workflow.basedir}/scripts/run_pbnpa.R {input} {output} {params.name1} {params.name2} {params.sam_names1} {params.sam_names2} 2> {log}
        """

rule run_PBNPA_control:
    input:  "05norm/control_{comb}.normalized.txt"
    output: "06stats/pbnpa_control_{comb}.csv"
    log: "00log/{comb}_pbnpa_control"
    threads: 1
    params:
        mem  = "5G",
        name1 = get_name1,
        name2 = get_name2,
        sam_names1 = get_sam_names,
        sam_names2 = get_sam_names2
    message: "PBNPA_control {input}: {output} : {threads} threads"
    shell:
        """
        Rscript {workflow.basedir}/scripts/run_pbnpa.R {input} {output} {params.name1} {params.name2} {params.sam_names1} {params.sam_names2} 2> {log}
        """



rule run_mageck_default:
    input: 
        counts = "04quant/mageck.count.txt"
    output: "06stats/mageck_default_{comb}.gene_summary.txt"
    log: "00log/{comb}_mageck_default"
    threads: 2
    params:
        mem  = "5G",
        name1 = get_name1,
        name2 = get_name2,
        sam_names1 = get_sam_names,
        sam_names2 = get_sam_names2
    message: "run_mageck_default {input}: {output} : {threads} threads"
    shell:
        """
        conda activate /home/hayerk/miniconda3/envs/snakemake_new2/envs/mageck-vispr
        mageck test --norm-method median -k {input.counts} -t {params.sam_names2} -c {params.sam_names1} -n 06stats/mageck_default_{params.name1}_vs_{params.name2}  2> {log}
        """

rule run_mageck_control:
    input: 
        counts = "05norm/control_{comb}.normalized.txt"
    output: "06stats/mageck_control_{comb}.gene_summary.txt"
    log: "00log/{comb}_mageck_control"
    threads: 2
    params:
        mem  = "5G",
        name1 = get_name1,
        name2 = get_name2,
        sam_names1 = get_sam_names,
        sam_names2 = get_sam_names2,
        control_sgRNA = control_sgRNA
    message: "run_mageck_control {input}: {output} : {threads} threads"
    shell:
        """
        conda activate /home/hayerk/miniconda3/envs/snakemake_new2/envs/mageck-vispr
        mageck test --control-sgrna {params.control_sgRNA} --norm-method none -k {input.counts} -t {params.sam_names2} -c {params.sam_names1} -n 06stats/mageck_control_{params.name1}_vs_{params.name2}  2> {log}
        """

rule run_mageck_total:
    input: 
        counts = "04quant/mageck_total_filtered.count_normalized.txt"
    output: "06stats/mageck_total_{comb}.gene_summary.txt"
    log: "00log/{comb}_mageck_total"
    threads: 2
    params:
        mem  = "5G",
        name1 = get_name1,
        name2 = get_name2,
        sam_names1 = get_sam_names,
        sam_names2 = get_sam_names2
    message: "run_mageck_total {input}: {output} : {threads} threads"
    shell:
        """
        conda activate /home/hayerk/miniconda3/envs/snakemake_new2/envs/mageck-vispr
        mageck test --norm-method none -k {input.counts} -t {params.sam_names2} -c {params.sam_names1} -n 06stats/mageck_total_{params.name1}_vs_{params.name2}  2> {log}
        """

#def get_help(wildcards):
#    print(wildcards)
#    return ALL_COMB[wildcards[0]][1]

rule run_make_genelists_pbnpa:
    input:  "06stats/pbnpa_{norm}_{comb1}.csv", "06stats/pbnpa_{norm}_{comb2}.csv"
    output: "06stats/pbnpa_{norm}_NORM_{comb1}_COMPARED_{comb2}_genelist.csv", "07vis/pbnpa_{norm}_NORM_{comb1}_COMPARED_{comb2}_negative.pdf", "07vis/pbnpa_{norm}_NORM_{comb1}_COMPARED_{comb2}_positve.pdf"
    log:    "00log/make_genelists_pbnpa_{norm}_{comb1}_{comb2}"
    threads: 2
    params:
        mem  = "5G",
        essential_genes = essential_genes,
        norm_method = "{norm}",
        comb2 = "{comb2}",
        all_genes = config["all_genes"],
        genes = genes

    message: "run_make_genelists_pbnpa {input}: {output} : {threads} threads" #"/ {params.mem}"
    shell: 
        """
        Rscript {workflow.basedir}/scripts/run_pbnpa_genelist.R {input[0]} {input[1]} {output[0]} {params.norm_method} {output[1]} {output[2]} {params.essential_genes} {params.all_genes} {params.genes}
        """

rule run_make_genelists_mageck:
    input:  "06stats/mageck_{norm}_{comb1}.gene_summary.txt", "06stats/mageck_{norm}_{comb2}.gene_summary.txt"
    output: "06stats/mageck_{norm}_NORM_{comb1}_COMPARED_{comb2}_genelist.csv", "07vis/mageck_{norm}_NORM_{comb1}_COMPARED_{comb2}_negative.pdf", "07vis/mageck_{norm}_NORM_{comb1}_COMPARED_{comb2}_positve.pdf"
    log:    "00log/make_genelists_mageck_{norm}_{comb1}_{comb2}"
    threads: 2
    params:
        mem  = "5G",
        essential_genes = essential_genes,
        norm_method = "{norm}",
        comb2 = "{comb2}",
        all_genes = config["all_genes"],
        genes = genes
        #nina = get_help,

    message: "run_make_genelists_mageck {input}: {output} : {threads} threads" #"/ {params.mem}"
    shell: 
        """
        Rscript {workflow.basedir}/scripts/run_mageck_genelist.R {input[0]} {input[1]} {output[0]} {params.norm_method} {output[1]} {output[2]} {params.essential_genes} {params.all_genes} {params.genes}
        """

rule run_make_genelists_screenbeam:
    input:  "06stats/{comb1}_ScreenBEAM_{norm}.csv", "06stats/{comb2}_ScreenBEAM_{norm}.csv"
    output: "06stats/ScreenBEAM_{norm}_NORM_{comb1}_COMPARED_{comb2}_genelist.csv", "07vis/ScreenBEAM_{norm}_NORM_{comb1}_COMPARED_{comb2}_negative.pdf", "07vis/ScreenBEAM_{norm}_NORM_{comb1}_COMPARED_{comb2}_positve.pdf"
    log:    "00log/make_genelists_screenbeam_{norm}_{comb1}_{comb2}"
    threads: 2
    params:
        mem  = "5G",
        essential_genes = essential_genes,
        norm_method = "{norm}",
        comb2 = "{comb2}",
        all_genes = config["all_genes"],
        genes = genes
        #nina = get_help,

    message: "run_make_genelists_screenbeam {input}: {output} : {threads} threads" #"/ {params.mem}"
    shell: 
        """
        Rscript {workflow.basedir}/scripts/run_screenbeam_genelist.R {input[0]} {input[1]} {output[0]} {params.norm_method} {output[1]} {output[2]} {params.essential_genes} {params.all_genes} {params.genes}
        """

rule run_make_genelists_stars:
    input:  
        "06stats/{norm}_NORM_{comb1}_Score_STARSOutput_P.txt",  # P is negative
        "06stats/{norm}_NORM_{comb1}_Score_STARSOutput_N.txt",  # N is positve
        "06stats/{norm}_NORM_{comb2}_Score_STARSOutput_P.txt", 
        "06stats/{norm}_NORM_{comb2}_Score_STARSOutput_N.txt"
    output: "06stats/STARS_{norm}_NORM_{comb1}_COMPARED_{comb2}_genelist.csv", "07vis/STARS_{norm}_NORM_{comb1}_COMPARED_{comb2}_negative.pdf", "07vis/STARS_{norm}_NORM_{comb1}_COMPARED_{comb2}_positve.pdf"
    log:    "00log/make_genelists_stars_control_{norm}_{comb1}_{comb2}"
    threads: 2
    params:
        mem  = "5G",
        essential_genes = essential_genes,
        norm_method = "{norm}",
        comb2 = "{comb2}",
        all_genes = config["all_genes"],
        genes = genes
        #nina = get_help,

    message: "run_make_genelists_stars {input}: {output} : {threads} threads" #"/ {params.mem}"
    shell: 
        """
        Rscript {workflow.basedir}/scripts/run_stars_genelist.R {input[0]} {input[1]} {input[2]} {input[3]} {output[0]} {params.norm_method} {output[1]} {output[2]} {params.essential_genes} {params.all_genes} {params.genes}
        """

