import os
import yaml
import re

# Load config
with open("config.yaml", "r") as f:
    config = yaml.safe_load(f)
TREE_TOPOLOGY = config["tree_topology"]
FOREGROUND_BRANCHES = config["foreground_branches"]
CHROMOSOMES = config["chromosomes"]
NUM_REPLICATES = config["num_replicates"]
MIN_BASES     = config["min_bases"]
MAX_KEEP      = config["max_keep"]

# Query basenames from select_good checkpoint
def get_query_base_names(wildcards):
    ck_out = checkpoints.select_good.get(**wildcards).output.goodlist
    with open(ck_out) as f:
        return [os.path.basename(l.strip()).rsplit(".", 1)[0] for l in f if l.strip()]

# Reference basenames from what DictGen actually created
def get_ref_base_names(wildcards):
    ck_out = checkpoints.extract_ref_basenames.get(**wildcards).output[0]
    with open(ck_out) as f:
        return [l.strip() for l in f if l.strip()]

def get_bf_names(wildcards):
    ck_out = checkpoints.make_bf_list.get(**wildcards).output[0]
    with open(ck_out) as f:
        return [l.strip().replace("bf_dir/", "").replace(".bf", "") for l in f if l.strip()]

def get_bf_for_chrom(wildcards):
    ck_out = checkpoints.make_bf_list.get(**wildcards).output[0]
    with open(ck_out) as f:
        all_bf = [l.strip() for l in f if l.strip()]
    return [bf for bf in all_bf if re.search(rf"\b{wildcards.chrom}\b", os.path.basename(bf))]

def get_query_fa_for_chrom(wildcards):
    ck_out = checkpoints.select_good.get(**wildcards).output.goodlist
    with open(ck_out) as f:
        all_bases = [os.path.basename(l.strip()).rsplit(".", 1)[0] for l in f if l.strip()]
    chrom_bases = [b for b in all_bases if re.search(rf"\b{wildcards.chrom}\b", b)]
    return [f"good_alignments/{b}.fa" for b in chrom_bases]

def get_ref_fa_for_chrom(wildcards):
    ck_out = checkpoints.extract_ref_basenames.get(**wildcards).output[0]
    with open(ck_out) as f:
        all_refs = [l.strip() for l in f if l.strip()]
    chrom_refs = [r for r in all_refs if re.search(rf"\b{wildcards.chrom}\b", r)]
    return expand("ref/{ref}.{replicate}.ref", ref=chrom_refs, replicate=range(NUM_REPLICATES))

####################################################################################################
rule all:
    input:
        "DONE_SUMMARY.txt"

rule init_dirs:
    output:
        "bf_dir/alt4-fgrnd_spec.model",
        "bf_dir/null4-fgrnd_spec.model"
    resources:
        mem_mb = 500
    shell:
        """
        mkdir -p bf_dir/ logs/
        cp scripts/null4-fgrnd_spec.model scripts/alt4-fgrnd_spec.model bf_dir/
        """

rule features:
    input:
        bed = config["windows"]
    output:
        "features/{chrom}.feat.bed"
    log:
        "logs/features_{chrom}.log"
    resources:
        mem_mb = 1000
    shell:
        """
        mkdir -p features
        grep -w {wildcards.chrom} {input.bed} | awk '{{printf "%s\\t%d\\t%s\\n", $1, ($2-1), $3}}' | sort -k1,1 -k2,2n > {output} 2> {log}
        [[ -s {output} ]] || touch {output}
        """

rule query:
    input:
        maf = config["maf_pattern"],
        fa = config["fa_pattern"],
        feat = "features/{chrom}.feat.bed"
    output:
        "query/{chrom}/done.txt"
    conda:
        "envs/biopython.yaml"
    log:
        "logs/query_{chrom}.log"
    resources:
        mem_mb = 10000
    shell:
        """
        mkdir -p query/{wildcards.chrom}
        msa_split {input.maf} --refseq {input.fa} --gap-strip ANY -q \
                  --features {input.feat} --for-features \
                  --out-root query/{wildcards.chrom} > {log} 2>&1
        touch {output}
        """

# =============================================================================
# Filtering good alignments, this creates a new directory with good alignments
# =============================================================================

checkpoint select_good:
    input:
        expand("query/{chrom}/done.txt", chrom=CHROMOSOMES)
    output:
        goodlist = "goodalignments.txt",
        dir = directory("good_alignments")
    params:
        min_bases = MIN_BASES,
        max_keep  = MAX_KEEP
    conda:
        "envs/biopython.yaml"
    log:
        "logs/select_and_filter.log"
    resources:
        mem_mb = 500
    shell:
        """
        mkdir -p good_alignments 
        python3 scripts/select_and_filter.py \
            --input query \
            --min-bases {params.min_bases} \
            --max-keep {params.max_keep} \
            --out good_alignments \
            --good-list {output.goodlist} > {log} 2>&1 || touch {output.goodlist}
        """

# =======================================================================================
# DictGen + completion flag: This uses a dictionary to create a reference for each query
# =======================================================================================
rule dict_gen:
    message: "!!!!!NOTE!!!!! if neutral == goodalignments then it is a local test"
    input:
        good = "goodalignments.txt",
        neutral = config["neutral_set"]
    output:
        dict = "global.dict",
        ref_list = "reference.list",
        ref_dir = directory("ref"),
        done = touch("ref/.dict_gen_complete")
    conda:
        "envs/biopython.yaml"
    log:
        "logs/dict_gen.log"
    resources:
        mem_mb = 1000
    shell:
        """
        mkdir -p ref
        python scripts/DictGen.py {input.neutral} > {log} 2>&1
        count=$(find ref -name "*.ref" | wc -l)
        echo "DictGen created $count reference alignments" >> {log}
        [[ $count -gt 0 ]] || {{ echo "ERROR: No .ref files!" >> {log}; exit 1; }}
        """

checkpoint extract_ref_basenames:
    input:
        "ref/.dict_gen_complete"
    output:
        "ref_basenames.txt"
    resources:
        mem_mb = 200
    shell:
        """
        find ref -name "*.ref" | sed -E 's/(.*\\/)|\\.[0-9]+\\.ref$//g' | sort -u > {output}
        """

# =============================================================================
# Adaptiphy branch tests (bf files creation to execution)
# =============================================================================


rule bfgenerator:
    input:
        ref_list = "reference.list",
        script = "scripts/bf_generator.py",
        model_alt = "bf_dir/alt4-fgrnd_spec.model",
        model_null = "bf_dir/null4-fgrnd_spec.model"
    output:
        touch("bf_dir/bf_generated.flag")
    conda:
        "envs/biopython.yaml"
    log:
        "logs/bfgenerator.log"
    resources:
        mem_mb = 500
    shell:
        """
        python {input.script} {input.ref_list} > {log} 2>&1
        """

checkpoint make_bf_list:
    input:
        "bf_dir/bf_generated.flag"
    output:
        "bf_list.txt"
    resources:
        mem_mb = 500
    shell:
        """
        find bf_dir -name "*.bf" > {output} || touch {output}
        """

rule run_hyphy:
    wildcard_constraints:
        chrom = "|".join(CHROMOSOMES)
    input:
        bfs = get_bf_for_chrom
    output:
        touch("HYPHY/{chrom}.done")
    conda:
        "envs/biopython.yaml"
    log:
        "logs/hyphy_{chrom}.log"
    resources:
        mem_mb = 250
    shell:
        """
        mkdir -p HYPHY res
        for bf in {input.bfs}; do
            name=$(basename "$bf" .bf)
            hyphy "$bf" > "HYPHY/$name.out" 2>> {log} || true
        done
        """

rule aggregate_hyphy_outputs:
    input:
        expand("HYPHY/{chrom}.done", chrom=CHROMOSOMES)
    output:
        touch("hyphy_outputs.done")
    resources:
        mem_mb = 500
    shell:
        "echo 'HyPhy branch-site tests complete'"

rule extract_results:
    input:
        "hyphy_outputs.done"
    output:
        "summary_table.txt"
    conda:
        "envs/biopython.yaml"
    log:
        "logs/extract_results.log"
    resources:
        mem_mb = 500
    shell:
        """
        python scripts/extract_res.py > {output} 2> {log}
        """

# =============================================================================
# phyloFit on query and reference
# =============================================================================
rule run_phylofit_query:
    wildcard_constraints:
        chrom = "|".join(CHROMOSOMES)
    input:
        fas = get_query_fa_for_chrom
    output:
        touch("MODELS_HKY85/query/{chrom}.done")
    conda:
        "envs/biopython.yaml"
    log:
        "logs/phylofit_query_{chrom}.log"
    resources:
        mem_mb = 400
    shell:
        """
        for fa in {input.fas}; do
            base=$(basename "$fa" .fa)
            phyloFit "$fa" --tree "{TREE_TOPOLOGY}" -i FASTA --subst-mod HKY85 \
                     --out-root MODELS_HKY85/query/"$base" >> {log} 2>&1
        done
        """

rule run_phylofit_ref:
    wildcard_constraints:
        chrom = "|".join(CHROMOSOMES)
    input:
        fas = get_ref_fa_for_chrom,
        done = "ref/.dict_gen_complete"
    output:
        touch("MODELS_HKY85/ref/{chrom}.done")
    conda:
        "envs/biopython.yaml"
    log:
        "logs/phylofit_ref_{chrom}.log"
    resources:
        mem_mb = 400
    shell:
        """
        for fa in {input.fas}; do
            ref_rep=$(basename "$fa" .ref)
            phyloFit "$fa" --tree "{TREE_TOPOLOGY}" -i FASTA --subst-mod HKY85 \
                     --out-root MODELS_HKY85/ref/"$ref_rep" >> {log} 2>&1
        done
        """

rule aggregate_query_mods:
    input:
        expand("MODELS_HKY85/query/{chrom}.done", chrom=CHROMOSOMES)
    output:
        touch("query_mods.done")
    resources:
        mem_mb = 200
    shell:
        "echo 'Query phyloFit done'"

rule aggregate_refs:
    input:
        expand("MODELS_HKY85/ref/{chrom}.done", chrom=CHROMOSOMES)
    output:
        touch("ref_mods.done")
    resources:
        mem_mb = 200
    shell:
        "echo 'Reference phyloFit done'"

rule calculate_zeta:
    input:
        ref_list = "reference.list",
        query_mods = "query_mods.done",
        ref_mods = "ref_mods.done"
    output:
        "summary_table_branch.txt"
    conda:
        "envs/biopython.yaml"
    log:
        "logs/calculate_zeta.log"
    resources:
        mem_mb = 500
    shell:
        """
        python scripts/calculate_zeta.py {input.ref_list} MODELS_HKY85/query MODELS_HKY85/ref {output} --log {log}
        """

rule merge_summaries:
    input:
        "summary_table_branch.txt",
        "summary_table.txt"
    output:
        "merged_summary_table.txt"
    conda:
        "envs/biopython.yaml"
    resources:
        mem_mb = 500
    shell:
        """
        python -c "import pandas as pd;
a = pd.read_csv('{input[0]}', sep='\t');
b = pd.read_csv('{input[1]}', sep='\t');
m = a.merge(b, on=['NAME','BRANCH','REPL'], how='left');
m.to_csv('{output}', sep='\t', index=False)"
        """

##############################################################
## FINAL CLEAN UP
#################



rule cleanup:
    input: "merged_summary_table.txt"
    output: "ADAPTIPHY_DONE"
    resources:
        mem_mb = 500
    shell:
        """
        mkdir -p OUTPUT_FINAL
        rm -r HYPHY *.done bf_list.txt bf_dir ref_basenames.txt reference.list global.dict goodalignments.txt selected_alignments features || true
        mv MODELS_HKY85 query good_alignments ref res logs summary_table*.txt merged_summary_table.txt OUTPUT_FINAL/ 2>/dev/null || true
        touch {output}
        """
