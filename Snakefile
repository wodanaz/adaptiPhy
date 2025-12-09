import os
import re
import yaml

# Load configuration

with open("config.yaml", "r") as f:
    config = yaml.safe_load(f)

TREE_TOPOLOGY = config["tree_topology"]
FOREGROUND_BRANCHES = config["foreground_branches"]
CHROMOSOMES = config["chromosomes"]
NUM_REPLICATES = config["num_replicates"]
MIN_FRAC = config["min_frac"]

# Helper functions

def get_query_fa_for_chrom(wildcards):
    ck_out = checkpoints.collect_good_alignments.get(**wildcards).output[0]
    with open(ck_out) as f:
        paths = [line.strip() for line in f if line.strip()]
    pattern = rf"\b{wildcards.chrom}\b"
    matching = [p for p in paths if re.search(pattern, os.path.basename(p))]
    return matching

def get_ref_fa_for_chrom(wildcards):
    ck_out = checkpoints.extract_ref_basenames.get(**wildcards).output[0]
    with open(ck_out) as f:
        refs = [line.strip() for line in f if line.strip()]
    pattern = rf"\b{wildcards.chrom}\b"
    matching = [r for r in refs if re.search(pattern, r)]
    return expand("ref/{ref}.{rep}.ref", ref=matching, rep=range(NUM_REPLICATES))

def get_bf_for_chrom(wildcards):
    ck_out = checkpoints.make_bf_list.get(**wildcards).output[0]
    with open(ck_out) as f:
        bfs = [line.strip() for line in f if line.strip()]
    pattern = rf"\b{wildcards.chrom}\b"
    return [bf for bf in bfs if re.search(pattern, os.path.basename(bf))]

# MAIN RULE : TARGET 
rule all:
    input: "ADAPTIPHY_DONE"

# 1. Setting up the hyphy working directory 

rule init_dirs:
    output:
        "bf_dir/alt4-fgrnd_spec.model",
        "bf_dir/null4-fgrnd_spec.model"
    shell:
        "mkdir -p bf_dir logs && cp scripts/*model bf_dir/"

# 2. split main query bed file into chromosome form
rule features:
    input: bed = config["windows"]
    output: "features/{chrom}.feat.bed"
    log: "logs/features_{chrom}.log"
    resources: mem_mb = 2000
    shell:
        """
        mkdir -p features
        grep -w {wildcards.chrom} {input.bed} | \
            awk '{{printf "%s\\t%d\\t%s\\n", $1, ($2-1), $3}}' | \
            sort -k1,1 -k2,2n > {output} 2> {log}
        [[ -s {output} ]] || touch {output}
        """

# 3. Split genome-wide alignment in MAF format into small fasta files
rule query:
    input:
        maf = config["maf_pattern"],
        fa = config["fa_pattern"],
        feat = "features/{chrom}.feat.bed"
    output:
        "query/{chrom}/done.txt"
    conda: "envs/biopython.yaml"
    log: "logs/query_{chrom}.log"
    resources: mem_mb = 15000
    shell:
        """
        mkdir -p query/{wildcards.chrom}
        msa_split {input.maf} --refseq {input.fa} --gap-strip ANY -q \
                  --features {input.feat} --for-features \
                  --out-root query/{wildcards.chrom}/{wildcards.chrom} > {log} 2>&1
        touch {output}
        """

# 4. FILTER GOOD ALIGNMENTS PER CHROMOSOME
rule filter_good_per_chrom:
    wildcard_constraints:
        chrom = "|".join(CHROMOSOMES)
    input:
        "query/{chrom}/done.txt"
    output:
        touch("filter/{chrom}.done")
    params:
        min_frac = MIN_FRAC
    conda: "envs/biopython.yaml"
    log: "logs/filter_{chrom}.log"
    resources: mem_mb = 10000
    shell:
        """
        # Create output dir
        mkdir -p good_alignments

        if [ -d "query/{wildcards.chrom}" ]; then
            python3 scripts/select_and_filter.py \
                --input query/{wildcards.chrom} \
                --min-frac {params.min_frac} \
                --out good_alignments > {log} 2>&1 || echo "[WARNING] select_and_filter failed for {wildcards.chrom}, continuing anyway" >> {log}
        else
            echo "[INFO] query/{wildcards.chrom} does not exist — nothing to filter (normal for small/no-window chromosomes)" > {log}
        fi

        # Always touch the done file
        touch {output}
        """

# 5. CHECKPOINT AND RULE TO FILTER OUT ALIGNMENTS WITH MISSING SEQUENCE (e.g. INDELS, MISSALIGNMENTS)

checkpoint collect_good_alignments:
    input:
        expand("filter/{chrom}.done", chrom=CHROMOSOMES)
    output:
        "goodalignments.txt"
    resources: mem_mb = 2000
    shell:
        """
        mkdir -p good_alignments
        find good_alignments -name "*.fa" -type f | sort > {output}
        total=$(wc -l < {output})
        echo "Found $total good alignments total" >&2
        [[ $total -gt 0 ]] || echo "WARNING: No good alignments found!" >&2
        touch {output} # always exists
        """
# 6.  DictGen: This makes a dictionary to create a concatenated reference for each query from neutral proxies from entire genome or local region
rule dict_gen:
    message: "!!!!!NOTE!!!!! if neutral == goodalignments.txt then it is a local test"
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
        mem_mb = 10000
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
        mem_mb = 2000
    shell:
        """
        find ref -name "*.ref" | sed -E 's/(.*\\/)|\\.[0-9]+\\.ref$//g' | sort -u > {output}
        """

# Adaptiphy branch tests (bf files creation to execution)
# 7. Generates batch files for  HYPHY input
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
        mem_mb = 2000
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
        mem_mb = 2000
    shell:
        """
        find bf_dir -name "*.bf" > {output} || touch {output}
        """
# 8. Runs Adaptiphy  on each query and respective replicate reference, for the null and the alternative hypothesis
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
        mem_mb = 10000
    shell:
        """
        mkdir -p HYPHY res
        for bf in {input.bfs}; do
            name=$(basename "$bf" .bf)
            hyphy "$bf" > "HYPHY/$name.out" 2>> {log} || true
        done
        """
# It's a type of checkpoint to maje sure adaptiphy has finished running
rule aggregate_hyphy_outputs:
    input:
        expand("HYPHY/{chrom}.done", chrom=CHROMOSOMES)
    output:
        touch("hyphy_outputs.done")
    resources:
        mem_mb = 2000
    shell:
        "echo 'HyPhy branch-site tests complete'"


# 9. Aggregates adaptiphy data using python function from result files
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
        mem_mb = 2000
    shell:
        """
        python scripts/extract_res.py > {output} 2> {log}
        """

# 10. Run phyloFit on query and reference to produce evolutionary model for the query and each replicated reference. It contain rules that act as checkpoints

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
        mem_mb = 5000
    shell:
        """
        for fa in {input.fas}; do
            base=$(basename "$fa" .fa)
            phyloFit "$fa" --tree "{TREE_TOPOLOGY}" -i FASTA --subst-mod HKY85 --init-random --precision HIGH \
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
        mem_mb = 5000
    shell:
        """
        for fa in {input.fas}; do
            ref_rep=$(basename "$fa" .ref)
            phyloFit "$fa" --tree "{TREE_TOPOLOGY}" -i FASTA --subst-mod HKY85 --init-random --precision HIGH \
                     --out-root MODELS_HKY85/ref/"$ref_rep" >> {log} 2>&1
        done
        """

rule aggregate_query_mods:
    input:
        expand("MODELS_HKY85/query/{chrom}.done", chrom=CHROMOSOMES)
    output:
        touch("query_mods.done")
    resources:
        mem_mb = 2000
    shell:
        "echo 'Query phyloFit done'"

rule aggregate_refs:
    input:
        expand("MODELS_HKY85/ref/{chrom}.done", chrom=CHROMOSOMES)
    output:
        touch("ref_mods.done")
    resources:
        mem_mb = 2000
    shell:
        "echo 'Reference phyloFit done'"


# 11. This function computes zeta values, which is the evolutionary ratio, which is the relative mutation rate. Similar to omega
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
        mem_mb = 2000
    shell:
        """
        python scripts/calculate_zeta.py {input.ref_list} MODELS_HKY85/query MODELS_HKY85/ref {output} --log {log}
        """
# 12. merges tabular results into a central table that can be processed in R
rule merge_summaries:
    input:
        "summary_table_branch.txt",
        "summary_table.txt"
    output:
        "merged_summary_table.txt"
    conda:
        "envs/biopython.yaml"
    resources:
        mem_mb = 2000
    shell:
        """
        python -c "import pandas as pd;
a = pd.read_csv('{input[0]}', sep='\t');
b = pd.read_csv('{input[1]}', sep='\t');
m = a.merge(b, on=['NAME','BRANCH','REPL'], how='left');
m.to_csv('{output}', sep='\t', index=False)"
        """

##############################################################
################# FINAL CLEAN UP  ############################
rule cleanup:
    input: "merged_summary_table.txt"
    output: "ADAPTIPHY_DONE"
    resources:
        mem_mb = 2000
    shell:
        """
        mkdir -p OUTPUT_FINAL
        rm -r HYPHY *.done bf_list.txt bf_dir ref_basenames.txt reference.list global.dict goodalignments.txt features filter || true
        mv MODELS_HKY85 query good_alignments ref res logs summary_table*.txt merged_summary_table.txt OUTPUT_FINAL/ 2>/dev/null || true
        touch {output}
        """
