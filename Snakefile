import os
import re
import yaml
###############################################################################
# Load config
###############################################################################
with open("config.yaml", "r") as f:
    config = yaml.safe_load(f)
TREE_TOPOLOGY = config["tree_topology"]
TREE_TOPOLOGY_NAMED = config.get("tree_topology_named", TREE_TOPOLOGY)  # USES unnamed if not defined, or commented out in config.yaml
FOREGROUND_BRANCHES = config["foreground_branches"]
CHROMOSOMES = config["chromosomes"]
NUM_REPLICATES = config["num_replicates"]
MIN_FRAC = config["min_frac"]
###############################################################################
# Input functions — depend on the final list of good .fa files
###############################################################################
def get_query_fa_for_chrom(wildcards):
    ck_out = checkpoints.collect_good_alignments.get(**wildcards).output[0]
    with open(ck_out) as f:
        paths = [line.strip() for line in f if line.strip()]
    pattern = rf"^{wildcards.chrom}\."
    matching = [p for p in paths if re.search(pattern, os.path.basename(p))]
    return matching
def get_ref_fa_for_chrom(wildcards):
    ck_out = checkpoints.extract_ref_basenames.get(**wildcards).output[0]
    with open(ck_out) as f:
        refs = [line.strip() for line in f if line.strip()]
    pattern = rf"^{wildcards.chrom}\."
    matching = [r for r in refs if re.search(pattern, r)]
    return expand("ref/{ref}.{rep}.ref", ref=matching, rep=range(NUM_REPLICATES))
def get_bf_for_chrom(wildcards):
    ck_out = checkpoints.make_bf_list.get(**wildcards).output[0]
    with open(ck_out) as f:
        bfs = [line.strip() for line in f if line.strip()]
    pattern = rf"^{wildcards.chrom}\."
    return [bf for bf in bfs if re.search(pattern, os.path.basename(bf))]
# Rules
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
        mkdir -p features logs
        grep -w {wildcards.chrom} {input.bed} | \
            awk '{{printf "%s\\t%d\\t%s\\n", $1, ($2-1), $3}}' | \
            sort -k1,1 -k2,2n > {output} 2> {log} || true
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
        mkdir -p query/{wildcards.chrom} logs
        msa_split {input.maf} --refseq {input.fa} --gap-strip ANY -q \
                  --features {input.feat} --for-features \
                  --out-root query/{wildcards.chrom}/{wildcards.chrom} > {log} 2>&1 || true
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
        mkdir -p good_alignments logs
        touch {log}
        if [ -d "query/{wildcards.chrom}" ]; then
            python3 scripts/select_and_filter.py \
                --input query/{wildcards.chrom} \
                --min-frac {params.min_frac} \
                --out good_alignments >> {log} 2>&1 || echo "[WARNING] select_and_filter failed for {wildcards.chrom}, continuing anyway" >> {log}
        else
            echo "[INFO] query/{wildcards.chrom} does not exist — nothing to filter (normal for small/no-window chromosomes)" >> {log}
        fi
        touch {output}
        """
# 5. CHECKPOINT AND RULE TO FILTER OUT ALIGNMENTS WITH MISSING SEQUENCE (e.g. INDELS, MISSALIGNMENTS)
rule collect_good_per_chrom:
    input:
        "filter/{chrom}.done"
    output:
        "goodalignments_{chrom}.txt"
    resources: mem_mb = 2000
    shell:
        """
        mkdir -p good_alignments
        find good_alignments -name "{wildcards.chrom}.*.fa" -type f | sort > {output} || touch {output}
        total=$(wc -l < {output})
        echo "Found $total good alignments for {wildcards.chrom}" >&2
        [[ $total -gt 0 ]] || echo "WARNING: No good alignments for {wildcards.chrom}" >&2
        """
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
# 6. DictGen: This makes a dictionary to create a concatenated reference for each query from neutral proxies from entire genome or local region

rule dict_gen_chrom:
    message: "!!!!!NOTE!!!!! if neutral == goodalignments_{wildcards.chrom}.txt then it is a local test for {wildcards.chrom}"
    input:
        good = "goodalignments_{chrom}.txt",
        neutral = config["neutral_set"]
    output:
        dict = "global_{chrom}.dict",
        ref_list = "reference_{chrom}.list",
        done = touch("ref/done_{chrom}")
    conda:
        "envs/biopython.yaml"
    log:
        "logs/dict_gen_{chrom}.log"
    resources:
        mem_mb = 10000
    shell:
        """
        mkdir -p ref logs
        touch {log} {output.dict} {output.ref_list} {output.done}
        if [[ -s {input.good} ]]; then
            python scripts/DictGen.py {input.good} {input.neutral} ref {output.dict} {output.ref_list} >> {log} 2>&1 || echo "WARNING: DictGen.py failed for {wildcards.chrom}, continuing" >> {log}
            count=$(find ref -name "{wildcards.chrom}.*.ref" | wc -l)
            echo "DictGen created $count reference alignments for {wildcards.chrom}" >> {log}
            [[ $count -gt 0 ]] || echo "WARNING: No .ref files created for {wildcards.chrom} (possibly empty or invalid input)" >> {log}
        else
            echo "[INFO] Empty goodalignments_{wildcards.chrom}.txt - skipping DictGen, no refs for {wildcards.chrom}" >> {log}
        fi
        true
        """


rule aggregate_dict_gen:
    input:
        expand("ref/done_{chrom}", chrom=CHROMOSOMES),
        expand("goodalignments_{chrom}.txt", chrom=CHROMOSOMES)
    output:
        ref_done = touch("ref/.dict_gen_complete"),
        ref_list = "reference.list"
    resources: mem_mb = 2000
    shell:
        """
        mkdir -p ref
        cat reference_*.list 2>/dev/null > {output.ref_list} || touch {output.ref_list}
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
        find ref -name "*.ref" | sed -E 's/(.*\\/)|\\.[0-9]+\\.ref$//g' | sort -u > {output} || touch {output}
        """
# Adaptiphy branch tests (bf files creation to execution)
# 7. Generates batch files for HYPHY input
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
        mkdir -p logs
        python {input.script} {input.ref_list} > {log} 2>&1 || true
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
# 8. Runs Adaptiphy on each query and respective replicate reference, for the null and the alternative hypothesis
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
        mkdir -p HYPHY res logs
        touch {log}
        if [ "{input.bfs}" != "" ]; then
            for bf in {input.bfs}; do
                name=$(basename "$bf" .bf)
                hyphy "$bf" > "HYPHY/$name.out" 2>> {log} || true
            done
        else
            echo "[INFO] No bfs for {wildcards.chrom} - skipping hyphy" >> {log}
        fi
        """
# It's a type of checkpoint to make sure adaptiphy has finished running
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
        mkdir -p logs
        python scripts/extract_res.py > {output} 2> {log} || touch {output}
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
        mkdir -p MODELS_HKY85/query logs
        touch {log}
        if [ "{input.fas}" != "" ]; then
            for fa in {input.fas}; do
                base=$(basename "$fa" .fa)
                phyloFit "$fa" --tree "{TREE_TOPOLOGY_NAMED}" -i FASTA --subst-mod HKY85 --init-random --precision HIGH \
                         --out-root MODELS_HKY85/query/"$base" >> {log} 2>&1 || echo "WARNING: phyloFit failed for $fa" >> {log}
            done
        else
            echo "[INFO] No query alignments for {wildcards.chrom} - skipping phyloFit" >> {log}
        fi
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
        mkdir -p MODELS_HKY85/ref logs
        touch {log}
        if [ "{input.fas}" != "" ]; then
            for fa in {input.fas}; do
                ref_rep=$(basename "$fa" .ref)
                phyloFit "$fa" --tree "{TREE_TOPOLOGY_NAMED}" -i FASTA --subst-mod HKY85 --init-random --precision HIGH \
                         --out-root MODELS_HKY85/ref/"$ref_rep" >> {log} 2>&1 || echo "WARNING: phyloFit failed for $fa" >> {log}
            done
        else
            echo "[INFO] No ref alignments for {wildcards.chrom} - skipping phyloFit" >> {log}
        fi
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
# 11. This rule computes zeta values, which is the evolutionary ratio, which is the relative mutation rate. Similar to omega
rule calculate_zeta_chrom:
    wildcard_constraints:
        chrom = "|".join(CHROMOSOMES)
    input:
        ref_list = "reference.list",
        query_mods = "query_mods.done",
        ref_mods = "ref_mods.done"
    output:
        "summary_table_branch_{chrom}.txt"
    conda:
        "envs/biopython.yaml"
    log:
        "logs/calculate_zeta_{chrom}.log"
    resources:
        mem_mb = 2000
    shell:
        """
        mkdir -p logs
        python scripts/calculate_zeta.py {input.ref_list} MODELS_HKY85/query MODELS_HKY85/ref {output} --chrom {wildcards.chrom} --log {log} || touch {output}
        """

rule aggregate_zeta:
    input:
        expand("summary_table_branch_{chrom}.txt", chrom=CHROMOSOMES)
    output:
        "summary_table_branch.txt"
    resources:
        mem_mb = 2000
    shell:
        """
        touch {output}
        first_file=$(find summary_table_branch_*.txt -size +0c | head -1)
        if [[ -n "$first_file" ]]; then
            cat "$first_file" > {output}  # Copy first non-empty file (with header)
            for f in summary_table_branch_*.txt; do
                if [[ "$f" != "$first_file" ]] && [[ -s "$f" ]]; then
                    tail -n +2 "$f" >> {output}  # Append data without header
                fi
            done
        else
            echo "No non-empty zeta files - creating empty output" >&2
        fi
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
try:
    a = pd.read_csv('{input[0]}', sep='\\t')
    b = pd.read_csv('{input[1]}', sep='\\t')
    if a.empty or b.empty:
        pd.DataFrame().to_csv('{output}', sep='\\t', index=False)
    else:
        m = a.merge(b, on=['NAME','BRANCH','REPL'], how='inner')
        m.to_csv('{output}', sep='\\t', index=False)
except pd.errors.EmptyDataError:
    print('One of the input files is empty - creating empty output')
    pd.DataFrame().to_csv('{output}', sep='\\t', index=False)
" || touch {output}
        """
##############################################################
################# FINAL CLEAN UP ############################
rule cleanup:
    input: "merged_summary_table.txt"
    output: "ADAPTIPHY_DONE"
    resources:
        mem_mb = 2000
    shell:
        """
        mkdir -p OUTPUT_FINAL
        rm -r *.done bf_list.txt query bf_dir ref ref_basenames.txt reference.list goodalignments.txt features filter summary_table_branch_*.txt || true
        rm -r global_*.dict reference_*.list goodalignments_*.txt || true
        mv HYPHY MODELS_HKY85 good_alignments res logs summary_table*.txt merged_summary_table.txt OUTPUT_FINAL/ 2>/dev/null || true
        touch {output}
        """
