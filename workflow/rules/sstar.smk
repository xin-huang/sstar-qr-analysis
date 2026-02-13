# Copyright 2026 Xin Huang and Andrea Koca
#
# GNU General Public License v3.0
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, please see
#
#    https://www.gnu.org/licenses/gpl-3.0.en.html


ruleorder: sstar_threshold > run_qr
ruleorder: get_sstar_inferred_tracts > get_qr_inferred_tracts
ruleorder: evaluate_sstar > evaluate_qr


rule sstar_score:
    input:
        vcf=rules.simulate_test_data.output.vcf,
        ref_list=rules.simulate_test_data.output.ref_list,
        tgt_list=rules.simulate_test_data.output.tgt_list,
    output:
        scores="results/simulation/test/rep_{test_rep}/sstar.phased.rep_{test_rep}.scores.tsv",
    params:
        win_len=50000,
        win_step=50000,
    resources:
        mem_gb=16, cpus=4,
    conda:
        "../envs/sstar.yaml",
    shell:
        """
        sstar score \
          --vcf {input.vcf} \
          --ref {input.ref_list} \
          --tgt {input.tgt_list} \
          --output {output.scores} \
          --thread {resources.cpus} \
          --win-len {params.win_len} \
          --win-step {params.win_step} \
          --phased \
        """


rule sstar_quantile:
    input:
        model="config/ArchIE_3D19_wo_introgression.yaml",
    output:
        quantile="results/sstar/rep_{test_rep}/quantile.summary.txt",
    params:
        ms_dir="resources/msdir",
        output_dir="results/sstar/rep_{test_rep}",
    resources:
        time=360, mem_gb=128, cpus=32,
    conda:
        "../envs/sstar.yaml",
    shell:
        """
        sstar quantile \
          --model {input.model} \
          --ms-dir {params.ms_dir} \
          --N0 1000 \
          --nsamp 200 \
          --nreps 10000 \
          --ref-index 3 \
          --ref-size 100 \
          --tgt-index 4 \
          --tgt-size 100 \
          --mut-rate 1.2e-8 \
          --rec-rate 1.0e-8 \
          --seq-len 50000 \
          --snp-num-range 50 350 5 \
          --output-dir {params.output_dir} \
          --thread {resources.cpus} \
        """


rule sstar_threshold:
    input:
        scores=rules.sstar_score.output.scores,
        quantile=rules.sstar_quantile.output.quantile,
    output:
        preds="results/sstar/rep_{test_rep}/sstar.phased.q_{quantile}.rep_{test_rep}.preds.tsv",
    conda:
        "../envs/sstar.yaml",
    shell:
        """
        sstar threshold \
          --score {input.scores} \
          --sim-data {input.quantile} \
          --recomb-rate 1.0e-8 \
          --quantile {wildcards.quantile} \
          --output {output.preds} \
          --k 8 \
          --phased
        """


rule get_sstar_inferred_tracts:
    input:
        preds=rules.sstar_threshold.output.preds,
    output:
        bed="results/sstar/rep_{test_rep}/sstar.phased.q_{quantile}.rep_{test_rep}.inferred.tracts.bed",
    shell:
        r"""
        awk 'BEGIN{{FS=OFS="\t"}} NR==1{{next}} $5>$6 {{print $1,$2,$3,$4}}' {input.preds} | sed 's/hap//' > {output.bed}
        """


rule evaluate_sstar:
    input:
        true_tracts=rules.simulate_test_data.output.bed,
        inferred_tracts=rules.get_sstar_inferred_tracts.output.bed,
    output:
        tsv="results/sstar/rep_{test_rep}/sstar.phased.q_{quantile}.rep_{test_rep}.perf.tsv",
    params:
        length_bp=200_000_000,
        cutoff="{quantile}",
    script:
        "../scripts/segment_based_evaluation.py"
