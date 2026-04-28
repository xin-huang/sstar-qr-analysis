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


rule run_qr:
    input:
        training_data=rules.merge_training_sstar_score.output.scores,
        test_data=rules.sstar_score.output.scores,
    output:
        preds="results/{qr_model}/{phase_state}/rep_{test_rep}/{qr_model}.{feature_set}.{phase_state}.q_{quantile}.rep_{test_rep}.preds.tsv",
    resources:
        time=720, mem_gb=16,
    conda:
        "../envs/qr.yaml"
    script:
        "../scripts/qr.py"


rule get_qr_inferred_tracts:
    input:
        preds=rules.run_qr.output.preds,
    output:
        bed="results/{qr_model}/{phase_state}/rep_{test_rep}/{qr_model}.{feature_set}.{phase_state}.q_{quantile}.rep_{test_rep}.inferred.tracts.bed",
    params:
        phased=lambda wildcards: wildcards.phase_state == "phased",
    shell:
        r"""
        awk 'BEGIN{{FS=OFS="\t"}} NR==1{{next}} $5>$9 {{print $1,$2,$3,$4}}' {input.preds} | awk 'BEGIN{{FS=OFS="\t"}} {{if ("{params.phased}"=="True") gsub(/hap/, "", $4); print}}' > {output.bed}
        """


rule evaluate_qr:
    input:
        true_tracts="results/simulation/test/rep_{test_rep}/simulation.rep_{test_rep}.true.tracts.{phase_state}.bed",
        inferred_tracts=rules.get_qr_inferred_tracts.output.bed,
    output:
        tsv="results/{qr_model}/{phase_state}/rep_{test_rep}/{qr_model}.{feature_set}.{phase_state}.q_{quantile}.rep_{test_rep}.perf.tsv",
    params:
        length_bp=200_000_000,
        cutoff="{quantile}",
    script:
        "../scripts/segment_based_evaluation.py"
