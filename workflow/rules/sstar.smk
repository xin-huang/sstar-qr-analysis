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


rule sstar_score:
    input:
        vcf="results/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/simulation/test/rep_{test_rep}/simulation.rep_{test_rep}.biallelic.snps.vcf.gz",
        ref_list="results/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/simulation/test/rep_{test_rep}/simulation.rep_{test_rep}.ref.list",
        tgt_list="results/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/simulation/test/rep_{test_rep}/simulation.rep_{test_rep}.tgt.list",
    output:
        scores="results/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/simulation/test/rep_{test_rep}/sstar.{phase_state}.rep_{test_rep}.scores.tsv",
    params:
        pop_config=get_pop_config,
        win_step=10000,
        phased_flag=lambda wildcards: "--phased" if wildcards.phase_state == "phased" else "",
    resources:
        time=360, mem_mb=16000, cpus=16,
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
          --win-len {params.pop_config.win_len} \
          --win-step {params.win_step} \
          {params.phased_flag} \
        """


rule sstar_quantile:
    input:
        model="config/demog_models/{demog_model}_wo_introgression.yaml",
    output:
        quantile="results/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/sstar/{phase_state}/rep_{test_rep}/quantile.summary.txt",
    params:
        ms_dir="resources/msdir",
        output_dir="results/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/sstar/{phase_state}/rep_{test_rep}",
        phased_flag=lambda wildcards: "--phased" if wildcards.phase_state == "phased" else "",
        pop_config=get_pop_config,
        sample_size=get_sample_size,
    resources:
        time=360, mem_mb=64000, cpus=32,
    conda:
        "../envs/sstar.yaml",
    shell:
        """
        sstar quantile \
          --model {input.model} \
          --ms-dir {params.ms_dir} \
          --N0 1000 \
          --nsamp {params.sample_size[total]} \
          --nreps 10000 \
          --ref-pop {params.pop_config.ref} \
          --ref-size {params.sample_size[ref]} \
          --tgt-pop {params.pop_config.tgt} \
          --tgt-size {params.sample_size[tgt]} \
          --mut-rate {params.pop_config.mut_rate} \
          --rec-rate {params.pop_config.rec_rate} \
          --seq-len {params.pop_config.win_len} \
          --snp-num-range 50 350 5 \
          --output-dir {params.output_dir} \
          --quantile-step 0.00001 \
          --thread {resources.cpus} \
          {params.phased_flag}
        """


rule sstar_threshold:
    input:
        scores="results/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/simulation/test/rep_{test_rep}/sstar.{phase_state}.rep_{test_rep}.scores.tsv",
        quantile="results/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/sstar/{phase_state}/rep_{test_rep}/quantile.summary.txt",
    output:
        pred=temp("results/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/sstar/{phase_state}/rep_{test_rep}/sstar.{phase_state}.q_{quantile}.rep_{test_rep}.pred.tsv"),
    params:
        phased_flag=lambda wildcards: "--phased" if wildcards.phase_state == "phased" else "",
    resources:
        mem_mb=16000,
    conda:
        "../envs/sstar.yaml",
    shell:
        """
        sstar threshold \
          --score {input.scores} \
          --sim-data {input.quantile} \
          --quantile {wildcards.quantile} \
          --output {output.pred} \
          --k 8 \
          {params.phased_flag}
        """


rule get_sstar_inferred_tracts:
    input:
        pred="results/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/sstar/{phase_state}/rep_{test_rep}/sstar.{phase_state}.q_{quantile}.rep_{test_rep}.pred.tsv",
    output:
        bed="results/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/sstar/{phase_state}/rep_{test_rep}/sstar.{phase_state}.q_{quantile}.rep_{test_rep}.inferred.tracts.bed",
    localrule: True,
    shell:
        r"""
        awk 'BEGIN{{FS=OFS="\t"}} NR==1{{next}} $5>$6 {{print $1,$2,$3,$4}}' {input.pred} | awk 'BEGIN{{FS=OFS="\t"}} {{if ("{wildcards.phase_state}"=="phased") gsub(/hap/, "", $4); print}}' > {output.bed}
        """


rule evaluate_sstar:
    input:
        true_tracts="results/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/simulation/test/rep_{test_rep}/simulation.rep_{test_rep}.true.tracts.{phase_state}.bed",
        inferred_tracts="results/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/sstar/{phase_state}/rep_{test_rep}/sstar.{phase_state}.q_{quantile}.rep_{test_rep}.inferred.tracts.bed",
    output:
        tsv="results/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/sstar/{phase_state}/rep_{test_rep}/sstar.{phase_state}.q_{quantile}.rep_{test_rep}.perf.tsv",
    params:
        length_bp=TEST_LENGTH_BP,
        cutoff="{quantile}",
    localrule: True,
    script:
        "../scripts/segment_based_evaluation.py"
