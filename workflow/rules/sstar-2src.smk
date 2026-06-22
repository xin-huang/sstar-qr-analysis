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

rule sstar_2src_all:
    input:
        expand(
            "results/2src/sstar/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "{phase_state}/rep_{test_rep}/"
            "sstar.{phase_state}.q_{quantile}.rep_{test_rep}.{source_name}.inferred.tracts.bed",
            demog_model=TWO_SOURCE_MODELS,
            n_ref=N_REFS,
            n_tgt=N_TGTS,
            n_src=N_SRCS,
            phase_state=["unphased"],
            test_rep=range(TEST_REP),
            quantile=cutoffs,
            source_name=["src1", "src2"],
        )


rule sstar_2src_score:
    input:
        vcf=(
            "results/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "simulation/test/rep_{test_rep}/"
            "simulation.rep_{test_rep}.biallelic.snps.vcf.gz"
        ),
        ref_list=(
            "results/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "simulation/test/rep_{test_rep}/"
            "simulation.rep_{test_rep}.ref.list"
        ),
        tgt_list=(
            "results/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "simulation/test/rep_{test_rep}/"
            "simulation.rep_{test_rep}.tgt.list"
        ),
    output:
        scores=(
            "results/2src/sstar/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "{phase_state}/rep_{test_rep}/"
            "sstar.{phase_state}.rep_{test_rep}.scores.tsv"
        ),
    wildcard_constraints:
        demog_model=TWO_SOURCE_MODELS_REGEX,
        phase_state="unphased",
    params:
        pop_config=get_pop_config,
        win_step=10000,
    resources:
        time=360,
        mem_mb=16000,
        cpus=16,
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
          --win-step {params.win_step}
        """


rule sstar_2src_quantile:
    input:
        model="config/demog_models/{demog_model}_wo_introgression.yaml",
    output:
        quantile=(
            "results/2src/sstar/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "{phase_state}/rep_{test_rep}/"
            "quantile.summary.txt"
        ),
    wildcard_constraints:
        demog_model=TWO_SOURCE_MODELS_REGEX,
        phase_state="unphased",
    params:
        ms_dir="resources/msdir",
        output_dir=(
            "results/2src/sstar/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "{phase_state}/rep_{test_rep}"
        ),
        pop_config=get_pop_config,
        sample_size=get_sample_size,
        quantile_start=float(cutoffs.min()),
        quantile_step=0.001,
    resources:
        time=360,
        mem_mb=64000,
        cpus=32,
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
          --quantile-start {params.quantile_start} \
          --quantile-step {params.quantile_step} \
          --output-dir {params.output_dir} \
          --thread {resources.cpus}
        """


rule sstar_2src_threshold:
    input:
        scores=(
            "results/2src/sstar/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "{phase_state}/rep_{test_rep}/"
            "sstar.{phase_state}.rep_{test_rep}.scores.tsv"
        ),
        quantile=(
            "results/2src/sstar/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "{phase_state}/rep_{test_rep}/"
            "quantile.summary.txt"
        ),
    output:
        threshold=(
            "results/2src/sstar/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "{phase_state}/rep_{test_rep}/"
            "sstar.{phase_state}.q_{quantile}.rep_{test_rep}.threshold.tsv"
        ),
    wildcard_constraints:
        demog_model=TWO_SOURCE_MODELS_REGEX,
        phase_state="unphased",
    conda:
        "../envs/sstar.yaml",
    shell:
        """
        sstar threshold \
          --score {input.scores} \
          --sim-data {input.quantile} \
          --quantile {wildcards.quantile} \
          --output {output.threshold} \
          --k 8
        """


rule sstar_2src_matchrate_src1:
    input:
        vcf=(
            "results/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "simulation/test/rep_{test_rep}/"
            "simulation.rep_{test_rep}.biallelic.snps.vcf.gz"
        ),
        ref_list=(
            "results/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "simulation/test/rep_{test_rep}/"
            "simulation.rep_{test_rep}.ref.list"
        ),
        tgt_list=(
            "results/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "simulation/test/rep_{test_rep}/"
            "simulation.rep_{test_rep}.tgt.list"
        ),
        src_list=(
            "results/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "simulation/test/rep_{test_rep}/"
            "simulation.rep_{test_rep}.src.list"
        ),
        scores=(
            "results/2src/sstar/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "{phase_state}/rep_{test_rep}/"
            "sstar.{phase_state}.rep_{test_rep}.scores.tsv"
        ),
    output:
        matchrate=(
            "results/2src/sstar/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "{phase_state}/rep_{test_rep}/"
            "sstar.{phase_state}.rep_{test_rep}.src1.matchrate.tsv"
        ),
    wildcard_constraints:
        demog_model=TWO_SOURCE_MODELS_REGEX,
        phase_state="unphased",
    resources:
        time=360,
        mem_mb=16000,
        cpus=16,
    conda:
        "../envs/sstar.yaml",
    shell:
        """
        sstar matchrate \
          --vcf {input.vcf} \
          --ref {input.ref_list} \
          --tgt {input.tgt_list} \
          --src {input.src_list} \
          --score {input.scores} \
          --output {output.matchrate} \
          --thread {resources.cpus}
        """


rule sstar_2src_matchrate_src2:
    input:
        vcf=(
            "results/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "simulation/test/rep_{test_rep}/"
            "simulation.rep_{test_rep}.biallelic.snps.vcf.gz"
        ),
        ref_list=(
            "results/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "simulation/test/rep_{test_rep}/"
            "simulation.rep_{test_rep}.ref.list"
        ),
        tgt_list=(
            "results/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "simulation/test/rep_{test_rep}/"
            "simulation.rep_{test_rep}.tgt.list"
        ),
        src_list=(
            "results/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "simulation/test/rep_{test_rep}/"
            "simulation.rep_{test_rep}.src2.list"
        ),
        scores=(
            "results/2src/sstar/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "{phase_state}/rep_{test_rep}/"
            "sstar.{phase_state}.rep_{test_rep}.scores.tsv"
        ),
    output:
        matchrate=(
            "results/2src/sstar/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "{phase_state}/rep_{test_rep}/"
            "sstar.{phase_state}.rep_{test_rep}.src2.matchrate.tsv"
        ),
    wildcard_constraints:
        demog_model=TWO_SOURCE_MODELS_REGEX,
        phase_state="unphased",
    resources:
        time=360,
        mem_mb=16000,
        cpus=16,
    conda:
        "../envs/sstar.yaml",
    shell:
        """
        sstar matchrate \
          --vcf {input.vcf} \
          --ref {input.ref_list} \
          --tgt {input.tgt_list} \
          --src {input.src_list} \
          --score {input.scores} \
          --output {output.matchrate} \
          --thread {resources.cpus}
        """


rule sstar_2src_tract:
    input:
        threshold=(
            "results/2src/sstar/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "{phase_state}/rep_{test_rep}/"
            "sstar.{phase_state}.q_{quantile}.rep_{test_rep}.threshold.tsv"
        ),
        src1_matchrate=(
            "results/2src/sstar/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "{phase_state}/rep_{test_rep}/"
            "sstar.{phase_state}.rep_{test_rep}.src1.matchrate.tsv"
        ),
        src2_matchrate=(
            "results/2src/sstar/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "{phase_state}/rep_{test_rep}/"
            "sstar.{phase_state}.rep_{test_rep}.src2.matchrate.tsv"
        ),
    output:
        src1=(
            "results/2src/sstar/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "{phase_state}/rep_{test_rep}/"
            "sstar.{phase_state}.q_{quantile}.rep_{test_rep}.src1.inferred.tracts.bed"
        ),
        src2=(
            "results/2src/sstar/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "{phase_state}/rep_{test_rep}/"
            "sstar.{phase_state}.q_{quantile}.rep_{test_rep}.src2.inferred.tracts.bed"
        ),
    wildcard_constraints:
        demog_model=TWO_SOURCE_MODELS_REGEX,
        phase_state="unphased",
    params:
        output_prefix=(
            "results/2src/sstar/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "{phase_state}/rep_{test_rep}/"
            "sstar.{phase_state}.q_{quantile}.rep_{test_rep}"
        ),
    conda:
        "../envs/sstar.yaml",
    shell:
        """
        sstar tract \
          --threshold {input.threshold} \
          --output-prefix {params.output_prefix} \
          --match-rate {input.src1_matchrate} {input.src2_matchrate}

        mv {params.output_prefix}.src1.bed {output.src1}
        mv {params.output_prefix}.src2.bed {output.src2}
        """
