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


def get_sstar2_2src_config_params(wildcards):
    pop = get_pop_config(wildcards)

    return {
        "ref_id": pop.ref,
        "tgt_id": pop.tgt,
        "mut_rate": pop.mut_rate,
        "rec_rate": pop.rec_rate,
        "win_len": pop.win_len,
        "nprocess": 32,
        "phase_state": wildcards.phase_state == "phased",
        "seed": int(seed_lists["test"][int(wildcards.test_rep)]),
        "ref_ind_file": (
            f"results/{wildcards.demog_model}/nref_{wildcards.n_ref}/"
            f"ntgt_{wildcards.n_tgt}/nsrc_{wildcards.n_src}/"
            f"simulation/test/rep_{wildcards.test_rep}/"
            f"simulation.rep_{wildcards.test_rep}.ref.list"
        ),
        "tgt_ind_file": (
            f"results/{wildcards.demog_model}/nref_{wildcards.n_ref}/"
            f"ntgt_{wildcards.n_tgt}/nsrc_{wildcards.n_src}/"
            f"simulation/test/rep_{wildcards.test_rep}/"
            f"simulation.rep_{wildcards.test_rep}.tgt.list"
        ),
    }


rule sstar2_2src_all:
    input:
        expand(
            "results/2src/sstar2/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "{phase_state}/rep_{test_rep}/q_{quantile}/"
            "sstar2.q_{quantile}.rep_{test_rep}.{source_name}.inferred.tracts.bed",
            demog_model=TWO_SOURCE_MODELS,
            n_ref=N_REFS,
            n_tgt=N_TGTS,
            n_src=N_SRCS,
            phase_state=PHASE_STATES,
            test_rep=range(TEST_REP),
            quantile=cutoffs,
            source_name=["src1", "src2"],
        )


rule render_sstar2_2src_config_template:
    input:
        template="config/sstar2.config.template.yaml",
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
        config=(
            "results/2src/sstar2/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "{phase_state}/rep_{test_rep}/q_{quantile}/"
            "sstar2.q_{quantile}.rep_{test_rep}.config.yaml"
        ),
    wildcard_constraints:
        demog_model=TWO_SOURCE_MODELS_REGEX,
    params:
        sstar2=get_sstar2_2src_config_params,
    template_engine:
        "yte"


rule run_sstar2_2src_train:
    input:
        demes="config/demog_models/{demog_model}_wo_introgression.yaml",
        config=(
            "results/2src/sstar2/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "{phase_state}/rep_{test_rep}/q_{quantile}/"
            "sstar2.q_{quantile}.rep_{test_rep}.config.yaml"
        ),
    output:
        model=(
            "results/2src/sstar2/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "{phase_state}/rep_{test_rep}/q_{quantile}/"
            "sstar2.q_{quantile}.rep_{test_rep}.onnx"
        ),
    wildcard_constraints:
        demog_model=TWO_SOURCE_MODELS_REGEX,
    resources:
        mem_mb=64000,
        cpus=32,
        time=360,
    conda:
        "../envs/sstar2.yaml",
    shell:
        """
        sstar2 train \
          --demes {input.demes} \
          --config {input.config} \
          --output {output.model}
        """


rule run_sstar2_2src_infer:
    input:
        model=(
            "results/2src/sstar2/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "{phase_state}/rep_{test_rep}/q_{quantile}/"
            "sstar2.q_{quantile}.rep_{test_rep}.onnx"
        ),
        config=(
            "results/2src/sstar2/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "{phase_state}/rep_{test_rep}/q_{quantile}/"
            "sstar2.q_{quantile}.rep_{test_rep}.config.yaml"
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
        feat=(
            "results/2src/sstar2/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "{phase_state}/rep_{test_rep}/q_{quantile}/"
            "sstar2.q_{quantile}.rep_{test_rep}.processed.features.tsv"
        ),
        pred=(
            "results/2src/sstar2/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "{phase_state}/rep_{test_rep}/q_{quantile}/"
            "sstar2.q_{quantile}.rep_{test_rep}.pred.tsv"
        ),
        tracts=(
            "results/2src/sstar2/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "{phase_state}/rep_{test_rep}/q_{quantile}/"
            "sstar2.q_{quantile}.rep_{test_rep}.candidate.tracts.bed"
        ),
    wildcard_constraints:
        demog_model=TWO_SOURCE_MODELS_REGEX,
    resources:
        mem_mb=64000,
        cpus=32,
    conda:
        "../envs/sstar2.yaml",
    shell:
        """
        sstar2 infer \
          --model {input.model} \
          --config {input.config} \
          --feat-file {output.feat} \
          --pred-file {output.pred} \
          --tract-file {output.tracts}
        """


rule sstar2_2src_match_src1:
    input:
        vcf=(
            "results/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "simulation/test/rep_{test_rep}/"
            "simulation.rep_{test_rep}.biallelic.snps.vcf.gz"
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
        tracts=(
            "results/2src/sstar2/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "{phase_state}/rep_{test_rep}/q_{quantile}/"
            "sstar2.q_{quantile}.rep_{test_rep}.candidate.tracts.bed"
        ),
    output:
        matchrate=(
            "results/2src/sstar2/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "{phase_state}/rep_{test_rep}/q_{quantile}/"
            "sstar2.q_{quantile}.rep_{test_rep}.src1.match.rate.bed"
        ),
    wildcard_constraints:
        demog_model=TWO_SOURCE_MODELS_REGEX,
    resources:
        mem_mb=16000,
        time=1440,
    conda:
        "../envs/sstar2.yaml",
    shell:
        """
        sstar2 match \
          --vcf {input.vcf} \
          --tgt {input.tgt_list} \
          --src {input.src_list} \
          --tract-file {input.tracts} \
          --output {output.matchrate}
        """


rule sstar2_2src_match_src2:
    input:
        vcf=(
            "results/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "simulation/test/rep_{test_rep}/"
            "simulation.rep_{test_rep}.biallelic.snps.vcf.gz"
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
        tracts=(
            "results/2src/sstar2/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "{phase_state}/rep_{test_rep}/q_{quantile}/"
            "sstar2.q_{quantile}.rep_{test_rep}.candidate.tracts.bed"
        ),
    output:
        matchrate=(
            "results/2src/sstar2/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "{phase_state}/rep_{test_rep}/q_{quantile}/"
            "sstar2.q_{quantile}.rep_{test_rep}.src2.match.rate.bed"
        ),
    wildcard_constraints:
        demog_model=TWO_SOURCE_MODELS_REGEX,
    resources:
        mem_mb=16000,
        time=1440,
    conda:
        "../envs/sstar2.yaml",
    shell:
        """
        sstar2 match \
          --vcf {input.vcf} \
          --tgt {input.tgt_list} \
          --src {input.src_list} \
          --tract-file {input.tracts} \
          --output {output.matchrate}
        """


rule sstar2_2src_assign:
    input:
        src1_matchrate=(
            "results/2src/sstar2/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "{phase_state}/rep_{test_rep}/q_{quantile}/"
            "sstar2.q_{quantile}.rep_{test_rep}.src1.match.rate.bed"
        ),
        src2_matchrate=(
            "results/2src/sstar2/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "{phase_state}/rep_{test_rep}/q_{quantile}/"
            "sstar2.q_{quantile}.rep_{test_rep}.src2.match.rate.bed"
        ),
    output:
        src1=(
            "results/2src/sstar2/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "{phase_state}/rep_{test_rep}/q_{quantile}/"
            "sstar2.q_{quantile}.rep_{test_rep}.src1.inferred.tracts.bed"
        ),
        src2=(
            "results/2src/sstar2/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "{phase_state}/rep_{test_rep}/q_{quantile}/"
            "sstar2.q_{quantile}.rep_{test_rep}.src2.inferred.tracts.bed"
        ),
    wildcard_constraints:
        demog_model=TWO_SOURCE_MODELS_REGEX,
    params:
        output_prefix=(
            "results/2src/sstar2/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "{phase_state}/rep_{test_rep}/q_{quantile}/"
            "sstar2.q_{quantile}.rep_{test_rep}"
        ),
    conda:
        "../envs/sstar2.yaml",
    shell:
        """
        sstar2 assign \
          --match-rate {input.src1_matchrate} {input.src2_matchrate} \
          --source-name src1 src2 \
          --output-prefix {params.output_prefix}
        """
