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


rule render_sstar2_config_template:
    input:
        template="config/sstar2.config.template.yaml",
        vcf="results/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/simulation/test/rep_{test_rep}/simulation.rep_{test_rep}.biallelic.snps.vcf.gz",
    output:
        config="results/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/sstar2/{phase_state}/rep_{test_rep}/q_{quantile}/sstar2.q_{quantile}.rep_{test_rep}.config.yaml",
    params:
        sstar2=get_sstar2_config_params,
    localrule: True,
    template_engine:
        "yte"


rule run_sstar2_train:
    input:
        demes="config/demog_models/{demog_model}_wo_introgression.yaml",
        config="results/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/sstar2/{phase_state}/rep_{test_rep}/q_{quantile}/sstar2.q_{quantile}.rep_{test_rep}.config.yaml",
    output:
        model="results/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/sstar2/{phase_state}/rep_{test_rep}/q_{quantile}/sstar2.q_{quantile}.rep_{test_rep}.onnx",
    resources:
        mem_mb=64000, cpus=32, time=360,
    conda:
        "../envs/sstar2.yaml",
    shell:
        """
        sstar2 train --demes {input.demes} --config {input.config} --output {output.model}
        """


rule run_sstar2_infer:
    input:
        model="results/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/sstar2/{phase_state}/rep_{test_rep}/q_{quantile}/sstar2.q_{quantile}.rep_{test_rep}.onnx",
        config="results/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/sstar2/{phase_state}/rep_{test_rep}/q_{quantile}/sstar2.q_{quantile}.rep_{test_rep}.config.yaml",
    output:
        feat=temp("results/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/sstar2/{phase_state}/rep_{test_rep}/q_{quantile}/sstar2.q_{quantile}.rep_{test_rep}.processed.features.tsv"),
        pred=temp("results/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/sstar2/{phase_state}/rep_{test_rep}/q_{quantile}/sstar2.q_{quantile}.rep_{test_rep}.pred.tsv"),
        tracts="results/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/sstar2/{phase_state}/rep_{test_rep}/q_{quantile}/sstar2.q_{quantile}.rep_{test_rep}.inferred.tracts.bed",
    resources:
        mem_mb=64000, cpus=32,
    conda:
        "../envs/sstar2.yaml",
    shell:
        """
        sstar2 infer --model {input.model} --config {input.config} --feat-file {output.feat} --pred-file {output.pred} --tract-file {output.tracts}
        """


rule evaluate_sstar2:
    input:
        true_tracts="results/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/simulation/test/rep_{test_rep}/simulation.rep_{test_rep}.true.tracts.{phase_state}.bed",
        inferred_tracts="results/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/sstar2/{phase_state}/rep_{test_rep}/q_{quantile}/sstar2.q_{quantile}.rep_{test_rep}.inferred.tracts.bed",
    output:
        tsv="results/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/sstar2/{phase_state}/rep_{test_rep}/sstar2.{phase_state}.q_{quantile}.rep_{test_rep}.perf.tsv",
    params:
        length_bp=TEST_LENGTH_BP,
        cutoff="{quantile}",
    localrule: True,
    script:
        "../scripts/segment_based_evaluation.py"
