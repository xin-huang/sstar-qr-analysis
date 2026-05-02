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


rule collect_qr_performance_across_cutoffs:
    input:
        perf=expand(
            rules.evaluate_qr.output.tsv,
            quantile=cutoffs,
            allow_missing=True,
        ),
    output:
        perf=temp("results/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/{qr_model}/{phase_state}/rep_{test_rep}/{qr_model}.{feature_set}.pred.perf.tsv"),
    shell:
        r"""
        cat {input.perf} | grep -v Cutoff | awk -v method="{wildcards.qr_model}_{wildcards.feature_set}" -v nref={wildcards.n_ref} -v ntgt={wildcards.n_tgt} '{{print method"\t"nref"\t"ntgt"\t"$0}}' > {output.perf}
        sed -i '1iMethod\tN_ref\tN_tgt\tCutoff\tPrecision\tRecall\tL_TT_sample\tL_IT_sample\tL_TP_sample\tL_FP_sample\tL_TN_sample\tL_FN_sample' {output.perf}
        """


rule collect_sstar_performance_across_cutoffs:
    input:
        perf=expand(
            "results/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/sstar/{phase_state}/rep_{test_rep}/sstar.{phase_state}.q_{quantile}.rep_{test_rep}.perf.tsv",
            quantile=cutoffs,
            allow_missing=True,
        ),
    output:
        perf=temp("results/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/sstar/{phase_state}/rep_{test_rep}/sstar.pred.perf.tsv"),
    shell:
        r"""
        cat {input.perf} | grep -v Cutoff | awk -v nref={wildcards.n_ref} -v ntgt={wildcards.n_tgt} '{{print "sstar\t"nref"\t"ntgt"\t"$0}}' > {output.perf}
        sed -i '1iMethod\tN_ref\tN_tgt\tCutoff\tPrecision\tRecall\tL_TT_sample\tL_IT_sample\tL_TP_sample\tL_FP_sample\tL_TN_sample\tL_FN_sample' {output.perf}
        """
