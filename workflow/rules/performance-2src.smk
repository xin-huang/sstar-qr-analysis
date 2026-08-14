TWO_SOURCE_NAMES = ["src1", "src2"]


rule collect_sstar_2src_performance_across_cutoffs:
    input:
        perf=expand(
            "results/2src/{version}/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "{phase_state}/rep_{test_rep}/"
            "{version}.{phase_state}.q_{quantile}.rep_{test_rep}.{source_name}.perf.tsv",
            quantile=cutoffs,
            allow_missing=True,
        ),
    output:
        perf=temp(
            "results/2src/{version}/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "{phase_state}/rep_{test_rep}/{version}.{source_name}.pred.perf.tsv"
        ),
    wildcard_constraints:
        demog_model=TWO_SOURCE_MODELS_REGEX,
        source_name="src1|src2",
        version="sstar|sstar2",
    shell:
        r"""
        cat {input.perf} | grep -v Cutoff | awk -v rep={wildcards.test_rep} -v version={wildcards.version} -v nref={wildcards.n_ref} -v ntgt={wildcards.n_tgt} '{{print rep"\t"version"\t"nref"\t"ntgt"\t"$0}}' > {output.perf}
        sed -i '1iReplicate\tMethod\tN_ref\tN_tgt\tCutoff\tPrecision\tRecall' {output.perf}
        """


rule collect_2src_performance_across_replicates:
    input:
        sstar=expand(
            "results/2src/{version}/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "{phase_state}/rep_{test_rep}/{version}.{source_name}.pred.perf.tsv",
            version=["sstar", "sstar2"],
            test_rep=range(TEST_REP),
            allow_missing=True,
        ),
    output:
        perf=(
            "results/2src/performance/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "{phase_state}/{source_name}/combined.pred.perf.tsv"
        ),
    wildcard_constraints:
        demog_model=TWO_SOURCE_MODELS_REGEX,
        source_name="src1|src2",
    shell:
        r"""
        head -n 1 {input.sstar[0]} > {output.perf}
        cat {input.sstar} | grep -v Cutoff >> {output.perf}
        """


rule plot_2src_pr_curve:
    input:
        perf=expand(
            "results/2src/performance/{demog_model}/nref_{n_ref}/ntgt_{n_tgt}/nsrc_{n_src}/"
            "{phase_state}/{source_name}/combined.pred.perf.tsv",
            demog_model=TWO_SOURCE_MODELS,
            phase_state=PHASE_STATES,
            source_name=TWO_SOURCE_NAMES,
            allow_missing=True,
        ),
    output:
        plot=(
            "results/2src/plots/"
            "pred.nref_{n_ref}.ntgt_{n_tgt}.nsrc_{n_src}.pr.curve.combined.png"
        ),
        summary_tsv=(
            "results/2src/plots/"
            "pred.nref_{n_ref}.ntgt_{n_tgt}.nsrc_{n_src}.pr.curve.combined.summary.tsv"
        ),
    script:
        "../scripts/plot_pr_curve_2src.py"
