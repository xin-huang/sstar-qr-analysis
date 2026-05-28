# Copyright 2025 Xin Huang
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


import numpy as np
import pandas as pd
import pyranges as pr


def calc_pr(
    nintro_tracts: int,
    ninferred_tracts: int,
    ntrue_positives: int,
) -> tuple[float, float]:

    if ninferred_tracts == 0:
        precision = np.nan
    else:
        precision = ntrue_positives / float(ninferred_tracts)

    if nintro_tracts == 0:
        recall = np.nan
    else:
        recall = ntrue_positives / float(nintro_tracts)

    return precision, recall


def tract_length(tracts) -> int:
    if tracts is None or len(tracts) == 0:
        return 0

    return int((tracts.End - tracts.Start).sum())


def read_tracts(path: str):
    try:
        tracts = pd.read_csv(
            path,
            sep="\t",
            header=None,
            names=["Chromosome", "Start", "End", "Sample"],
            dtype={"Chromosome": str, "Sample": str},
        )
    except pd.errors.EmptyDataError:
        return None, np.array([])

    if tracts.empty:
        return None, np.array([])

    tracts["Start"] = tracts["Start"].astype(int)
    tracts["End"] = tracts["End"].astype(int)

    samples = tracts["Sample"].unique()
    tracts = pr.PyRanges(tracts).merge(by="Sample")

    return tracts, samples


def evaluate(
    true_tract_file: str,
    inferred_tract_file: str,
    seq_len: int,
    cutoff: float,
    output: str,
) -> None:

    intro_tracts, intro_tracts_samples = read_tracts(true_tract_file)
    inferred_tracts, inferred_tracts_samples = read_tracts(inferred_tract_file)

    all_samples = np.union1d(intro_tracts_samples, inferred_tracts_samples)
    sample_size = len(all_samples)

    sum_nintro_tracts = 0
    sum_ninferred_tracts = 0
    sum_ntrue_positives = 0
    sum_nfalse_positives = 0
    sum_ntrue_negatives = 0
    sum_nfalse_negatives = 0

    for s in all_samples:
        if intro_tracts is None:
            ind_intro_tracts = None
        else:
            ind_intro_tracts = intro_tracts[intro_tracts.Sample == s]

        if inferred_tracts is None:
            ind_inferred_tracts = None
        else:
            ind_inferred_tracts = inferred_tracts[inferred_tracts.Sample == s]

        nintro_tracts = tract_length(ind_intro_tracts)
        ninferred_tracts = tract_length(ind_inferred_tracts)

        if (
            ind_intro_tracts is None
            or ind_inferred_tracts is None
            or len(ind_intro_tracts) == 0
            or len(ind_inferred_tracts) == 0
        ):
            ntrue_positives = 0
        else:
            ind_overlaps = ind_intro_tracts.intersect(ind_inferred_tracts)
            ntrue_positives = tract_length(ind_overlaps)

        nfalse_positives = ninferred_tracts - ntrue_positives
        nfalse_negatives = nintro_tracts - ntrue_positives

        # true nonintro = FP + TN
        # TN = true nonintro - FP
        #    = seq_len - nintro_tracts - nfalse_positives
        #    = seq_len - nintro_tracts - ninferred_tracts + ntrue_positives
        ntrue_negatives = (
            seq_len
            - nintro_tracts
            - ninferred_tracts
            + ntrue_positives
        )

        sum_nintro_tracts += nintro_tracts
        sum_ninferred_tracts += ninferred_tracts
        sum_ntrue_positives += ntrue_positives
        sum_nfalse_positives += nfalse_positives
        sum_nfalse_negatives += nfalse_negatives
        sum_ntrue_negatives += ntrue_negatives

    total_precision, total_recall = calc_pr(
        sum_nintro_tracts,
        sum_ninferred_tracts,
        sum_ntrue_positives,
    )

    if sample_size == 0:
        per_sample_metrics = [np.nan] * 6
    else:
        per_sample_metrics = [
            sum_nintro_tracts / sample_size,
            sum_ninferred_tracts / sample_size,
            sum_ntrue_positives / sample_size,
            sum_nfalse_positives / sample_size,
            sum_ntrue_negatives / sample_size,
            sum_nfalse_negatives / sample_size,
        ]

    res = pd.DataFrame(
        [[
            cutoff,
            total_precision,
            total_recall,
            *per_sample_metrics,
        ]],
        columns=[
            "Cutoff",
            "Precision",
            "Recall",
            "L_TT_sample",
            "L_IT_sample",
            "L_TP_sample",
            "L_FP_sample",
            "L_TN_sample",
            "L_FN_sample",
        ],
    )

    res.fillna("NaN").to_csv(output, sep="\t", index=False)


evaluate(
    true_tract_file=snakemake.input.true_tracts,
    inferred_tract_file=snakemake.input.inferred_tracts,
    seq_len=int(snakemake.params.length_bp),
    cutoff=float(snakemake.params.cutoff),
    output=snakemake.output.tsv,
)
