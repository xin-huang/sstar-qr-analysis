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


import numpy as np
import pandas as pd
import pyranges as pr


def calc_pr(
    ntrue_tracts: int, ninferred_tracts: int, ntrue_positives: int
) -> tuple[float, float]:
    """
    Calculate segment-based precision and recall.

    Parameters
    ----------
    ntrue_tracts : int
        Length of true introgressed fragments.
    ninferred_tracts : int
        Length of inferred introgressed fragments.
    ntrue_positives : int
        Length of fragments belonging to true positives.

    Returns
    -------
    precision : float
        Estimated precision.
    recall : float
        Estimated recall.
    """
    if float(ninferred_tracts) == 0:
        precision = np.nan
    else:
        precision = ntrue_positives / float(ninferred_tracts)
    if float(ntrue_tracts) == 0:
        recall = np.nan
    else:
        recall = ntrue_positives / float(ntrue_tracts)

    return precision, recall


def evaluate(
    true_tract_file: str,
    inferred_tract_file: str,
    seq_len: int,
    cutoff: float,
    output: str,
) -> None:
    """
    Evaluate model performance using segment-based precision and recall.

    Parameters
    ----------
    true_tract_file : str
        Path to the file containing true introgressed fragments.
    inferred_tract_file : str
        Path to the file containing inferred introgressed fragments.
    seq_len : int
        Total length of the sequence.
    cutoff : float
        Probability threshold used to classify a fragment as introgressed.
    output : str
        Path to the output file where stores model performance metrics.
    """
    try:
        true_tracts = pd.read_csv(
            true_tract_file,
            sep="\t",
            header=None,
            names=["Chromosome", "Start", "End", "Sample"],
            dtype={"Chromosome": str, "Start": int, "End": int, "Sample": str},
        )
    except pd.errors.EmptyDataError:
        true_tracts_samples = []
    else:
        true_tracts_samples = true_tracts["Sample"].unique()
        true_tracts = pr.PyRanges(true_tracts).merge(by="Sample")

    try:
        inferred_tracts = pd.read_csv(
            inferred_tract_file,
            sep="\t",
            header=None,
            names=["Chromosome", "Start", "End", "Sample"],
            dtype={"Chromosome": str, "Start": int, "End": int, "Sample": str},
        )
    except pd.errors.EmptyDataError:
        inferred_tracts_samples = []
    else:
        inferred_tracts_samples = inferred_tracts["Sample"].unique()
        inferred_tracts = pr.PyRanges(inferred_tracts).merge(by="Sample")

    res = pd.DataFrame(
        columns=[
            "Cutoff",
            "Precision",
            "Recall",
        ]
    )

    precisions = []
    recalls = []

    overlap = np.intersect1d(true_tracts_samples, inferred_tracts_samples)
    true_tracts_only = np.setdiff1d(true_tracts_samples, inferred_tracts_samples)
    inferred_tracts_only = np.setdiff1d(inferred_tracts_samples, true_tracts_samples)

    # Samples exist in both true_tracts and inferred_tracts
    for s in overlap:
        ind_true_tracts = true_tracts[true_tracts.Sample == s]
        ind_inferred_tracts = inferred_tracts[inferred_tracts.Sample == s]
        ind_overlaps = ind_true_tracts.intersect(ind_inferred_tracts)

        ntrue_tracts = (ind_true_tracts.End - ind_true_tracts.Start).sum()
        ninferred_tracts = (ind_inferred_tracts.End - ind_inferred_tracts.Start).sum()
        ntrue_positives = (
            (ind_overlaps.End - ind_overlaps.Start).sum() if len(ind_overlaps) > 0 else 0
        )

        precision, recall = calc_pr(ntrue_tracts, ninferred_tracts, ntrue_positives)
        precisions.append(precision)
        recalls.append(recall)

    # Samples only exist in true_tracts
    for s in true_tracts_only:
        ind_true_tracts = true_tracts[true_tracts.Sample == s]

        ninferred_tracts = 0
        ntrue_positives = 0
        ntrue_tracts = (ind_true_tracts.End - ind_true_tracts.Start).sum()

        precision, recall = calc_pr(ntrue_tracts, ninferred_tracts, ntrue_positives)
        precisions.append(precision)
        recalls.append(recall)

    # Samples only exist in inferred_tracts
    for s in inferred_tracts_only:
        ind_inferred_tracts = inferred_tracts[inferred_tracts.Sample == s]

        ntrue_tracts = 0
        ntrue_positives = 0
        ninferred_tracts = (ind_inferred_tracts.End - ind_inferred_tracts.Start).sum()

        precision, recall = calc_pr(ntrue_tracts, ninferred_tracts, ntrue_positives)
        precisions.append(precision)
        recalls.append(recall)

    avg_precision = np.nanmean(precisions)
    avg_recall = np.nanmean(recalls)

    res.loc[len(res.index)] = [
        cutoff,
        avg_precision,
        avg_recall,
    ]

    res.fillna("NaN").to_csv(output, sep="\t", index=False)


evaluate(
    true_tract_file=snakemake.input.true_tracts,
    inferred_tract_file=snakemake.input.inferred_tracts,
    seq_len=int(snakemake.params.length_bp),
    cutoff=float(snakemake.params.cutoff),
    output=snakemake.output.tsv,
)
