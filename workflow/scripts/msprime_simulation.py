"""
Copyright 2025 Xin Huang

GNU General Public License v3.0

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, please see

   https://www.gnu.org/licenses/gpl-3.0.en.html
"""

import demes
import msprime
import tskit
import pyranges as pr
from typing import Optional


def get_haplotype_index(
    ts: tskit.TreeSequence,
    node_id: int,
    expected_ploidy: Optional[int] = None,
) -> int:
    """
    Get 1-based haplotype index of a node within its associated individual.

    Parameters
    ----------
    ts : tskit.TreeSequence
        Tree sequence containing node and individual tables.
    node_id : int
        Node identifier whose haplotype index will be computed.
    expected_ploidy : int, optional
        Expected number of nodes associated with the node's individual. If
        provided, the function validates that ``len(individual.nodes)`` equals
        this value.

    Returns
    -------
    int
        One-based position of ``node_id`` in ``ts.individual(ind_id).nodes``.

    Raises
    ------
    ValueError
        If ``node_id`` has no associated individual, if the individual's node
        count does not match ``expected_ploidy``, or if ``node_id`` is not found
        in the associated individual's node list.
    """
    ind_id = ts.node(node_id).individual
    if ind_id == tskit.NULL:
        raise ValueError(f"Node {node_id} has no associated individual.")

    ind_nodes = list(ts.individual(ind_id).nodes)
    if expected_ploidy is not None and len(ind_nodes) != expected_ploidy:
        raise ValueError(
            f"Individual {ind_id} has {len(ind_nodes)} nodes, expected {expected_ploidy}."
        )
    try:
        return ind_nodes.index(node_id) + 1
    except ValueError as e:
        raise ValueError(
            f"Node {node_id} is not listed in individual {ind_id} nodes."
        ) from e


def simulate(
    demog: str,
    nref: int,
    ntgt: int,
    nsrc: int,
    ref_id: str,
    tgt_id: str,
    src_id: str,
    seq_len: int,
    mut_rate: float,
    rec_rate: float,
    seed: int,
    ploidy: int = 2,
) -> tskit.TreeSequence:
    """
    Simulate ancestry and mutations under a demography specified in a demes model.

    Parameters
    ----------
    demog : str
        Demes model specification.
    nref : int
        Number of reference samples.
    ntgt : ini
        Number of target samples.
    nsrc : int
        Number of source samples.
    ref_id : str
        Population identifier in the demography for the reference population.
    tgt_id : str
        Population identifier in the demography for the target population.
    src_id : str
        Population identifier in the demography for the source population.
    seq_len : float
        Simulated sequence length.
    mut_rate : float
        Per-base mutation rate used for ``msprime.sim_mutations``.
    rec_rate : float
        Per-base recombination rate used for ``msprime.sim_ancestry``.
    seed : int
        Random seed used for both ancestry and mutation simulations.
    ploidy : int, optional
        Ploidy of samples in all populations, by default 2.

    Returns
    -------
    tskit.TreeSequence
        Simulated tree sequence with ancestry and mutations.
    """
    demo_graph = demes.load(demog)
    demography = msprime.Demography.from_demes(demo_graph)
    samples = [
        msprime.SampleSet(nref, ploidy=ploidy, population=ref_id),
        msprime.SampleSet(ntgt, ploidy=ploidy, population=tgt_id),
        msprime.SampleSet(nsrc, ploidy=ploidy, population=src_id),
    ]

    ts = msprime.sim_ancestry(
        recombination_rate=rec_rate,
        sequence_length=seq_len,
        samples=samples,
        demography=demography,
        record_migrations=True,
        random_seed=seed,
    )
    ts = msprime.sim_mutations(
        ts,
        rate=mut_rate,
        random_seed=seed,
        model=msprime.JC69(),
    )

    return ts


def create_sample_lists(
    nref: int,
    ntgt: int,
    nsrc: int,
    ref_list: str,
    tgt_list: str,
    src_list: str,
    identifier: str = "tsk",
) -> None:
    """
    Create text files with sample identifiers for each population.

    Parameters
    ----------
    nref : int
        Number of reference samples.
    ntgt : int
        Number of target samples.
    nsrc : int
        Number of source samples.
    ref_list : str
        Path to the output file containing reference sample identifiers.
    tgt_list : str
        Path to the output file containing target sample identifiers.
    src_list : str
        Path to the output file containing source sample identifiers.
    out_list : str
        Path to the output file containing outgroup sample identifiers.
    identifier : str, optional
        Prefix used for all sample identifiers, by default "tsk".
    """
    ref_range = nref
    tgt_range = ref_range + ntgt
    src_range = tgt_range + nsrc

    with open(ref_list, "w") as f:
        for i in range(ref_range):
            f.write(f"{identifier}_{i}\n")

    with open(tgt_list, "w") as f:
        for i in range(ref_range, tgt_range):
            f.write(f"{identifier}_{i}\n")

    with open(src_list, "w") as f:
        for i in range(tgt_range, src_range):
            f.write(f"{identifier}_{i}\n")


def get_true_tracts(
    ts: tskit.TreeSequence,
    tgt_id: str,
    src_id: str,
    is_phased: bool = True,
    ploidy: int = 2,
) -> str:
    """
    Extract introgressed ancestry tracts for target samples from a tree sequence.

    For all migration events between the specified source and target populations,
    this function identifies target sample nodes that are descendants of the
    migrated node in overlapping tree intervals and records the corresponding
    genomic tracts. The output is a tab-delimited string with a header line:

        Chromosome  Start  End  Sample

    where `Sample` is formatted as:
      - `tsk_{individual_id}_{hap_index}` when `is_phased=True`;
      - `tsk_{individual_id}` when `is_phased=False` (unphased, i.e., union
        across haplotypes per individual after merge).

    Parameters
    ----------
    ts : tskit.TreeSequence
        Input tree sequence containing population metadata and migration records.
    tgt_id : str
        Name of the target population (as stored in ts.populations().metadata["name"]).
    src_id : str
        Name of the source population (as stored in ts.populations().metadata["name"]).
    is_phased : bool, optional
        Whether to output haplotype-level sample identifiers. Default is True.
    ploidy : int, optional
        Ploidy used to infer haplotype index in phased mode. Default is 2.
    Returns
    -------
    str
        A tab-delimited string listing inferred introgressed tracts with columns
        Chromosome, Start, End, and Sample.
    """
    tracts = "Chromosome\tStart\tEnd\tSample\n"

    try:
        src_id = [p.id for p in ts.populations() if p.metadata["name"] == src_id][0]
    except IndexError:
        raise ValueError(f"Population {src_id} is not found.")

    try:
        tgt_id = [p.id for p in ts.populations() if p.metadata["name"] == tgt_id][0]
    except IndexError:
        raise ValueError(f"Population {tgt_id} is not found.")

    for m in ts.migrations():
        if (m.dest == src_id) and (m.source == tgt_id):
            for t in ts.trees():
                if m.left >= t.interval.right:
                    continue
                if m.right <= t.interval.left:
                    break  # [l, r)
                for n in ts.samples(tgt_id):
                    if t.is_descendant(n, m.node):
                        left = m.left if m.left > t.interval.left else t.interval.left
                        right = (
                            m.right if m.right < t.interval.right else t.interval.right
                        )
                        if is_phased:
                            hap_index = get_haplotype_index(ts, n, ploidy)
                            sample_id = f"tsk_{ts.node(n).individual}_{hap_index}"
                        else:
                            sample_id = f"tsk_{ts.node(n).individual}"
                        tracts += f"1\t{int(left)}\t{int(right)}\t{sample_id}\n"

    return tracts


with open(snakemake.output.seed_file, "w") as o:
    o.write(f"{snakemake.params.seed}\n")

ts = simulate(
    demog=snakemake.input.demes,
    nref=int(snakemake.params.n_ref),
    ntgt=int(snakemake.params.n_tgt),
    nsrc=int(snakemake.params.n_src),
    ref_id=snakemake.params.ref_id,
    tgt_id=snakemake.params.tgt_id,
    src_id=snakemake.params.src_id,
    seq_len=int(snakemake.params.length_bp),
    mut_rate=float(snakemake.params.mu),
    rec_rate=float(snakemake.params.rho),
    seed=int(snakemake.params.seed),
)

ts.dump(snakemake.output.ts)
with open(snakemake.output.vcf, "w") as o:
    # See https://github.com/tskit-dev/tskit/issues/2838
    # msprime is 0-based
    ts.write_vcf(o, allow_position_zero=True)

create_sample_lists(
    nref=int(snakemake.params.n_ref),
    ntgt=int(snakemake.params.n_tgt),
    nsrc=int(snakemake.params.n_src),
    ref_list=snakemake.output.ref_list,
    tgt_list=snakemake.output.tgt_list,
    src_list=snakemake.output.src_list,
)

true_tract_output = {
    "phased": snakemake.output.bed_phased,
    "unphased": snakemake.output.bed_unphased,
}

for phased_status in ["phased", "unphased"]:
    true_tracts = get_true_tracts(
        ts=ts,
        tgt_id=snakemake.params.tgt_id,
        src_id=snakemake.params.src_id,
        is_phased=phased_status == "phased",
        ploidy=int(snakemake.params.ploidy),
    )

    true_tracts = pr.from_string(true_tracts).merge(by="Sample")
    if true_tracts.empty:
        open(true_tract_output[phased_status], "w").close()
    else:
        true_tracts.to_csv(true_tract_output[phased_status], sep="\t", header=False)
