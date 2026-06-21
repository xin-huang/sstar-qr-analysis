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


import demes
import msprime
import numpy as np
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
        provided, the function validates that `len(individual.nodes)` equals
        this value.

    Returns
    -------
    int
        One-based position of `node_id` in `ts.individual(ind_id).nodes`.

    Raises
    ------
    ValueError
        If `node_id` has no associated individual, if the individual's node
        count does not match `expected_ploidy`, or if `node_id` is not found
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
        Per-base mutation rate used for `msprime.sim_mutations`.
    rec_rate : float
        Per-base recombination rate used for `msprime.sim_ancestry`.
    seed : int
        Random seed used for both ancestry and mutation simulations.
    ploidy : int, optional
        Ploidy of samples in all populations. Default: 2.

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
        Prefix used for all sample identifiers. Default: "tsk".
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
    this function identifies target sample nodes that are descendants of each
    migrated node in overlapping tree intervals and records the corresponding
    genomic tracts. Migration records are first filtered with tskit column arrays,
    then mapped to overlapping tree-index ranges with `numpy.searchsorted`.
    During tree traversal, an active migration set and tracked target samples are
    used to avoid scanning every tree and every target sample for each migration.

    The output is a tab-delimited string with a header line:

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
        Name of the target population, as stored in
        `ts.populations().metadata["name"]`.
    src_id : str
        Name of the source population, as stored in
        `ts.populations().metadata["name"]`.
    is_phased : bool, optional
        Whether to output haplotype-level sample identifiers. Default: True.
    ploidy : int, optional
        Ploidy used to infer haplotype index in phased mode. Default: 2.

    Returns
    -------
    str
        A tab-delimited string listing inferred introgressed tracts with columns
        Chromosome, Start, End, and Sample.

    Raises
    ------
    ValueError
        If either the source or target population name is not found, or if phased
        sample-name construction encounters inconsistent individual/ploidy data.
    """
    header = "Chromosome\tStart\tEnd\tSample\n"

    try:
        src_pop = [p.id for p in ts.populations() if p.metadata["name"] == src_id][0]
    except IndexError:
        raise ValueError(f"Population {src_id} is not found.")

    try:
        tgt_pop = [p.id for p in ts.populations() if p.metadata["name"] == tgt_id][0]
    except IndexError:
        raise ValueError(f"Population {tgt_id} is not found.")

    # Use NumPy column arrays to select matching migration records once.
    keep = (ts.migrations_dest == src_pop) & (ts.migrations_source == tgt_pop)

    if not np.any(keep):
        return header

    mig_left = ts.migrations_left[keep]
    mig_right = ts.migrations_right[keep]
    mig_node = ts.migrations_node[keep]

    target_samples = np.asarray(ts.samples(population=tgt_pop))
    target_sample_list = target_samples.tolist()

    # Use an O(1) node-indexed mask to test whether a descendant sample is in
    # the target population.
    is_target_sample = np.zeros(ts.num_nodes, dtype=bool)
    is_target_sample[target_samples] = True

    # Compute sample names once for target samples.
    sample_names = {}
    nodes_individual = ts.nodes_individual

    for node in target_sample_list:
        individual = int(nodes_individual[node])

        if is_phased:
            hap_index = get_haplotype_index(ts, node, ploidy)
            sample_names[node] = f"tsk_{individual}_{hap_index}"
        else:
            sample_names[node] = f"tsk_{individual}"

    # Map each migration interval to the range of tree indexes it overlaps.
    breakpoints = ts.breakpoints(as_array=True)

    start_tree = np.searchsorted(breakpoints, mig_left, side="right") - 1
    end_tree = np.searchsorted(breakpoints, mig_right, side="left")

    valid = start_tree < end_tree
    if not np.any(valid):
        return header

    mig_left = mig_left[valid]
    mig_right = mig_right[valid]
    mig_node = mig_node[valid]
    start_tree = start_tree[valid]
    end_tree = end_tree[valid]

    # Migration records are sorted by time, not necessarily by genomic position.
    start_order = np.argsort(start_tree, kind="stable")
    end_order = np.argsort(end_tree, kind="stable")

    start_ptr = 0
    end_ptr = 0
    num_migrations = len(mig_node)

    # Use an insertion-ordered dict as the active migration set.
    active = {}
    output = [header]

    trees = ts.trees(tracked_samples=target_sample_list, sample_lists=True)

    for tree_index, tree in enumerate(trees):
        # Remove migrations that ended before the current tree.
        while end_ptr < num_migrations and end_tree[end_order[end_ptr]] <= tree_index:
            migration_index = int(end_order[end_ptr])
            active.pop(migration_index, None)
            end_ptr += 1

        # Add migrations that become active at or before the current tree.
        while (
            start_ptr < num_migrations
            and start_tree[start_order[start_ptr]] <= tree_index
        ):
            migration_index = int(start_order[start_ptr])
            active[migration_index] = None
            start_ptr += 1

        if not active:
            continue

        tree_left = tree.interval.left
        tree_right = tree.interval.right

        for migration_index in active:
            ancestor = int(mig_node[migration_index])

            # Skip migration nodes with no tracked target-sample descendants in O(1).
            if tree.num_tracked_samples(ancestor) == 0:
                continue

            left = max(tree_left, mig_left[migration_index])
            right = min(tree_right, mig_right[migration_index])

            if left >= right:
                continue

            # Iterate over descendant samples instead of testing every target sample
            # with is_descendant.
            for sample_node in tree.samples(ancestor):
                if is_target_sample[sample_node]:
                    output.append(
                        f"1\t{int(left)}\t{int(right)}\t"
                        f"{sample_names[sample_node]}\n"
                    )

    return "".join(output)


with open(snakemake.output.seed_file, "w") as o:
    o.write(f"{snakemake.params.sim['seed']}\n")

ts = simulate(
    demog=snakemake.input.demes,
    nref=int(snakemake.wildcards.n_ref),
    ntgt=int(snakemake.wildcards.n_tgt),
    nsrc=int(snakemake.wildcards.n_src),
    ref_id=snakemake.params.sim["ref_id"],
    tgt_id=snakemake.params.sim["tgt_id"],
    src_id=snakemake.params.sim["src_id"],
    seq_len=int(snakemake.params.sim["length_bp"]),
    mut_rate=float(snakemake.params.sim["mu"]),
    rec_rate=float(snakemake.params.sim["rho"]),
    seed=int(snakemake.params.sim["seed"]),
)

ts.dump(snakemake.output.ts)
with open(snakemake.output.vcf, "w") as o:
    # See https://github.com/tskit-dev/tskit/issues/2838
    # msprime is 0-based
    ts.write_vcf(o, allow_position_zero=True)

create_sample_lists(
    nref=int(snakemake.wildcards.n_ref),
    ntgt=int(snakemake.wildcards.n_tgt),
    nsrc=int(snakemake.wildcards.n_src),
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
        tgt_id=snakemake.params.sim["tgt_id"],
        src_id=snakemake.params.sim["src_id"],
        is_phased=phased_status == "phased",
        ploidy=int(snakemake.params.sim["ploidy"]),
    )

    true_tracts = pr.from_string(true_tracts).merge(by="Sample")
    if true_tracts.empty:
        open(true_tract_output[phased_status], "w").close()
    else:
        true_tracts = pr.PyRanges(
            true_tracts.df.sort_values(["Sample", "Chromosome", "Start", "End"])
        )
        true_tracts.to_csv(true_tract_output[phased_status], sep="\t", header=False)
