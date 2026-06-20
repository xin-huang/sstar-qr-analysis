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
import tskit
import pyranges as pr
from typing import Optional, Union


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
    ref_id: str,
    tgt_id: str,
    src1_id: str,
    src2_id: str,
    seq_len: int,
    mut_rate: float,
    rec_rate: float,
    seed: int,
    ploidy: int = 2,
    nsrc1: int = 1,
    nsrc2: int = 1,
    src1_sampling_time: float = None,
    src2_sampling_time: float = None,
) -> tskit.TreeSequence:
    """
    Simulate ancestry and mutations under a demography specified in a demes model.

    Parameters
    ----------
    demog : str
        Demes model specification.
    nref : int
        Number of reference samples.
    ntgt : int
        Number of target samples.
    ref_id : str
        Population identifier in the demography for the reference population.
    tgt_id : str
        Population identifier in the demography for the target population.
    src1_id : str
        Population identifier in the demography for the first source population.
    src2_id : str, optional
        Population identifier in the demography for the optional second source
        population. If None, only the first source population is sampled.
    seq_len : int
        Simulated sequence length.
    mut_rate : float
        Per-base mutation rate used for ``msprime.sim_mutations``.
    rec_rate : float
        Per-base recombination rate used for ``msprime.sim_ancestry``.
    seed : int
        Random seed used for both ancestry and mutation simulations.
    ploidy : int, optional
        Ploidy of samples in all populations, by default 2.
    nsrc1 : int, optional
        Number of samples from the first source population, by default 1.
    nsrc2 : int, optional
        Number of samples from the second source population, by default 1.
    src1_sampling_time : float, optional
        Sampling time for the first source population.
    src2_sampling_time : float, optional
        Sampling time for the second source population.

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
    ]

    if nsrc1 > 0:
        samples.append(
            msprime.SampleSet(
                nsrc1,
                ploidy=ploidy,
                population=src1_id,
                time=src1_sampling_time,
            )
        )

    if src2_id is not None and nsrc2 > 0:
        samples.append(
            msprime.SampleSet(
                nsrc2,
                ploidy=ploidy,
                population=src2_id,
                time=src2_sampling_time,
            )
        )

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
    ref_list: str,
    tgt_list: str,
    src_list: str,
    src2_list: str = None,
    nsrc1: int = 1,
    nsrc2: int = 1,
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
    ref_list : str
        Path to the output file containing reference sample identifiers.
    tgt_list : str
        Path to the output file containing target sample identifiers.
    src_list : str
        Path to the output file containing first source sample identifiers.
    src2_list : str, optional
        Path to the output file containing second source sample identifiers.
        If None, no second source list is written.
    nsrc1 : int, optional
        Number of samples from the first source population, by default 1.
    nsrc2 : int, optional
        Number of samples from the second source population, by default 1.
    identifier : str, optional
        Prefix used for all sample identifiers, by default "tsk".
    """
    ref_range = nref
    tgt_range = ref_range + ntgt
    src1_range = tgt_range + nsrc1
    src2_range = src1_range + nsrc2

    with open(ref_list, "w") as f:
        for i in range(ref_range):
            f.write(f"{identifier}_{i}\n")

    with open(tgt_list, "w") as f:
        for i in range(ref_range, tgt_range):
            f.write(f"{identifier}_{i}\n")

    with open(src_list, "w") as f:
        for i in range(tgt_range, src1_range):
            f.write(f"{identifier}_{i}\n")

    if src2_list is not None:
        with open(src2_list, "w") as f:
            for i in range(src1_range, src2_range):
                f.write(f"{identifier}_{i}\n")


def get_true_tracts(
    ts: tskit.TreeSequence,
    tgt_id: str,
    src_id: Union[str, tuple[str, ...], list[str]],
    output: str,
    is_phased: bool = True,
    ploidy: int = 2,
) -> None:
    """
    Extract introgressed ancestry tracts for target samples from a tree sequence.

    For all migration events between the specified source population(s) and
    target population, this function identifies target sample nodes that are
    descendants of the migrated node in overlapping tree intervals and records
    the corresponding genomic tracts. The output is written as a BED-like
    tab-delimited file with the columns:

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
        ``ts.populations().metadata["name"]``.
    src_id : str or tuple/list of str
        Name(s) of the source population(s), as stored in
        ``ts.populations().metadata["name"]``. Multiple source populations are
        merged into one output BED.
    output : str
        Path to the output BED-like file.
    is_phased : bool, optional
        Whether to output haplotype-level sample identifiers. Default is True.
    ploidy : int, optional
        Ploidy used to infer haplotype index in phased mode. Default is 2.

    Raises
    ------
    ValueError
        If any source or target population name is not found in the tree sequence.
    """
    src_names = (src_id,) if isinstance(src_id, str) else tuple(src_id)
    tracts = "Chromosome\tStart\tEnd\tSample\n"

    try:
        tgt_pop_id = [p.id for p in ts.populations() if p.metadata["name"] == tgt_id][0]
    except IndexError:
        raise ValueError(f"Population {tgt_id} is not found.")

    for src_name in src_names:
        try:
            src_pop_id = [
                p.id for p in ts.populations() if p.metadata["name"] == src_name
            ][0]
        except IndexError:
            raise ValueError(f"Population {src_name} is not found.")

        for m in ts.migrations():
            if (m.dest == src_pop_id) and (m.source == tgt_pop_id):
                for t in ts.trees():
                    if m.left >= t.interval.right:
                        continue
                    if m.right <= t.interval.left:
                        break  # [l, r)

                    for n in ts.samples(tgt_pop_id):
                        if t.is_descendant(n, m.node):
                            left = max(m.left, t.interval.left)
                            right = min(m.right, t.interval.right)

                            if is_phased:
                                hap_index = get_haplotype_index(ts, n, ploidy)
                                sample_id = f"tsk_{ts.node(n).individual}_{hap_index}"
                            else:
                                sample_id = f"tsk_{ts.node(n).individual}"

                            tracts += f"1\t{int(left)}\t{int(right)}\t{sample_id}\n"

    true_tracts = pr.from_string(tracts).merge(by="Sample")

    if true_tracts.empty:
        open(output, "w").close()
    else:
        true_tracts.to_csv(output, sep="\t", header=False)


with open(snakemake.output.seed_file, "w") as o:
    o.write(f"{snakemake.params.sim['seed']}\n")


if snakemake.params.sim["src2_id"] is None:
    nsrc1 = 0
    nsrc2 = 0
else:
    nsrc1 = 1
    nsrc2 = 1


ts = simulate(
    demog=snakemake.input.demes,
    nref=int(snakemake.wildcards.n_ref),
    ntgt=int(snakemake.wildcards.n_tgt),
    ref_id=snakemake.params.sim["ref_id"],
    tgt_id=snakemake.params.sim["tgt_id"],
    src1_id=snakemake.params.sim["src1_id"],
    src2_id=snakemake.params.sim["src2_id"],
    seq_len=int(snakemake.params.sim["length_bp"]),
    mut_rate=float(snakemake.params.sim["mu"]),
    rec_rate=float(snakemake.params.sim["rho"]),
    seed=int(snakemake.params.sim["seed"]),
    ploidy=int(snakemake.params.sim["ploidy"]),
    nsrc1=nsrc1,
    nsrc2=nsrc2,
    src1_sampling_time=snakemake.params.sim["src1_sampling_time"],
    src2_sampling_time=snakemake.params.sim["src2_sampling_time"],
)


ts.dump(snakemake.output.ts)

with open(snakemake.output.vcf, "w") as o:
    # See https://github.com/tskit-dev/tskit/issues/2838
    # msprime is 0-based
    ts.write_vcf(o, allow_position_zero=True)


create_sample_lists(
    nref=int(snakemake.wildcards.n_ref),
    ntgt=int(snakemake.wildcards.n_tgt),
    ref_list=snakemake.output.ref_list,
    tgt_list=snakemake.output.tgt_list,
    src_list=snakemake.output.src_list,
    src2_list=snakemake.output.src2_list,
    nsrc1=nsrc1,
    nsrc2=nsrc2,
)


if hasattr(snakemake.output, "bed_phased"):
    for phased_status, output_bed in {
        "phased": snakemake.output.bed_phased,
        "unphased": snakemake.output.bed_unphased,
    }.items():
        get_true_tracts(
            ts=ts,
            tgt_id=snakemake.params.sim["tgt_id"],
            src_id=snakemake.params.sim["src1_tract_sources"],
            output=output_bed,
            is_phased=phased_status == "phased",
            ploidy=int(snakemake.params.sim["ploidy"]),
        )
else:
    true_tract_output = {
        "src1": {
            "phased": snakemake.output.bed_src1_phased,
            "unphased": snakemake.output.bed_src1_unphased,
        },
        "src2": {
            "phased": snakemake.output.bed_src2_phased,
            "unphased": snakemake.output.bed_src2_unphased,
        },
    }

    for phased_status in ["phased", "unphased"]:
        get_true_tracts(
            ts=ts,
            tgt_id=snakemake.params.sim["tgt_id"],
            src_id=snakemake.params.sim["src1_tract_sources"],
            output=true_tract_output["src1"][phased_status],
            is_phased=phased_status == "phased",
            ploidy=int(snakemake.params.sim["ploidy"]),
        )
        get_true_tracts(
            ts=ts,
            tgt_id=snakemake.params.sim["tgt_id"],
            src_id=snakemake.params.sim["src2_tract_sources"],
            output=true_tract_output["src2"][phased_status],
            is_phased=phased_status == "phased",
            ploidy=int(snakemake.params.sim["ploidy"]),
        )
