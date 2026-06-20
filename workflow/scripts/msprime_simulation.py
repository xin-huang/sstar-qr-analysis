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
from typing import Optional


def get_haplotype_index(
    ts: tskit.TreeSequence,
    node_id: int,
    expected_ploidy: Optional[int] = None,
) -> int:
    """
    Get 1-based haplotype index of a node within its associated individual.
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
    """
    demo_graph = demes.load(demog)
    demography = msprime.Demography.from_demes(demo_graph)

    samples = [
        msprime.SampleSet(nref, ploidy=ploidy, population=ref_id),
        msprime.SampleSet(ntgt, ploidy=ploidy, population=tgt_id),
        msprime.SampleSet(
            nsrc1,
            ploidy=ploidy,
            population=src1_id,
            time=src1_sampling_time,
        ),
    ]

    if src2_id is not None:
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
    src_id: str,
    output: str,
    is_phased: bool = True,
    ploidy: int = 2,
) -> None:
    """
    Extract true introgressed ancestry tracts for target samples from a tree sequence.

    The function writes a BED-like file with columns:

        Chromosome  Start  End  Sample

    Sample names are:
      - tsk_{individual_id}_{hap_index} when is_phased=True
      - tsk_{individual_id} when is_phased=False
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
                    break

                for n in ts.samples(tgt_id):
                    if t.is_descendant(n, m.node):
                        left = m.left if m.left > t.interval.left else t.interval.left
                        right = m.right if m.right < t.interval.right else t.interval.right

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


nsrc1 = 1
nsrc2 = 1 if snakemake.params.sim["src2_id"] is not None else 0


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
        src_id=snakemake.params.sim["src1_id"],
        output=true_tract_output["src1"][phased_status],
        is_phased=phased_status == "phased",
        ploidy=int(snakemake.params.sim["ploidy"]),
    )

    if snakemake.params.sim["src2_id"] is not None:
        get_true_tracts(
            ts=ts,
            tgt_id=snakemake.params.sim["tgt_id"],
            src_id=snakemake.params.sim["src2_id"],
            output=true_tract_output["src2"][phased_status],
            is_phased=phased_status == "phased",
            ploidy=int(snakemake.params.sim["ploidy"]),
        )
    else:
        open(true_tract_output["src2"][phased_status], "w").close()
