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


def get_pop_config(wildcards):
    return POP_CONFIGS[wildcards.demog_model]


def get_simulation_params(wildcards, split):
    pop = get_pop_config(wildcards)
    if split == "training":
        seed = seed_lists[split][int(wildcards.test_rep)][int(wildcards.training_rep)]
    else:
        seed = seed_lists[split][int(wildcards.test_rep)]

    return {
        "length_bp": LENGTH_BPS[split],
        "ploidy": 2,
        "seed": seed,
        "mu": pop.mut_rate,
        "rho": pop.rec_rate,
        "ref_id": pop.ref,
        "tgt_id": pop.tgt,
        "src_id": pop.src,
    }


def get_sample_size(wildcards):
    return {
        "total": 2*(int(wildcards.n_ref)+int(wildcards.n_tgt)),
        "ref": 2*(int(wildcards.n_ref)),
        "tgt": 2*(int(wildcards.n_tgt)),
    }
