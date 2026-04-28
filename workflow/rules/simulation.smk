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


rule simulate_training_data:
    input:
        demes="config/ArchIE_3D19_wo_introgression.yaml",
    output:
        ts=temp("results/simulation/training/rep_{test_rep}/simulation.rep_{training_rep}.ts"),
        vcf=temp("results/simulation/training/rep_{test_rep}/simulation.rep_{training_rep}.vcf"),
        bed_phased=temp("results/simulation/training/rep_{test_rep}/simulation.rep_{training_rep}.true.tracts.phased.bed"),
        bed_unphased=temp("results/simulation/training/rep_{test_rep}/simulation.rep_{training_rep}.true.tracts.unphased.bed"),
        ref_list=temp("results/simulation/training/rep_{test_rep}/simulation.rep_{training_rep}.ref.list"),
        tgt_list=temp("results/simulation/training/rep_{test_rep}/simulation.rep_{training_rep}.tgt.list"),
        src_list=temp("results/simulation/training/rep_{test_rep}/simulation.rep_{training_rep}.src.list"),
        seed_file=temp("results/simulation/training/rep_{test_rep}/simulation.rep_{training_rep}.seedmsprime"),
    params:
        n_ref=50,
        n_tgt=50,
        n_src=0,
        length_bp=50_000,
        mu=1.2e-8,
        rho=1e-8,
        ref_id="Reference",
        tgt_id="Target",
        src_id="Source",
        ploidy=2,
        seed=lambda wildcards: training_seed_list[int(wildcards.test_rep)][int(wildcards.training_rep)],
    resources:
        mem_gb=16,
    script:
        "../scripts/msprime_simulation.py"


rule simulate_test_data:
    input:
        demes="config/ArchIE_3D19.yaml",
    output:
        ts=temp("results/simulation/test/rep_{test_rep}/simulation.rep_{test_rep}.ts"),
        vcf=temp("results/simulation/test/rep_{test_rep}/simulation.rep_{test_rep}.vcf"),
        bed_phased=temp("results/simulation/test/rep_{test_rep}/simulation.rep_{test_rep}.true.tracts.phased.bed"),
        bed_unphased=temp("results/simulation/test/rep_{test_rep}/simulation.rep_{test_rep}.true.tracts.unphased.bed"),
        ref_list=temp("results/simulation/test/rep_{test_rep}/simulation.rep_{test_rep}.ref.list"),
        tgt_list=temp("results/simulation/test/rep_{test_rep}/simulation.rep_{test_rep}.tgt.list"),
        src_list=temp("results/simulation/test/rep_{test_rep}/simulation.rep_{test_rep}.src.list"),
        seed_file=temp("results/simulation/test/rep_{test_rep}/simulation.rep_{test_rep}.seedmsprime"),
    params:
        n_ref=50,
        n_tgt=50,
        n_src=0,
        length_bp=200_000_000,
        mu=1.2e-8,
        rho=1e-8,
        ref_id="Reference",
        tgt_id="Target",
        src_id="Source",
        ploidy=2,
        seed=lambda wildcards: test_seed_list[int(wildcards.test_rep)],
    resources:
        time=360, mem_gb=16,
    script:
        "../scripts/msprime_simulation.py"


rule extract_training_biallelic_snps:
    input:
        vcf=rules.simulate_training_data.output.vcf,
    output:
        vcf=temp("results/simulation/training/rep_{test_rep}/simulation.rep_{training_rep}.biallelic.snps.vcf.gz"),
    shell:
        """
        bcftools view {input.vcf} -v snps -m 2 -M 2 -g ^miss | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule extract_test_biallelic_snps:
    input:
        vcf=rules.simulate_test_data.output.vcf,
    output:
        vcf=temp("results/simulation/test/rep_{test_rep}/simulation.rep_{test_rep}.biallelic.snps.vcf.gz"),
    shell:
        """
        bcftools view {input.vcf} -v snps -m 2 -M 2 -g ^miss | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule calc_training_sstar_score:
    input:
        vcf=rules.simulate_training_data.output.vcf,
        ref_list=rules.simulate_training_data.output.ref_list,
        tgt_list=rules.simulate_training_data.output.tgt_list,
    output:
        score=temp("results/simulation/training/rep_{test_rep}/sstar.{phase_state}.rep_{training_rep}.scores.tsv"),
    params:
        win_len=50000,
        win_step=50000,
        phased_flag=lambda wildcards: "--phased" if wildcards.phase_state == "phased" else "",
    resources:
        mem_gb=16, cpus=4,
    conda:
        "../envs/sstar.yaml",
    shell:
        """
        sstar score \
          --vcf {input.vcf} \
          --ref {input.ref_list} \
          --tgt {input.tgt_list} \
          --output {output.score} \
          --thread {resources.cpus} \
          --win-len {params.win_len} \
          --win-step {params.win_step} \
          {params.phased_flag} \
        """


rule merge_training_sstar_score:
    input:
        scores=expand(
            "results/simulation/training/rep_{test_rep}/sstar.{phase_state}.rep_{training_rep}.scores.tsv",
            training_rep=range(TRAINING_REP),
            allow_missing=True,
        ),
    output:
        scores="results/simulation/training/rep_{test_rep}/sstar.{phase_state}.rep_{test_rep}.training.scores.tsv",
    shell:
        """
        awk 'FNR==1 && NR!=1 {{next}} {{print}}' {input.scores} > {output.scores}
        """
