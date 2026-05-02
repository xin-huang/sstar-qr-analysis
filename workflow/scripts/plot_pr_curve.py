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


from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def parse_metadata(path: str) -> dict[str, str]:
    parts = Path(path).parts
    return {
        "demog_model": parts[1],
        "n_ref": parts[2],
        "n_tgt": parts[3],
        "n_src": parts[4],
        "phase_state": parts[6],
    }


frames = []
for perf_path in snakemake.input.perf:
    meta = parse_metadata(perf_path)
    df = pd.read_csv(perf_path, sep="\t")
    for key, value in meta.items():
        df[key] = value
    frames.append(df)

all_df = pd.concat(frames, ignore_index=True)

summary = (
    all_df.groupby(["Method", "phase_state", "demog_model", "Cutoff"], as_index=False)
    .agg(
        Precision_mean=("Precision", "mean"),
        Precision_std=("Precision", "std"),
        Recall_mean=("Recall", "mean"),
        Recall_std=("Recall", "std"),
        N=("Replicate", "count"),
    )
    .sort_values(["Method", "phase_state", "demog_model", "Cutoff"])
)

summary.to_csv(snakemake.output.summary_tsv, sep="\t", index=False, na_rep="NaN")

phase_states = summary["phase_state"].drop_duplicates().tolist()
demog_models = summary["demog_model"].drop_duplicates().tolist()
methods = summary["Method"].drop_duplicates().tolist()

fig, axes = plt.subplots(
    nrows=len(phase_states),
    ncols=len(demog_models),
    figsize=(5 * len(demog_models), 4 * len(phase_states)),
    dpi=300,
    squeeze=False,
)

color_map = {method: plt.cm.tab10(i % 10) for i, method in enumerate(methods)}

for row_idx, phase_state in enumerate(phase_states):
    for col_idx, demog_model in enumerate(demog_models):
        ax = axes[row_idx][col_idx]
        panel_df = summary.query(
            "phase_state == @phase_state and demog_model == @demog_model"
        )
        for method in methods:
            curve_df = panel_df.query("Method == @method")
            if curve_df.empty:
                continue
            ax.plot(
                curve_df["Recall_mean"],
                curve_df["Precision_mean"],
                marker="o",
                markersize=3,
                linewidth=1.5,
                color=color_map[method],
                label=method,
            )

        ax.set_xlim([0, 1])
        ax.set_ylim([0, 1])
        ax.set_xlabel("Recall")
        ax.set_ylabel("Precision")
        ax.set_title(f"{phase_state} | {demog_model}")
        ax.grid(True, alpha=0.3)

handles, labels = axes[0][0].get_legend_handles_labels()
if handles:
    fig.legend(handles, labels, loc="lower center", ncol=len(methods), frameon=False)

plt.tight_layout(rect=[0, 0.06, 1, 0.95])
plt.savefig(snakemake.output.plot, bbox_inches="tight")
plt.close()
