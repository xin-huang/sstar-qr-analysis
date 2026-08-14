from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def parse_metadata(path: str) -> dict[str, str]:
    parts = Path(path).parts
    idx = parts.index("performance")

    return {
        "demog_model": parts[idx + 1],
        "n_ref": parts[idx + 2],
        "n_tgt": parts[idx + 3],
        "n_src": parts[idx + 4],
        "phase_state": parts[idx + 5],
        "source_name": parts[idx + 6],
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
    all_df.groupby(
        [
            "Method",
            "source_name",
            "phase_state",
            "demog_model",
            "Cutoff",
        ],
        as_index=False,
    )
    .agg(
        Precision_mean=("Precision", "mean"),
        Precision_std=("Precision", "std"),
        Recall_mean=("Recall", "mean"),
        Recall_std=("Recall", "std"),
        N=("Replicate", "count"),
    )
    .sort_values(
        [
            "Method",
            "source_name",
            "phase_state",
            "demog_model",
            "Cutoff",
        ]
    )
)


summary.to_csv(
    snakemake.output.summary_tsv,
    sep="\t",
    index=False,
    na_rep="NaN",
)


phase_states = list(dict.fromkeys(all_df["phase_state"]))
demog_models = list(dict.fromkeys(all_df["demog_model"]))

methods = ["sstar", "sstar2"]
source_names = ["src1", "src2"]


fig, axes = plt.subplots(
    nrows=len(phase_states),
    ncols=len(demog_models),
    figsize=(5 * len(demog_models), 4 * len(phase_states)),
    dpi=300,
    squeeze=False,
)


color_map = {
    "sstar": plt.cm.tab10(0),
    "sstar2": plt.cm.tab10(1),
}

linestyle_map = {
    "src1": "-",
    "src2": "--",
}


for row_idx, phase_state in enumerate(phase_states):

    for col_idx, demog_model in enumerate(demog_models):

        ax = axes[row_idx][col_idx]

        panel_df = summary.query(
            "phase_state == @phase_state and "
            "demog_model == @demog_model"
        )

        for method in methods:

            for source_name in source_names:

                curve_df = panel_df.query(
                    "Method == @method and "
                    "source_name == @source_name"
                )

                if curve_df.empty:
                    continue

                color = color_map[method]

                if source_name == "src1":
                    marker_facecolor = color
                else:
                    marker_facecolor = "white"

                ax.plot(
                    curve_df["Recall_mean"],
                    curve_df["Precision_mean"],
                    marker="o",
                    markersize=4,
                    linewidth=1.5,
                    color=color,
                    linestyle=linestyle_map[source_name],
                    markerfacecolor=marker_facecolor,
                    markeredgecolor=color,
                    markeredgewidth=1,
                    label=f"{method} {source_name}",
                )

        ax.set_xlim([0, 1])
        ax.set_ylim([0, 1])

        ax.set_xlabel("Recall")
        ax.set_ylabel("Precision")

        ax.set_title(
            f"{phase_state} | {demog_model}"
        )

        ax.grid(True, alpha=0.3)


handles, labels = axes[0][0].get_legend_handles_labels()

if handles:
    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=4,
        frameon=False,
    )


plt.tight_layout(
    rect=[0, 0.07, 1, 1]
)

plt.savefig(
    snakemake.output.plot,
    bbox_inches="tight",
)

plt.close()
