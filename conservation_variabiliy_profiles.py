import os
from typing import Union
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


def create_profiles(
    conservation_df: pd.DataFrame,
    path_to_save_plots: Union[str, Path],
    profile_type: str,
):
    assert (
        type(conservation_df) == pd.DataFrame
    ), "conservation_df has to be a pd.DataFrame"
    assert type(path_to_save_plots) in [
        str,
        Path,
    ], "path_to_save_plots has to be a string or Path"
    assert profile_type in [
        "amino_acid_variability",
        "conservation",
    ], "profile_type has to be wither 'amino_acid_variability' or 'conservation'"

    grouped_helices = conservation_df.groupby("helix_number")
    axes_variability = []
    for i, (g, v) in enumerate(grouped_helices):
        fig, ax = plt.subplots(figsize=(30, 15))

        ax.plot(range(len(v["NS_number"])), v[profile_type], color="black", linewidth=4)
        ax.scatter(range(len(v["NS_number"])), v[profile_type], color="black", s=300)

        conserved_ind = list(v["location_id"]).index(50)
        ax.scatter(
            conserved_ind, list(v[profile_type])[conserved_ind], color="#0542f7", s=1000
        )

        ax.set_xticks(range(len(v["NS_number"])))
        y_lab = (
            range(1, int(conservation_df[profile_type].max()) + 1, 2)
            if profile_type == "amino_acid_variability"
            else [i / 10 for i in range(0, 11, 1)]
        )
        ax.set_yticks(y_lab)

        x_ticks_labels = [
            f"{ns}-{val}"
            for ns, val in zip(v["NS_number"], v["most_frequent_amino_acid"])
        ]
        ax.set_xticklabels(x_ticks_labels, rotation=90, fontsize=50, va="baseline")
        ax.tick_params(axis="x", bottom=True, labelbottom=True, pad=190)

        ax.set_title(f"TM {i+1}", fontsize=52, loc="left")
        ax.set_yticklabels(y_lab, fontsize=50)
        ylabel = (
            "Amino acid variablitiy"
            if profile_type == "amino_acid_variability"
            else "Normalized frequency"
        )
        ax.set_ylabel(ylabel, fontsize=52)

        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)

        fig.tight_layout()
        axes_variability.append(ax)
        fig.savefig(
            os.path.join(path_to_save_plots, f"{profile_type}_profil_TM_{i+1}.png")
        )
