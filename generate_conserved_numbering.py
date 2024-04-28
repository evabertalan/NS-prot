from typing import Union
import pandas as pd
import numpy as np
import os
from pathlib import Path
import MDAnalysis as _mda


def assign_conserved_numbers(alignment_df: pd.DataFrame) -> pd.DataFrame:
    assert (
        type(alignment_df) == pd.DataFrame
    ), "alignment_df, the argument of assign_conserved_numbers() has to be a Pandas DataFrame"
    # assign helix numbers
    struct = [col for col in alignment_df.columns]

    helix_numbers = {}
    is_block = False
    block = 0
    block_number = 0
    window = 4
    for i, s in enumerate(struct):
        if s.startswith("H"):
            non_helix = [
                part for part in struct[i : i + window] if part.startswith("-")
            ]

            if is_block and len(non_helix) >= 3:
                is_block = False
                block_number += 1
            helix = [part for part in struct[i : i + window] if part.startswith("H")]
            if len(helix) >= 3:
                is_block = True
                block = block_number
                block += 1
            helix_numbers.update({s: block})
        else:
            helix_numbers.update({s: None})

    # calculate most frequent values and conservation
    most_frequent_amino_acid = alignment_df.mode(axis=0).iloc[0, :].to_frame()
    occurrence = alignment_df.apply(lambda x: x.value_counts().max(), axis=0).to_frame()
    conservation = alignment_df.apply(
        lambda x: x.value_counts().max() / x.value_counts().sum(), axis=0
    ).to_frame()
    variability = alignment_df.replace({"-": None}).nunique(axis=0).to_frame()
    alignment_df = alignment_df.append(
        pd.DataFrame(data=helix_numbers, index=[0]), ignore_index=True
    )

    sequence_conserv = pd.concat(
        [
            alignment_df,
            variability.T,
            most_frequent_amino_acid.T,
            occurrence.T,
            conservation.T,
        ],
        axis=0,
    )
    sequence_conserv.iloc[-1, 0] = "conservation"
    sequence_conserv.iloc[-2, 0] = "occurrence"
    sequence_conserv.iloc[-3, 0] = "most_frequent_amino_acid"
    sequence_conserv.iloc[-4, 0] = "amino_acid_variability"
    sequence_conserv.iloc[-5, 0] = "helix_number"

    # transpose and re-index table --> PDB structures as columns, amino acid sequence as rows
    sequence_conserv.set_index("TM_helix", inplace=True)
    sequence_conserv.index.name = None
    sequence_conserv_tr = sequence_conserv.transpose().rename_axis("col")
    sequence_conserv_tr.reset_index(inplace=True)
    sequence_conserv_tr = sequence_conserv_tr.drop("col", axis=1)
    sequence_conserv_tr.loc[
        sequence_conserv_tr["helix_number"] == 0, "helix_number"
    ] = None
    sequence_conserv_tr["IDs"] = sequence_conserv_tr.index
    # If there are multiple residues with the same conservation on the helix,
    # select the one which index is closest to the middle of the helix
    middle_ids = (
        sequence_conserv_tr.groupby("helix_number")["IDs"]
        .apply(lambda x: x.iloc[(len(x) + 1) // 2])
        .to_frame()
    )
    for i, mid_id in middle_ids.iterrows():
        sequence_conserv_tr.loc[
            sequence_conserv_tr["helix_number"] == i, "middle_id"
        ] = mid_id["IDs"]

    idx = (
        sequence_conserv_tr.groupby("helix_number")["conservation"].transform(max)
        == sequence_conserv_tr["conservation"]
    )
    max_values_idx = sequence_conserv_tr[idx].index

    def select_most_conserved_closer_to_helix_middle(helix_group):
        distance_to_middle = list(
            abs(helix_group.index - helix_group["middle_id"].values[0])
        )
        return helix_group.iloc[distance_to_middle.index(min(distance_to_middle))]

    middle_most_conserved = (
        sequence_conserv_tr.iloc[max_values_idx]
        .groupby("helix_number")
        .apply(select_most_conserved_closer_to_helix_middle)
    )
    sequence_conserv_tr["location_id"] = 0
    # assign number 50 to the most conserved amino acid
    sequence_conserv_tr.loc[middle_most_conserved["IDs"], "location_id"] = 50
    reference_row_indexes = sequence_conserv_tr.groupby("helix_number")[
        "location_id"
    ].idxmax()

    def add_value(group):
        group["location_id"] = list(
            50 + group.index - reference_row_indexes[group["helix_number"].values[0]]
        )
        return group

    assigned_groups = sequence_conserv_tr.groupby("helix_number")[
        [
            "helix_number",
            "amino_acid_variability",
            "most_frequent_amino_acid",
            "occurrence",
            "conservation",
            "location_id",
        ]
    ].apply(add_value)
    sequence_conserv_tr = sequence_conserv_tr.iloc[:, :-8].join(assigned_groups)
    sequence_conserv_tr["NS_number"] = sequence_conserv_tr.apply(
        lambda row: f"{int(row['helix_number'])}.{int(row['location_id'])}"
        if row["helix_number"] > 0 and row["location_id"] > 0
        else "",
        axis=1,
    )

    return sequence_conserv_tr


def _renumber_pdb(
    sequence_conservation: pd.DataFrame,
    pdb_file: Union[str, Path],
    pdb_id: str,
    protein_name_in_alignment: str,
    folder_to_save_results: Union[str, Path],
) -> pd.DataFrame:
    assert (
        type(sequence_conservation) == pd.DataFrame
    ), "sequence_conservation has to be Pandas DataFrame"
    assert type(pdb_file) in [str, Path], "pdb_file has to be a string or Path"
    assert type(pdb_id) == str, "pdb_id has to be a string"
    assert (
        type(protein_name_in_alignment) == str
    ), "protein_name_in_alignment has to be a string"
    assert type(folder_to_save_results) in [
        str,
        Path,
    ], "folder_to_save_results has to be a string or Path"

    np.warnings.filterwarnings(
        "ignore", category=np.VisibleDeprecationWarning
    )  # the error is coming from MD Analysis, not form this script

    amino_dict = {
        "CYS": "C",
        "ASP": "D",
        "SER": "S",
        "GLN": "Q",
        "LYS": "K",
        "ILE": "I",
        "PRO": "P",
        "THR": "T",
        "PHE": "F",
        "ASN": "N",
        "GLY": "G",
        "HIS": "H",
        "LEU": "L",
        "ARG": "R",
        "TRP": "W",
        "ALA": "A",
        "VAL": "V",
        "GLU": "E",
        "TYR": "Y",
        "MET": "M",
        "HSD": "H",
        "HSE": "H",
    }

    df = sequence_conservation[
        [
            protein_name_in_alignment,
            "helix_number",
            "most_frequent_amino_acid",
            "location_id",
        ]
    ]
    protein_struct = _mda.Universe(pdb_file)

    res_dict = []
    for res in protein_struct.atoms.residues:
        if res.resname in amino_dict.keys() and res.resid >= 1:
            res_code = amino_dict[res.resname]
            res_dict.append((res_code, res.resid))
    res_range = list(range(res_dict[0][1], res_dict[-1][1]))

    # as rolling window is not supported on string in pandas, this is an alternative
    # method has to be improved but for testing ok
    seq_start_in_pdb = "".join([r[0] for r in res_dict[0:5]])
    for i, row in df.iterrows():
        if row[protein_name_in_alignment] != "-":
            start_row_index = i
            end_row_index = start_row_index + 120
            seq_part = "".join(
                df[protein_name_in_alignment]
                .iloc[start_row_index:end_row_index]
                .apply(str)
            ).replace("-", "")
            if seq_part[0:5] == seq_start_in_pdb:
                start_row_id = i
                break

    def match_res_id_with_row(row, protein_name_in_alignment, res_range):
        if (
            row[protein_name_in_alignment] != "-"
            and row[protein_name_in_alignment] != "1"
            and len(res_range)
        ):
            row["res_id"] = res_range.pop(0)
        else:
            row["res_id"] = "-"
        return row

    pdb_sequence_df = df.iloc[start_row_id:].apply(
        lambda r: match_res_id_with_row(r, protein_name_in_alignment, res_range), axis=1
    )
    pdb_sequence_df["NS_number"] = pdb_sequence_df.apply(
        lambda row: f"{int(row['helix_number'])}.{int(row['location_id'])}"
        if row["helix_number"] > 0 and row["location_id"] > 0
        else "",
        axis=1,
    )
    pdb_sequence_df.to_csv(
        os.path.join(folder_to_save_results, f"{pdb_id.lower()}_renumbered.csv")
    )

    pdb_sequence_df["helix_number"] = (
        pdb_sequence_df["helix_number"].fillna(0).astype(int)
    )
    return pdb_sequence_df


def renumber_pdb(
    sequence_conservation: pd.DataFrame,
    pdb_file: Union[str, Path],
    pdb_id: str,
    protein_name_in_alignment: str,
    folder_to_save_results: Union[str, Path],
) -> None:
    pdb_sequence_df = _renumber_pdb(
        sequence_conservation,
        pdb_file,
        pdb_id,
        protein_name_in_alignment,
        folder_to_save_results,
    )
    protein_struct = _mda.Universe(pdb_file)
    og_struct = _mda.Universe(pdb_file)
    for i, res in enumerate(og_struct.atoms):
        row = pdb_sequence_df.loc[pdb_sequence_df["res_id"] == res.resid]

        if len(row) and row["helix_number"].values[0] in [1, 2, 3, 4, 5, 6, 7]:
            NS_number = f"{row['helix_number'].to_numpy()[0]}{int(row['location_id'].to_numpy()[0])}"
            protein_struct.atoms[i].residue.resid = int(NS_number)
            protein_struct.atoms[i].residue.resname = "BWX"

    protein_struct.atoms.write(
        os.path.join(folder_to_save_results, f"{pdb_id.lower()}_renum.pdb")
    )


def add_NS_numbers_to_pdb(
    sequence_conservation: pd.DataFrame,
    pdb_file: Union[str, Path],
    pdb_id: str,
    protein_name_in_alignment: str,
    folder_to_save_results: Union[str, Path],
) -> None:
    pdb_sequence_df = _renumber_pdb(
        sequence_conservation,
        pdb_file,
        pdb_id,
        protein_name_in_alignment,
        folder_to_save_results,
    )
    protein_struct = _mda.Universe(pdb_file)

    for i, res in enumerate(protein_struct.atoms):
        row = pdb_sequence_df.loc[pdb_sequence_df["res_id"] == res.resid]
        if len(row) and row["helix_number"].values[0] in [1, 2, 3, 4, 5, 6, 7]:
            NS_number = f"{row['helix_number'].to_numpy()[0]}.{int(row['location_id'].to_numpy()[0])}"
            res.tempfactor = float(NS_number)
        else:
            res.tempfactor = 0

    protein_struct.atoms.write(
        os.path.join(folder_to_save_results, f"{pdb_id.lower()}_with_NS.pdb")
    )
