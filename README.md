This repository contains scripts to create a custom conserved numbering scheme for protein sequence alignments.
An additional Jupyter Notebook and example files were added to illustrate the usage of the script.


## Main functionality: generate_conserved_numbering
The main generate_conserved_numbering.py creates the numbering scheme and renumbers PDB files with the following steps:

#### 1. Prepare the alignment file:
* The sequence alignment file has to be saved as a csv file in which the rows are the aligned protein sequences. With '-' indicating the gab character.
* The first row contains the structural information. The TM regions marked with 'H' and the rest as '-' (hyphen) for the entire length of the alignment, including the gaps.
* The first column contains the protein name. It can be any name or the PDB ID of a corresponding protein structure, this column will be later used to map the aligned sequences to renumber the PDB structures.

#### 2.Assigne conserved numbers:
* Call the *assign_conserved_numbers* function from the *generate_conserved_numbering.py* script with the loaded alignment csv as a Pandas DataFrame.
* The function returns a DataFrame with the assigned *Numbering Scheme* (*NS*)-numbers. This DataFrame can be saved for further analysis.

#### 3. Renumber PDB structures with conserved numbering:
* In this step the NS-numbers can be assigned to PDB structures whose sequences are in the alignment. 
* Create a mapping between the PDB files and the protein name in the sequence alignment file. In the Jupyter Notebook there is an example how to do this, but the mapping can be created in any other way.
* The same alignment can be assigned to multiple PDB structures, e.g: when there are multiple crystal structures of the same proteins or for structures with point mutations.
* Call the *renumber_pdb* function to rewrite the PDB file. Residue IDs will be replaced with the NS-prot numbers.

## Additional functionalities:
#### 1. Add the NS-prot numbers to the original PDB file:
* In this option, the original PBD file will be kept, and the NS-prot numbers are added to the B-factor column. This way the PDB file remans with the initial residue IDs next to the assigned numbers.

#### 2. Plot amino acid conservation and variability profiles:
* A separate script *conservation_variabiliy_profiles.py* contains a function to plot the most conserved amino acids and the amino acid variability profiles of the TM helices
  
![amino acid variablity profile](https://github.com/evabertalan/NS-prot/blob/main/example_files/PDBs/plots/amino_acid_variability_profil_TM_1.png)
