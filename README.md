## 1. Prepare the alignment file:
* The alignment file has to be saved as a csv file in which the rows are the aligned protein sequences. With '-' indicating the gab character.
* The first row contains the structural information. The TM regions marked with 'H' and the rest as '-' (hyphen)” for the entire length of the alignment, including the gaps.
* The first column contains the protein name. It can be any name or the PDB ID of a corresponding protein structure, this column will be later used to map the aligned sequences to renumber the PDB structures.

## 2.Assigne conserved numbers:
* Call the *assign_conserved_numbers* function from the *generate_conserved_numbering.py* script with the loaded alignment csv as a Pandas DataFrame.
* The function returns a DataFrame with the assigned *Numbering Scheme* (*NS*)-numbers. This DataFrame can be saved for furher analysis.

## 3. Renumber PDB structures with conserved numbering:
* In this step the NS-numbers can be assigned to PDB structures whose sequences are in the alignment. 
* The PDB file has to be mapped to the protein name in the alignment file. Below an example is shown, but the mapping can be created in any other way.
* The same alignment can be assigned to multiple PDB structures, e.g: when there are multiple crystal structures of the same proteins or for structures with point mutations.
