## Description of the structural data

Unzipping the data tarball with command `tar -xzvf hiqbind.tar.gz` will yield two directories corresponding to to the "small molecule" and "polymer" subset of HiQBind.

```bash
-- raw_data_hiq_sm/
-- raw_data_hiq_poly/
```

In each of the directory, you will see a file structure like this: 
```bash
-- 1a69/
   |-- 1a69_FMB_A_240/
       |-- 1a69_FMB_A_240_ligand.pdb
       |-- 1a69_FMB_A_240_protein.pdb
       |-- 1a69_FMB_A_240_protein_hetatm.pdb
       |-- 1a69_FMB_A_240_hetatm.pdb
       |-- 1a69_FMB_A_240_ligand_refined.sdf
       |-- 1a69_FMB_A_240_protein_refined.pdb
   |-- 1a4m_FMB_B_240/
   |-- 1a4m_FMB_C_240/
-- 1a85/
```

Description of the naming conventions:

+ `1a69`: 4-letter PDB ID
+ `FMB`: Name of the ligand. If the ligand is a polymer, it will be format like "ACE-DIP", where "ACE" is the name of the first residue and "DIP" is the name of the last residue.
+ `A`: Ligand chain ID.
+ `240`: Ligand residue number. If the ligand is a polymer, it will be format like "1-3", where "1" is the residue number of the first residue and "3" is the number of the last residue. Note the residue number may contain insertion code, or be a negative integer or zero.

Description of the files:
+ `*_ligand.pdb`: ligand structure extracted from the original PDB (not processed)
+ `*_protein.pdb`: protein structure extracted from the original PDB (not processed). A protein is defined as chains within 10 angstrom of the ligand structure. 
+ `*_protein_hetatm.pdb`: protein structure with additives (solvents, ions) extracted from the original PDB (not processed). Additives are specified with "HETATM" atoms that are within 4 angstroms of the protein chains.
+ `*_hetatm.pdb`: additives' structure extracted from the original PDB (not processed)
+ `*_ligand_refined.sdf`: refined ligand structures (hydrogen added, correct bond order, better tautomer states/protonataion states) with PDBBind-Opt workflow. 
+ `*_protein_refined.pdb`: refined protein structures (hydrogen added, missing atoms/residues added) with PDBBind-Opt workflow. 


## Description of columns in the metadata csv file:
| Field                         | Type            | Description                                                                                      |
|-------------------------------|-----------------|--------------------------------------------------------------------------------------------------|
| PDBID                      | *string*        | Four-letter PDB code.                                                                            |
| Resolution                 | *string/float*  | Resolution of the structure or "NMR".                                                            |
| Year                        | *int*           | Initial deposit year.                                                                            |
| Ligand Name                 | *string*        | Ligand name. For polymers, format "ACE-DIP" (first-last residue names).                          |
| Ligand Chain                | *string*        | Chain ID.                                                                                        |
| Ligand Residue Number       | *string*        | Residue number. For polymers, format "1-3". May include insertion code or be negative/zero.       |
| Binding Affinity Measurement| *string*        | Assay type: "kd", "ki", "ic50". ("Ka"/"Kb" are converted to Kd via Ka=1/Kd.)                     |
| Binding Affinity Sign       | *string*        | Sign: "=", ">=", "<=" or "~".                                                                    |
| Binding Affinity Value      | *float*         | Affinity value.                                                                                  |
| Binding Affinity Unit       | *string*        | Unit: "fM", "pM", "nM", "uM", "mM", "M".                                                         |
| Log Binding Affinity        | *float*         | Affinity in log unit.                                                                            |
| Binding Affinity Source     | *string*        | Source: "BindingMOAD", "BindingDB", or "BioLiP".                                                 |
| Binding Affinity Annotation | *string*        | Original annotation.                                                                             |
| Protein UniProtID           | *string*        | UniProtID(s), comma-separated if multiple.                                                       |
| Protein UniProtName         | *string*        | Protein name(s), comma-separated if multiple.                                                    |
| Ligand SMILES               | *string*        | SMILES notation.                                                                                 |
| Ligand MW                   | *float*         | Molecular weight.                                                                                |
| Ligand LogP                 | *float*         | LogP value computed by RDKit.                                                                    |
| Ligand TPSA                 | *float*         | TPSA value computed by RDKit.                                                                    |
| Ligand NumRotBond           | *int*           | Number of rotatable bonds.                                                                       |
| Ligand NumHeavyAtoms        | *int*           | Number of heavy atoms.                                                                           |
| Ligand NumHDon              | *int*           | Number of hydrogen bond donors.                                                                  |
| Ligand NumHAcc              | *int*           | Number of hydrogen bond acceptors.                                                               |
| Ligand QED                  | *float*         | QED value computed by RDKit.                                                                     |
