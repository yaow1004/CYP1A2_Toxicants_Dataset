import os
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs


input_file = '/home/users/yao.wei/Desktop/cyp_toxicogenomics/zinc/inhibitor_toxic_dataset.txt'


output_file = '/home/users/yao.wei/Desktop/cyp_toxicogenomics/zinc/inhibitor_toxic_dataset_09.txt'

# Inhibitor/toxicant SMILES list
inhibitor_smiles_list = [
    "CC1=NC2=C(N1)C(=O)N(C(=O)N2CC3=CC=CO3)C",        # Furafylline_i1
    "CC(C1=CC2=CC=CC=C2S1)N(C(=O)N)O", # Zileuton_i2
    "CS(=O)(=O)C1=CC=C(C=C1)C2=C(C(=O)OC2)C3=CC=CC=C3",    # rofecoxib_i3
    "COC1=C2C(=CC3=C1OC=C3)C=CC(=O)O2",  # Methoxsalen_i4
    "C1CCC2=NC3=CC=CC=C3C(=C2C1)N",  # Tacrine_i5
    "C1=CN=CC=C1C(=O)NN",  # Isoniazid_i6
    "C[C@H]1COC2=C3N1C=C(C(=O)C3=C(C(=C2N4CCN(CC4)C)F)N)C(=O)O",   # Antofloxacin_i7
    "COC1=C(C2=C(C[C@H]3C4=CC(=C(C=C4CCN3C2)OC)OC)C=C1)OC",  # Tetrahydropalmatine_i8
    "C1CC1N2C=C(C(=O)C3=CC(=C(C=C32)N4CCNCC4)F)C(=O)O",   # Ciprofloxacin_i9
    "COC1=CC(=CC2=C1OCO2)CC=C", # Myristicin_i10
]

# Convert inhibitor SMILES to fingerprints
toxic_fps = []
for smi in inhibitor_smiles_list:
    mol = Chem.MolFromSmiles(smi)
    if mol:
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, 1024)
        toxic_fps.append(fp)

with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
    f_out.write("SMILES\tZINC_ID\tSimilarity\n")
    next(f_in)  # Skip header line

    for i, line in enumerate(f_in):
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            smi, zinc_id = parts[0], parts[1]
            mol = Chem.MolFromSmiles(smi)
            if mol:
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, 1024)
                max_sim = max(DataStructs.TanimotoSimilarity(fp, tox_fp) for tox_fp in toxic_fps)
                if max_sim >= 0.9:
                    f_out.write(f"{smi}\t{zinc_id}\t{max_sim:.3f}\n")

        if i % 10000 == 0:
            print(f"Processed {i} molecules")

print(f"Finished processing. Output saved to {output_file}")
