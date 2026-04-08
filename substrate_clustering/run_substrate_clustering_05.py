import os
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs


input_file = '/home/users/yao.wei/Desktop/cyp_toxicogenomics/zinc/substrate_toxic_dataset.txt'

output_file = '/home/users/yao.wei/Desktop/cyp_toxicogenomics/zinc/substrate_toxic_dataset_05.txt'

# Substrate/toxicant SMILES list
substrate_smiles_list = [
    "C1=CC=C2C(=C1)C(=NN=C2NN)NN",        
    "COC1=C2C3=C(C(=O)CC3)C(=O)OC2=C4[C@@H]5C=CO[C@@H]5OC4=C1",
    "CCC1=CC=C(C=C1)C(=O)C2=CC(=C(C(=C2)O)O)[N+](=O)[O-]",
    "CC(=O)NC1=CC=C(C=C1)O",
    "CCOC1=CC=C(C=C1)NC(=O)C",
    "C1[C@@H]2[C@@H](C2N)CN1C3=C(C=C4C(=O)C(=CN(C4=N3)C5=C(C=C(C=C5)F)F)C(=O)O)F",
    "CC(=O)NC1=CC2=C(C=C1)C3=CC=CC=C3C2",
    "C1=CC=C(C=C1)C2=CC=C(C=C2)N",
    "C1=CC=C2C3=C4C(=CC2=C1)C=CC5=C4C(=CC=C5)C=C3",
    "CN1C2=CC(=C3C(=C2N=C1N)C=CC=N3)O",
    "CC(C)(C1=CC=C(C=C1)O)C2=CC=C(C=C2)O",
    "C1C[C@H](N(C1)N=O)C2=CN=CC=C2"
]

# Convert substrate SMILES to fingerprints
toxic_fps = []
for smi in substrate_smiles_list:
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
                if max_sim >= 0.5:
                    f_out.write(f"{smi}\t{zinc_id}\t{max_sim:.3f}\n")

        if i % 10000 == 0:
            print(f"Processed {i} molecules")

print(f"Finished processing. Output saved to {output_file}")
