from rdkit_lib.forcefields import MMFF94Strategy, UFFStrategy
from rdkit_lib.generate_conformers import ConformerGenerator
from rdkit_lib.utils_rdkit import select_low_energy, write_xyz_files
from rdkit_lib.serializer import serializer_json
from pathlib import Path
from rdkit import Chem
import logging
import json

logging.basicConfig(
    filename='conformers_process.log',
    filemode='a',
    format='%(asctime)s - %(levelname)s - %(message)s',
    level=logging.INFO
)

path_mol = Path("./mol_data.json")
output_dir = Path("./xyz")
output_dir.mkdir(parents=True, exist_ok=True)
path_json_output = Path("conformers_backup.jsonl")

strategies = [MMFF94Strategy(),UFFStrategy()]

with open(path_mol) as json_file:
    all_mol = json.load(json_file)

with open(path_json_output, 'a') as file_output:

    for mol_data in all_mol:
        molecule_id = mol_data["id"]
        molecule = Chem.MolFromSmiles(mol_data["canonical_smiles"])

        generator = ConformerGenerator(mol=molecule)
        success = False

        for strategy in strategies:
            try:
                mol_3d, res = generator.run_optimization(strategy)
                best_energies = select_low_energy(res, max_confs=5)
                
                xyz_data = write_xyz_files(
                    mol_3d,
                    best_energies,
                    output_dir,
                    molecule_id
                    )
                
                serializer_json(
                    xyz_data=xyz_data,
                    file_output=file_output,
                    molecule_id=molecule_id,
                    strategy=strategy
                    )
                
                logging.info(f"ID {molecule_id}: Sucesso com {strategy.name}")
                success = True
                break 

            except Exception:
                continue 

        if not success:
            logging.error(f"FALHA: {molecule_id}")

with open(path_json_output, 'r') as file_input:
    final_list = [json.loads(row) for row in file_input]

with open("conformers.json", 'w') as file_output:
    json.dump(final_list, file_output, indent=4)