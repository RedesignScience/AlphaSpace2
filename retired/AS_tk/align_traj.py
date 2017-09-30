"""
Usage:

python align_traj.py reference_pdb_file target_folder


target_folder is optional, if left blank will use current working directory


"""

import os
from ..PDB import PDBParser, PDBIO,Superimposer


parser = PDBParser(QUIET=True)
io = PDBIO()
imposer = Superimposer()

if __name__ == '__main__':
    import sys,os

    reference_file = sys.argv[1]
    try:
        target_folder = sys.argv[2]
    except:
        target_folder = os.getcwd()

reference_structure = parser.get_structure('ref',reference_file)


def get_protein_bb(structure):
    protein = [m for m in structure.get_models()][0]
    if len([m for m in structure.get_models()]) > 1:
        for model in structure.get_models():
            if len([residue for residue in model.get_residues()]) > len(
                    [residue for residue in protein.get_residues()]):
                protein = model
    return sorted([atom for atom in protein.get_atoms() if atom.backbone],key=lambda a: a.get_serial_number())


ref_bb = get_protein_bb(reference_structure)

for traj_dir in [os.path.join(target_folder,d) for d in os.listdir(target_folder) if
                 os.path.isdir(os.path.join(target_folder,d)) and d[:5] == "Prod_"]:
    for traj_file in [os.path.join(traj_dir,d) for d in os.listdir(traj_dir) if d[:8] == "prod.pdb"]:
        new_file_name = traj_file[:traj_file.find('prod.pdb') + 5] + traj_file[traj_file.find('prod.pdb') + 9:] + '.pdb'
        target_structure = parser.get_structure(id='pro',file=traj_file)
        target_protein_bb = get_protein_bb(target_structure)
        imposer.set_atoms(fixed=ref_bb,moving=target_protein_bb)
        imposer.apply(target_structure)
        for model in target_structure.get_models():
            io.set_structure(model)
            io.save(file=new_file_name,append=True,write_ter=True)
        os.remove(traj_file)
    print(traj_dir)
