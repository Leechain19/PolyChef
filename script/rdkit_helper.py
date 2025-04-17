#!/home/public/miniconda3/envs/develop/bin/python

from rdkit import Chem
from rdkit.Chem import AllChem

SINGLE_BOND = "1"
DOUBLE_BOND = "2"
TRIPLE_BOND = "3"
AROMATIC_BOND = "ar"

def read_smi(smi):
    mol = Chem.MolFromSmiles(smi)
    # convert to 3D
    m2 = Chem.AddHs(mol)
    AllChem.EmbedMolecule(m2)
    AllChem.MMFFOptimizeMolecule(m2)
    # AllChem.UFFOptimizeMolecule(m2)
    # m3 = Chem.RemoveHs(m2)
    m3 = m2
    atom_num, bond_num = m3.GetNumAtoms(), m3.GetNumBonds() # no Hydrogen atoms connected with 13C

    # get Atom and process with C13
    # poly 0: nope 1: poly
    poly = [0] * atom_num

    for i in range(atom_num):
        atom = m3.GetAtomWithIdx(i)
        if atom.GetIsotope() == 13:
            poly[i] = 1

    # make hashing
    hashing = [-1] * atom_num
    sz = 0
    for i in range(atom_num):
        if poly[i] == 0: # normal atom
            hashing[i] = sz
            sz += 1

    atom_cache = [tuple() for _ in range(sz)]
    edge_cache = []
    poly_cache = []
    bond_ring_cache = []

    # get atom
    for i, atom in enumerate(m3.GetAtoms()):
        if poly[i]: continue
        symbol = atom.GetSymbol()
        x, y, z = m3.GetConformer().GetAtomPosition(i)
        atom_cache[hashing[i]] = (symbol, x, y, z, int(atom.GetIsAromatic()))
        # g.add_atom(Atom(atom=symbol, x=x, y=y, z=z), mono=0, idx = hashing[i], ar=atom.GetIsAromatic())

    # get bond
    for i, bond in enumerate(m3.GetBonds()):
        start, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if poly[start] or poly[end]:
            if poly[end]:
                start, end = end, start
            x, y, z = m3.GetConformer().GetAtomPosition(start)
            assert hashing[end] >= 0
            # g.add_poly(x, y, z, hashing[end])
            poly_cache.append((hashing[end], x, y, z))
        else:
            tp = bond.GetBondType()
            Type = DOUBLE_BOND if tp == Chem.rdchem.BondType.DOUBLE else (TRIPLE_BOND if tp == Chem.rdchem.BondType.TRIPLE else (AROMATIC_BOND if tp == Chem.rdchem.BondType.AROMATIC else SINGLE_BOND))
            # g.add_edge(From=hashing[start], To=hashing[end], Type=Type)
            # g.add_edge(From=hashing[end], To=hashing[start], Type=Type)
            edge_cache.append((hashing[start], hashing[end], Type))

            if bond.IsInRing():
                bond_ring_cache.append((hashing[start], hashing[end]))
                # g.add_ring_edge(From=hashing[start], To=hashing[end])

    # print('sz: ', sz)
    # print('atom_cache:', atom_cache)
    # print('poly_cache: ', poly_cache)
    return atom_cache, edge_cache, poly_cache, bond_ring_cache

def parse_mol2_file(mol2_file_path) -> dict:
    poly_dict = {}
    in_atom_section = False

    with open(mol2_file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('@<TRIPOS>ATOM'):
                in_atom_section = True
                continue
            if line.startswith('@<TRIPOS>BOND'):
                in_atom_section = False
                break
            if in_atom_section:
                parts = line.split()
                atom_index = int(parts[0]) - 1
                atom_name = parts[1]
                # x = float(parts[2])
                # y = float(parts[3])
                # z = float(parts[4])
                # atom_type = parts[5]
                # residue_number = int(parts[6])
                # residue_name = parts[7]
                # charge = float(parts[8])
                if atom_name.startswith('*'):
                    poly_dict[atom_index] = atom_name

    return poly_dict


def read_mol2(path):
    m3 = Chem.MolFromMol2File(path, removeHs = False)
    atom_num, bond_num = m3.GetNumAtoms(), m3.GetNumBonds() # no Hydrogen atoms connected with 13C
    poly_dict = parse_mol2_file(path)

    # get Atom and process with *
    # poly 0: nope 1: poly
    poly = [0] * atom_num

    for i in range(atom_num):
        atom = m3.GetAtomWithIdx(i)
        if atom.GetSymbol() == '*':
            poly[i] = 1

    # make hashing
    hashing = [-1] * atom_num
    sz = 0
    for i in range(atom_num):
        if poly[i] == 0: # normal atom
            hashing[i] = sz
            sz += 1

    atom_cache = [tuple() for _ in range(sz)]
    edge_cache = []
    poly_cache = []
    bond_ring_cache = []

    # get atom
    for i, atom in enumerate(m3.GetAtoms()):
        if poly[i]: continue
        symbol = atom.GetSymbol()
        x, y, z = m3.GetConformer().GetAtomPosition(i)
        atom_cache[hashing[i]] = (symbol, x, y, z, int(atom.GetIsAromatic()))
        # g.add_atom(Atom(atom=symbol, x=x, y=y, z=z), mono=0, idx = hashing[i], ar=atom.GetIsAromatic())

    # get bond
    for i, bond in enumerate(m3.GetBonds()):
        start, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if poly[start] or poly[end]:
            if poly[end]:
                start, end = end, start
            x, y, z = m3.GetConformer().GetAtomPosition(start)
            assert hashing[end] >= 0
            poly_cache.append((start, (hashing[end], x, y, z)))
        else:
            tp = bond.GetBondType()
            Type = '2' if tp == Chem.rdchem.BondType.DOUBLE else ('3' if tp == Chem.rdchem.BondType.TRIPLE else ('4' if tp == Chem.rdchem.BondType.AROMATIC else '1'))
            edge_cache.append((hashing[start], hashing[end], Type))

            if bond.IsInRing():
                bond_ring_cache.append((hashing[start], hashing[end]))

    poly_cache = [y[1] for y in sorted(poly_cache, key=lambda x : poly_dict[x[0]])]

    return atom_cache, edge_cache, poly_cache, bond_ring_cache

if __name__ == '__main__':
    atom_cache, edge_cache, poly_cache, bond_ring_cache = read_mol2('inputs/test_mol.mol2')
    print('atom:', atom_cache)
    print('edge:', edge_cache)
    print('poly:', poly_cache)
    print('bond_ring:', bond_ring_cache)