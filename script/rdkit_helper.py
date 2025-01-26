#!/home/public/miniconda3/envs/develop/bin/python

from rdkit import Chem
from rdkit.Chem import AllChem

def rdhelper(smi):
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

    # print('sz : ', sz)
    # print('atom_num:', atom_num)
    # print('hashing:', hashing)
    # print('poly:', poly)

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
            Type = '2' if tp == Chem.rdchem.BondType.DOUBLE else ('3' if tp == Chem.rdchem.BondType.TRIPLE else ('ar' if tp == Chem.rdchem.BondType.AROMATIC else '1'))
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

# if __name__ == '__main__':
#     smiles = 'O=C(NCCCC(CCCNC(N[*])=O)CCNC(N[*])=O)N[*]'
#     atom_cache, edge_cache, poly_cache, bond_ring_cache = rdhelper(smiles)
#     print('atom_cache:', atom_cache)
#     print('edge_cache:', edge_cache)
#     print('poly_cache:', poly_cache)
#     print('bond_ring_cache:', bond_ring_cache)
