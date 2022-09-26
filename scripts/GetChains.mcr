#
# Yanaconda script to determine the number of chains and their respective first and last residue number
#
#

LoadPDB (MacroTarget)

# This is the number of Molecules (protein subunits, ligands, co-factors, etc)
MoleculeNum = CountMol OBJ 1

MoleculeNames() = NameMol OBJ 1

# This is the number of protein chains (subunits)
ChainNum = CountMol OBJ 1 protein


MakeTab ChainRes, dimensions=2, columns=3

for i = 1 to (ChainNum)
  print (MoleculeNames(i))
  Reslist() = ListRes OBJ 1 MOL (MoleculeNames(i)) protein, format=RESNUM
  Tabulate (Reslist(1))
  Tabulate (Reslist(count (Reslist())))
  Tabulate
  ListMol protein OBJ 1 MOL (MoleculeNames(i)), format=MOLNAME

SaveTab ChainRes,ChainRes,NumFormat=-6.0f,Start  End    Name

for i = 1 to (ChainNum)
  MakeTab ResiduesChain(i), dimensions=2, columns=1
  SelectTab ResiduesChain(i)
  Tabulate
  ListRes protein OBJ 1 MOL (MoleculeNames(i)),format=RESNUM
  SaveTab ResiduesChain(i),ResiduesChain(MoleculeNames(i)),NumFormat=-0.0f
Exit


