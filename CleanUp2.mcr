#
# Yanaconda script to Clean up a .pdb file to prepare for MD simulation with Yasara (case 1) or for Rosettaddg (case2)
#
#
LoadPDB (MacroTarget)
DelRes (DelResidues)

AllMols() = ListMol protein Obj 1, format=MOLNAME

for i = 1 to count AllMols
  RenumberRes protein Obj 1 Mol (AllMols(i)), first=(ResidueOne)
SavePDB Obj 1, (MacroTarget)_renum
DelRes !protein and !HOH
CleanAll
Opthydall
SavePDB OBJ 1, (MacroTarget)_cleaned2

AllRes()=ListRes protein Obj 1 Mol A,format=RESNUM
First=(AllRes(1))
ResAmount=Count (AllRes)
Last=(First) + (ResAmount)
Shell echo (First) > FirstLast.tab
Shell echo (Last) >> FirstLast.tab

DelRes HOH
DelRes !protein
RenameMol All, A
JoinRes protein
RenumberRes Obj 1, first=1
SavePDB Obj 1, (MacroTarget)_forRosetta2
Exit