#
# Yanaconda script to Clean up a .pdb file to prepare for MD simulation with Yasara (case 1) or for Rosettaddg (case2)
#
#

LoadPDB (MacroTarget)
DelRes !protein and !HOH
CleanAll
OptHydAll
SavePDB OBJ 1, (MacroTarget)_cleaned

DelRes !protein
DelRes (DelResidues)
SavePDB OBJ 1, (MacroTarget)_test

AllMols() = ListMol protein Obj 1, format=MOLNAME
for i = 1 to count AllMols
  RenameMol (AllMols(i)),(AllMols(1))
JoinRes protein
SavePDB OBJ 1, (MacroTarget)_forRosetta

  
Exit
