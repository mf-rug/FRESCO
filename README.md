Please also see the [discussion forum](https://groups.google.com/g/fresco-stabilization-of-proteins) for FRESCO and the [official download](https://www.dropbox.com/sh/j0tlsi2d3lg3jx2/AADZfPvSKeucfBtlWt44wZYBa?dl=0) of the main scripts.

TicDCRiy is an attempt to largely automate FRESCO, some parts (especially the submission to the high-performance computer cluster) will have to be modified for the user's particular environment. The additional scripts on which TicDCRiy relies are in this github repository. The folder structure from the TicDCRiy script is identical to the one you get if you follow all the commands in [our paper](https://link.springer.com/protocol/10.1007/978-1-4939-7366-8_5) and is essential to make the YASARA plugin for the mutant inspection work. The directory structure should look like this:

```
frescoXXXX
  |
  |___designsMD
        |
        |___XXXX_cleaned.pdb
        |
        |___list_SelectedMutations.tab
        |
        |___NamedPdbFiles
        |      |
        |      |___Subdir_XNY
        |      |      |
        |      |      |___XXXX_cleaned_XNY_WOW.pdb
        |      |      |
        |      |      |___XXXX_cleaned_XNY_WOW_N_20000fs_50000fs_01_LSOn_Avg.yob            # (N=1 to 5)
        |      |      |
        |      |      ...
        |      |
        |      |
        |      |___Subdir_template
        |      |
        |      ...
        |
        |
        |___UniqueDisulfides
               |
               |___Subdir_XNY
               |
               |___Subdir_templates
               |
               ...
```
