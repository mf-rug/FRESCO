# YASARA PLUGIN
# TOPIC:       Thermostable Mutations
# TITLE:       Inspect Mutated Residues suggested by FRESCO
# AUTHOR:      Maximilian Fuerst
# LICENSE:     GPL
# DESCRIPTION: This Plugin can be used to conveniently go through all the Mutations suggested by the FRESCO software
             # (which uses FoldX and Rosettaddg algorithms and a Yasara script for disulfide bond discovery),
             # which is intended to predict mutations that increase thermostability of a protein
#
# This is a YASARA plugin to be placed in the yasara/plg subdirectory
# Go to www.yasara.org/plugins for documentation and downloads
#

"""
MainMenu: Analyze
  PullDownMenu after Superpose: FRESC_O_
    SubMenu: _S_tart Inspection of Mutants
      Request: All
    SubMenu: Prepare _E_xcel file from Mutations list
      FileSelectionWindow: Select tab file with predicted Mutations and Energies
        MultipleSelections: No
        Filename: ./*.tab
      Request: Excel
"""

# import the yasara module (which is in the yasara/plg/ folder, therefore this script has to be there too)
import yasara, os, sys, re, shutil

global Amino_acid_dict
Amino_acid_dict = {"G": ("Gly", -0.4), "A": ("Ala", 1.8), "V": ("Val", 4.2), "L": ("Leu", 3.8), "I": ("Ile", 4.5),
                   "M": ("Met", 1.9), "W": ("Trp", -0.9), "F": ("Phe", 2.8), "P": ("Pro", -1.6), "S": ("Ser", -0.8),
                   "T": ("Thr", -0.7), "C": ("Cys", 2.5), "Y": ("Tyr", -1.3), "N": ("Asn", -3.5), "Q": ("Gln", -3.5),
                   "D": ("Asp", -3.5), "E": ("Glu", -3.5), "K": ("Lys", -3.9), "R": ("Arg", -4.5), "H": ("His", -3.2)}


# A function to add a slash (Unix) or backslash (Windows) to directory paths
def add_sep(string):
    return string if string.endswith(os.sep) else string + os.sep


# A function to get information like residue number from a string
def fresco_grab(string,request="residueplus"):
    res_search = re.findall("(?:^|_)([A-Z](\d+)[A-Z])(?:_|-)",string)
    if request == "mutation":
        return res_search[0][0]
    if request == "residueplus":
        return '{0:04}'.format(int(string.split("-")[0][1:-1])) + string.split("-")[0][-1]
    if request == "residue":
        try:
            return int(res_search[0][1])
        except IndexError:
            try:
                return int(re.sub(".*[A-Z](\d*)[A-Z].*",r"\1",string))
            except:
                return None
    if request == "disulfide":
        item = re.findall("((?:(?:\d)*ps_)?[A-Z]_[A-Z]\d*[A-Z]_[A-Z]_[A-Z]\d*[A-Z]_[^_]*)", string)[0]
        if item.split("_")[0].endswith("ps"):
            return item.split("_")[2] + "-" + item.split("_")[4] + ", SU " + item.split("_")[1] + "-" + item.split("_")[3] + ", at MD " + item.split("_")[0] + ", conf: " + item.split("_")[-1]
        else:
            return item.split("_")[1] + "-" + item.split("_")[3] + ", SU " + item.split("_")[0] + "-" + item.split("_")[2]  + ", conf: " + item.split("_")[-1]


# If the user starts any of the FRESCO menu choices, this will be initiated
if yasara.request != "CheckIfDisabled":

    # Path and file definitons
    # ------------------------------------------------------

    fresco_skip_uniques = 0

    # try to find the correct folder automatically by checking yasara's current working directory
    if ("fresco" in yasara.workdir or "designsMD" in yasara.workdir) and not len(re.findall("fresco", yasara.workdir)) > 1:
        fresco_design_dir = add_sep(re.sub("fresco([^\\\/]*).*", "fresco\\1" + os.sep + "designsMD", yasara.workdir))
    # if that fails, ask the user for input where the files are
    else:
        fresco_design_dir = add_sep(os.path.dirname(yasara.ShowWin("FileSelection","Select any file in the designsMD folder of your FRESCO files", "No","*.pdb")[1]))    # Check if the _cleaned.pdb file exists, exit otherwise
    if os.path.isfile(add_sep(fresco_design_dir) + [file for file in os.listdir(fresco_design_dir) if file.endswith("_cleaned.pdb")][0]) is True:
        Fresco_wildtype_file = [file for file in os.listdir(fresco_design_dir) if file.endswith("_cleaned.pdb")][0]
    else:
        yasara.RaiseError("Couldn't find a '_cleaned.pdb' file in this folder. Exiting")
        yasara.plugin.end()

    print "Designated folder for designs: " + fresco_design_dir

    # Define and check the subfolders
    fresco_namedpdbs_dir = fresco_design_dir + "NamedPdbFiles" + os.sep
    if not os.path.exists(fresco_namedpdbs_dir):
        yasara.RaiseError("Couldn't find the NamedPdbFiles folder. Exiting.")
        yasara.plugin.end()

    fresco_unique_dsb_dir = fresco_design_dir + "UniqueDisulfides" + os.sep
    if not os.path.exists(fresco_unique_dsb_dir):
        print '################################Information######################################'
        print '#################################################################################'
        print "Couldn't find the UniqueDisulfides folder. Going to continue with single mutants."
        print '#################################################################################'
        print '#################################################################################'
        fresco_skip_uniques = 1
        
    # Loop through the Subfolders to find all mutants and disulfide bonds, store them in a dictionary format: { Mutant: Path to mutant yob file, ... }
    Fresco_yob_dict = {}
    # First check the NamedPdbFiles folder, since os.listdir includes files, check that only Subdirectories are analyzed
    for folder in sorted([item for item in os.listdir(fresco_namedpdbs_dir) if (os.path.isdir(os.path.join(fresco_namedpdbs_dir,item)) and item.startswith("Subdir"))]):
        # Make sure the folder contains all the MD files by confirming the existance of the fifth .yob file
        file = [item for item in os.listdir(os.path.join(fresco_namedpdbs_dir, folder)) if item.endswith("05_LSOn_Avg.yob")]
        # ignore the template folder
        if len(file) > 0 and "template" not in folder:
            file = file[0]
            # Add the mutation and the path to the dictionary
            Fresco_yob_dict[fresco_grab(file, "mutation")] = os.path.join(fresco_namedpdbs_dir,folder) + os.sep + file

    if fresco_skip_uniques != 1:
        for folder in sorted([item for item in os.listdir(fresco_unique_dsb_dir) if os.path.isdir(os.path.join(fresco_unique_dsb_dir,item))]):
            file = [item for item in os.listdir(os.path.join(fresco_unique_dsb_dir, folder)) if item.endswith("05_LSOn_Avg.yob")]
            if len(file) > 0 and "template" not in folder:
                for item in file:
                    Fresco_yob_dict[fresco_grab(item, "disulfide")] = os.path.join(fresco_unique_dsb_dir,folder) + os.sep + item

    if Fresco_yob_dict == {}:
        yasara.RaiseError("Couldn't find the .yob files of the mutant's MD simulations. Was looking in " + fresco_namedpdbs_dir + ". Exiting")
        yasara.plugin.end()

# Request All, the user starts the Inspection from the menu
# ---------------------------------------------------------

if yasara.request == "All":

    # Turning off the console increases speed and keeps the Console clean of spam
    yasara.Console("OFF")

    # Here the buttons are created.
    # ------------------------------------------------------
    yasara.MakeImage("Image1",  width=154, height=365)
    yasara.Font("Monospaced", height=14, color="black", spacing=1.5)
    yasara.ShowButton("Prev" ,x=47,y=5, color="red", action="none")
    yasara.ShowButton("Next",x=106,y=5, color="red", action="none")
    yasara.ShowButton("Go to...",x=75,y=45, color="red", action="none")
    yasara.ShowButton("HUD", x=38, y=85, color="white", action="none")
    yasara.ShowButton("Next Subunit",x=75,y=205, color="white", action="none")
    yasara.ShowButton("Water Visible",x=75,y=125, color="white", action="none")
    yasara.ShowButton("Sidechains",x=75,y=165, color="white", action="none")
    yasara.ShowButton("Status",x=101,y=85, color="white",action="none")
    yasara.ShowButton("Fog ON/OFF",x=70,y=245, color="white",action="none")
    yasara.ShowButton("Fog Fur",x=115,y=285, color="white", action="none")
    yasara.ShowButton("Fog Clo",x=38,y=285, color="white", action="none")
    yasara.ShowButton("Raytrace",x=65,y=325, color="white", action="none")
    yasara.Font("Arial", height=9, color="red", spacing=1)
    yasara.ShowButton("x",x=143,y=335, color="black", action="none")
    yasara.ShowImage("Image1",x=-10,y=-545,width=2000,height=1000,alpha=75,priority=1)

    # Initialize the Yasara view, declare visualization variables initial state
    # ------------------------------------------------------
    yasara.ColorFog("white")
    yasara.Fog(0)
    yasara.PrintHUD()
    yasara.PosText(1,20)
    HudStatus = False
    HideStatus = False
    BackboneStatus = False
    WaterVisibility = True
    FogStart = 0.5
    CurrentFogStatus = False
    yasara.ColorFog("white")
    yasara.HUD("Off")
    yasara.Fog(density="0%")

    ########################################################
    # Start with showing the first mutant, the next ~90 lines of code are the representation of one mutant and are more
    # or less redundant for the "Next" and "Prev" buttons further below
    # ------------------------------------------------------
    i = 0
    NameMutation = sorted(Fresco_yob_dict, key=fresco_grab)[i]
    MutantFile = Fresco_yob_dict[sorted(Fresco_yob_dict, key=fresco_grab)[i]]
    MutantPDB = re.sub(r"_N_\d.*",".pdb", MutantFile)
    print "MutantFile = " + MutantFile
    print "MutantPDB = " + MutantPDB

    # this is true if its a point mutation, not a disulfide bond
    if "-" not in NameMutation:
        WTFile = fresco_design_dir + Fresco_wildtype_file
        MutatedResidue = fresco_grab(NameMutation, "residue")
        MutatedSU = "A"

    # Otherwise we're dealing with a disulfide bond; in this case the WTFile could be an MD snapshot
    else:
        # check the pdb files in the Subdir_template folder containing *ps*, these are the snapshots.
        WTList = [item for item in os.listdir(fresco_unique_dsb_dir + "Subdir_templates") if re.sub(".*MD ((?:\d)*ps).*",r"\1", NameMutation) in item and item.endswith(".pdb")]
        if len(WTList) == 0:
            if "ps" in NameMutation:
                print "Warning: Couldn't find the MD snapshot file upon which this disulfide bond was predicted. Expected file in: " + os.path.realpath(fresco_namedpdbs_dir + "Subdir_template") + os.sep 
            WTFile = fresco_design_dir + Fresco_wildtype_file
        else:
            # If there is a match with a <number>ps occurance in the disulfide bond name (NameMuatation), this is the WT
            WTFile = os.path.realpath(fresco_unique_dsb_dir + "Subdir_templates") + os.sep + [item for item in os.listdir(fresco_unique_dsb_dir + "Subdir_templates") if re.sub(".*MD ((?:\d)*ps).*",r"\1", NameMutation) in item and item.endswith(".pdb")][0]
        # For disulfide bonds the two mutant residues are simply together, separated by space, as complying with Yasara syntax
        MutatedResidue = " ".join(re.findall("[A-Z](\d+)[A-Z]", NameMutation))
        # Also store the subunit in which the mutation was predicted
        MutatedSU = re.findall(", SU ([A-Z])", NameMutation)[0]

    # Clear the view and delete all Objects that were loaded
    yasara.DelObj("All")

    # this is for loading the pdb files of the WT
    # ------------------------------------------------------
    print "Loading wildtype: " + WTFile
    yasara.LoadPDB(WTFile)
    yasara.NameObj("Obj 1", "WT_static")
    # Make sure there is an MD of Wildtype
    if not os.path.isfile(fresco_namedpdbs_dir + "Subdir_template" + os.sep + Fresco_wildtype_file[:-4] + "_N_20000fs_50000fs_05_LSOn_Avg.yob"):
        yasara.RaiseError("Template MD of " + Fresco_wildtype_file + " not found. Exiting. Expected path: " + fresco_namedpdbs_dir + os.sep + "Subdir_template" + os.sep + Fresco_wildtype_file[:-4] + "_N_20000fs_50000fs_05_LSOn_Avg.yob")
        yasara.plugin.end()
    # Load all the yob files of the MD
    for x in range(1, 6):
        if not "-" in NameMutation:
            yasara.LoadYOb(fresco_namedpdbs_dir + "Subdir_template" + os.sep + Fresco_wildtype_file[:-4] + "_N_20000fs_50000fs_0" + str(x) + "_LSOn_Avg.yob")
        else:
            yasara.LoadYOb(fresco_unique_dsb_dir + "Subdir_templates" + os.sep + Fresco_wildtype_file[:-4] + "_N_20000fs_50000fs_0" + str(x) + "_LSOn_Avg.yob")
        yasara.NameObj("Obj " + str(x + 1), "WT_" + str(x) + "_MD")

    TotalWT = yasara.CountObj("All")
    yasara.ColorAtom("Obj 1 - " + str(TotalWT) + " element C", first="057000")

    # This if for loading the files of the mutant.
    # ------------------------------------------------------
    print "Loading mutant: " + MutantFile
    yasara.LoadPDB(MutantPDB)
    yasara.NameObj("Obj " + str(TotalWT + 1), "Mut_static")
    # No need to check if Yob files exist, the script wouldn't be at this point if they weren't all there
    for x in range(1, 6):
        yasara.LoadYOb(re.sub("_05_","_0" + str(x) + "_", MutantFile))
        yasara.NameObj("Obj " + str(x + TotalWT + 1), "Mutant_" + str(x) + "_MD")

    TotalMutant = str(yasara.CountObj("All"))
    yasara.ColorAtom("Obj " + str(1 + TotalWT) + " - " + str(TotalMutant) + " element C", first="0002350")

    # this is  all to improve visualization.
    # ------------------------------------------------------
    yasara.SupObj("2-12", "1", match="yes")
    yasara.Style("tube")
    yasara.ShowRes("All with distance < 8 from protein res " + str(MutatedResidue))
    yasara.ColorAtom("Obj " + str(1 + TotalWT) + " - " + str(TotalMutant) + " protein res " + str(MutatedResidue) + " element C", "magenta")
    yasara.ColorAtom("Obj " + str(1 + TotalWT) + " - " + str(TotalMutant) + " protein res " + str(MutatedResidue) + " element C", "magenta")
    yasara.HideAtom("element H with bond to atom element C")
    ListSubunits = yasara.ListMol("Obj 1 Res Protein res " + str(MutatedResidue), format="MOLNAME")
    CurrentSubunit = 1
    NumberOfSubunits = len(ListSubunits)

    yasara.CenterAtom("Obj 1 Mol " + MutatedSU + " protein res " + str(MutatedResidue if "-" not in NameMutation else MutatedResidue.split(" ")[0]) + " atom CA", coordsys="global")
    yasara.ZoomRes("Obj 1 Mol " + MutatedSU + " protein res " + str(MutatedResidue), steps=5)
    yasara.LabelPar("Arial")
    if not "-" in NameMutation:
        yasara.LabelRes("Obj " + str(1 + TotalWT) + " protein res " + str(MutatedResidue), format=NameMutation, height=0.8, color="black")
    else:
        yasara.LabelRes("Obj " + str(1 + TotalWT) + " protein res " + str(MutatedResidue.split(' ')[0]) + " Mol " + MutatedSU, format=NameMutation.split('-')[0], height=0.8, color="black")
        yasara.LabelRes("Obj " + str(1 + TotalWT) + " protein res " + str(MutatedResidue.split(' ')[1]) + " Mol " + MutatedSU, format=NameMutation.split('-')[1].split(",")[0], height=0.8, color="black")
    yasara.ShowHBoRes("visible")
    yasara.BallStickRes("!protein and !HOH")
    yasara.ColorBG("white")

    # here the multiple panels are created and the content of what is visible is set.
    # ------------------------------------------------------
    yasara.DuplicateView("Main", "WT_MD")
    yasara.DuplicateView("Main", "Mutant_MD")
    yasara.ShowView("2")
    yasara.SwitchObj("1 " + str(TotalWT + 1) + " - " + str(TotalMutant), "off")
    yasara.ShowView("3")
    yasara.SwitchObj("1 - " + str(TotalWT + 1), "off")
    yasara.ShowView("1")
    yasara.SwitchObj("2 - " + str(TotalWT), "off")
    yasara.SwitchObj(str(TotalWT + 2) + " - " + str(TotalMutant), "off")
    yasara.ShowMessage(NameMutation + " Subunit " + ListSubunits[CurrentSubunit -1])

    ########################################################
    # The definitions of what the buttons do
    # ------------------------------------------------------
    while i <= len(Fresco_yob_dict):
        button = yasara.Wait("UserButton")

        #################################################################
        # The next ~90 lines of code are nearly exactly the same as above and define the "Prev" button
        # ---------------------------------------------------------------
        if button == "Prev":
            i = max(0, i -1)
            NameMutation = sorted(Fresco_yob_dict, key=fresco_grab)[i]
            MutantFile = Fresco_yob_dict[sorted(Fresco_yob_dict, key=fresco_grab)[i]]
            MutantPDB = re.sub(r"_N_\d.*",".pdb", MutantFile)

            if not "-" in NameMutation:
                WTFile = fresco_design_dir + Fresco_wildtype_file
                MutatedResidue = fresco_grab(NameMutation,"residue")
                MutatedSU = "A"

            else:
                WTList = [item for item in os.listdir(fresco_unique_dsb_dir + "Subdir_templates") if re.sub(".*MD ((?:\d)*ps).*",r"\1", NameMutation) in item and item.endswith(".pdb")]
                if len(WTList) == 0:
                    if "ps" in NameMutation:
                        print "Warning: Couldn't find the MD snapshot file upon which this disulfide bond was predicted. Expected file in: " + os.path.realpath(fresco_namedpdbs_dir + "Subdir_template") + os.sep
                    WTFile = fresco_design_dir + Fresco_wildtype_file
                else:
                    WTFile = os.path.realpath(fresco_unique_dsb_dir + "Subdir_templates") + os.sep + [item for item in os.listdir(fresco_unique_dsb_dir + "Subdir_templates") if re.sub(".*MD ((?:\d)*ps).*",r"\1", NameMutation) in item and item.endswith(".pdb")][0]
                MutatedResidue = " ".join(re.findall("[A-Z](\d+)[A-Z]", NameMutation))
                MutatedSU = re.findall(", SU ([A-Z])", NameMutation)[0]

            yasara.DelObj("All")
            yasara.DelView("2")
            yasara.DelView("2")

            print "Loading wildtype: " + WTFile
            yasara.LoadPDB(WTFile)
            yasara.NameObj("Obj 1", "WT_static")
            if not os.path.isfile(fresco_namedpdbs_dir + "Subdir_template" + os.sep + Fresco_wildtype_file[:-4] + "_N_20000fs_50000fs_05_LSOn_Avg.yob"):
                yasara.RaiseError("Template MD of " + Fresco_wildtype_file + " not found. Exiting. Expected path: " + fresco_namedpdbs_dir + os.sep + "Subdir_template" + os.sep + Fresco_wildtype_file[:-4] + "_N_20000fs_50000fs_05_LSOn_Avg.yob")
                yasara.plugin.end()
            for x in range(1, 6):
                if not "-" in NameMutation:
                    yasara.LoadYOb(fresco_namedpdbs_dir + "Subdir_template" + os.sep + Fresco_wildtype_file[:-4] + "_N_20000fs_50000fs_0" + str(x) + "_LSOn_Avg.yob")
                else:
                    yasara.LoadYOb(fresco_unique_dsb_dir + "Subdir_templates" + os.sep + Fresco_wildtype_file[:-4] + "_N_20000fs_50000fs_0" + str(x) + "_LSOn_Avg.yob")
                yasara.NameObj("Obj " + str(x + 1), "WT_" + str(x) + "_MD")

            TotalWT = yasara.CountObj("All")
            yasara.ColorAtom("Obj 1 - " + str(TotalWT) + " element C", first="057000")

            print "Loading mutant: " + MutantFile
            yasara.LoadPDB(MutantPDB)
            yasara.NameObj("Obj " + str(TotalWT + 1), "Mut_static")

            for x in range(1, 6):
                yasara.LoadYOb(re.sub("_05_","_0" + str(x) + "_", MutantFile))
                yasara.NameObj("Obj " + str(x + TotalWT + 1), "Mutant_" + str(x) + "_MD")

            TotalMutant = str(yasara.CountObj("All"))
            yasara.ColorAtom("Obj " + str(1 + TotalWT) + " - " + str(TotalMutant) + " element C", first="0002350")

            yasara.SupObj("2-12", "1", match="yes")
            yasara.Style("tube")
            yasara.ShowRes("All with distance < 8 from protein res " + str(MutatedResidue))
            yasara.ColorAtom("Obj " + str(1 + TotalWT) + " - " + str(TotalMutant) + " protein res " + str(MutatedResidue) + " element C", "magenta")
            yasara.ColorAtom("Obj " + str(1 + TotalWT) + " - " + str(TotalMutant) + " protein res " + str(MutatedResidue) + " element C", "magenta")
            yasara.HideAtom("element H with bond to atom element C")
            ListSubunits = yasara.ListMol("Obj 1 Res Protein res " + str(MutatedResidue), format="MOLNAME")
            CurrentSubunit = 1
            NumberOfSubunits = len(ListSubunits)
            yasara.CenterAtom("Obj 1 Mol " + MutatedSU + " protein res " + str(MutatedResidue if not "-" in NameMutation else MutatedResidue.split(" ")[0]) + " atom CA", coordsys="global")
            yasara.ZoomRes("Obj 1 Mol " + MutatedSU + " protein res " + str(MutatedResidue), steps=5)
            yasara.LabelPar("Arial")
            if not "-" in NameMutation:
                yasara.LabelRes("Obj " + str(1 + TotalWT) + " protein res " + str(MutatedResidue), format=NameMutation, height=0.8, color="black")
            else:
                yasara.LabelRes("Obj " + str(1 + TotalWT) + " protein res " + str(MutatedResidue.split(' ')[0]) + " Mol " + MutatedSU, format=NameMutation.split('-')[0], height=0.8, color="black")
                yasara.LabelRes("Obj " + str(1 + TotalWT) + " protein res " + str(MutatedResidue.split(' ')[1]) + " Mol " + MutatedSU, format=NameMutation.split('-')[1].split(",")[0], height=0.8, color="black")
            yasara.ShowHBoRes("visible")
            yasara.BallStickRes("!protein and !HOH")
            yasara.ColorBG("white")

            yasara.DuplicateView("Main", "WT_MD")
            yasara.DuplicateView("Main", "Mutant_MD")
            yasara.ShowView("2")
            yasara.SwitchObj("1 " + str(TotalWT + 1) + " - " + str(TotalMutant), "off")
            yasara.ShowView("3")
            yasara.SwitchObj("1 - " + str(TotalWT + 1), "off")
            yasara.ShowView("1")
            yasara.SwitchObj("2 - " + str(TotalWT), "off")
            yasara.SwitchObj(str(TotalWT + 2) + " - " + str(TotalMutant), "off")
            yasara.ShowMessage(NameMutation + " Subunit " + ListSubunits[CurrentSubunit -1])


        #################################################################
        # The next ~90 lines of code are nearly exactly the same as above and define the "Next" button
        # ---------------------------------------------------------------
        elif button == "Next":
            i += 1
            if i == len(Fresco_yob_dict):
                yasara.ShowWin("Custom", "Finished!",580,145,
                               "Text",               80, 48, "You made it successfully through the visual inspection",
                               "Text",               20, 75, "of the single mutants and disulfide bonds suggested by FRESCO.",
                               "Button",           285, 99, "_O_K")
                i -= 1
            else:
                NameMutation = sorted(Fresco_yob_dict, key=fresco_grab)[i]
                MutantFile = Fresco_yob_dict[sorted(Fresco_yob_dict, key=fresco_grab)[i]]
                MutantPDB = re.sub(r"_N_\d.*",".pdb", MutantFile)

                if not "-" in NameMutation:
                    WTFile = fresco_design_dir + Fresco_wildtype_file
                    MutatedResidue = fresco_grab(NameMutation,"residue")
                    MutatedSU = "A"

                else:
                    WTList = [item for item in os.listdir(fresco_unique_dsb_dir + "Subdir_templates") if re.sub(".*MD ((?:\d)*ps).*",r"\1", NameMutation) in item and item.endswith(".pdb")]
                    if len(WTList) == 0:
                        if "ps" in NameMutation:
                            print "Warning: Couldn't find the MD snapshot file upon which this disulfide bond was predicted. Expected file in: " + os.path.realpath(fresco_namedpdbs_dir + "Subdir_template") + os.sep
                        WTFile = fresco_design_dir + Fresco_wildtype_file
                    else:
                        WTFile = os.path.realpath(fresco_unique_dsb_dir + "Subdir_templates") + os.sep + [item for item in os.listdir(fresco_unique_dsb_dir + "Subdir_templates") if re.sub(".*MD ((?:\d)*ps).*",r"\1", NameMutation) in item and item.endswith(".pdb")][0]
                    MutatedResidue = " ".join(re.findall("[A-Z](\d+)[A-Z]", NameMutation))
                    print MutatedResidue
                    MutatedSU = re.findall(", SU ([A-Z])", NameMutation)[0]

                yasara.DelObj("All")
                yasara.DelView("2")
                yasara.DelView("2")

                print "Loading wildtype: " + WTFile
                yasara.LoadPDB(WTFile)
                yasara.NameObj("Obj 1", "WT_static")
                if not os.path.isfile(fresco_namedpdbs_dir + "Subdir_template" + os.sep + Fresco_wildtype_file[:-4] + "_N_20000fs_50000fs_05_LSOn_Avg.yob"):
                    yasara.RaiseError("Template MD of " + Fresco_wildtype_file + " not found. Exiting. Expected path: " + fresco_namedpdbs_dir + os.sep + "Subdir_template" + os.sep + Fresco_wildtype_file[:-4] + "_N_20000fs_50000fs_05_LSOn_Avg.yob")
                    yasara.plugin.end()

                for x in range(1, 6):
                    if not "-" in NameMutation:
                        yasara.LoadYOb(fresco_namedpdbs_dir + "Subdir_template" + os.sep + Fresco_wildtype_file[:-4] + "_N_20000fs_50000fs_0" + str(x) + "_LSOn_Avg.yob")
                    else:
                        yasara.LoadYOb(fresco_unique_dsb_dir + "Subdir_templates" + os.sep + Fresco_wildtype_file[:-4] + "_N_20000fs_50000fs_0" + str(x) + "_LSOn_Avg.yob")
                    yasara.NameObj("Obj " + str(x + 1), "WT_" + str(x) + "_MD")


                TotalWT = yasara.CountObj("All")
                yasara.ColorAtom("Obj 1 - " + str(TotalWT) + " element C", first="057000")

                print "Loading mutant: " + MutantFile
                yasara.LoadPDB(MutantPDB)
                yasara.NameObj("Obj " + str(TotalWT + 1), "Mut_static")
                # No need to check if Yob files exist, the script wouldn't be at this point if they weren't all there
                for x in range(1, 6):
                    yasara.LoadYOb(re.sub("_05_","_0" + str(x) + "_", MutantFile))
                    yasara.NameObj("Obj " + str(x + TotalWT + 1), "Mutant_" + str(x) + "_MD")

                TotalMutant = str(yasara.CountObj("All"))
                yasara.ColorAtom("Obj " + str(1 + TotalWT) + " - " + str(TotalMutant) + " element C", first="0002350")

                yasara.SupObj("2-12", "1", match="yes")
                yasara.Style("tube")
                yasara.ShowRes("All with distance < 8 from protein res " + str(MutatedResidue))
                yasara.ColorAtom("Obj " + str(1 + TotalWT) + " - " + str(TotalMutant) + " protein res " + str(MutatedResidue) + " element C", "magenta")
                yasara.ColorAtom("Obj " + str(1 + TotalWT) + " - " + str(TotalMutant) + " protein res " + str(MutatedResidue) + " element C", "magenta")
                yasara.HideAtom("element H with bond to atom element C")
                ListSubunits = yasara.ListMol("Obj 1 Res Protein res " + str(MutatedResidue), format="MOLNAME")
                CurrentSubunit = 1
                NumberOfSubunits = len(ListSubunits)
                yasara.CenterAtom("Obj 1 Mol " + MutatedSU + " protein res " + str(MutatedResidue if not "-" in NameMutation else MutatedResidue.split(" ")[0]) + " atom CA", coordsys="global")
                yasara.ZoomRes("Obj 1 Mol " + MutatedSU + " protein res " + str(MutatedResidue), steps=5)
                yasara.LabelPar("Arial")
                if not "-" in NameMutation:
                    yasara.LabelRes("Obj " + str(1 + TotalWT) + " protein res " + str(MutatedResidue), format=NameMutation, height=0.8, color="black")
                else:
                    yasara.LabelRes("Obj " + str(1 + TotalWT) + " protein res " + str(MutatedResidue.split(' ')[0]) + " Mol " + MutatedSU, format=NameMutation.split('-')[0], height=0.8, color="black")
                    yasara.LabelRes("Obj " + str(1 + TotalWT) + " protein res " + str(MutatedResidue.split(' ')[1]) + " Mol " + MutatedSU, format=NameMutation.split('-')[1].split(",")[0], height=0.8, color="black")
                yasara.ShowHBoRes("visible")
                yasara.BallStickRes("!protein and !HOH")
                yasara.ColorBG("white")

                yasara.DuplicateView("Main", "WT_MD")
                yasara.DuplicateView("Main", "Mutant_MD")
                yasara.ShowView("2")
                yasara.SwitchObj("1 " + str(TotalWT + 1) + " - " + str(TotalMutant), "off")
                yasara.ShowView("3")
                yasara.SwitchObj("1 - " + str(TotalWT + 1), "off")
                yasara.ShowView("1")
                yasara.SwitchObj("2 - " + str(TotalWT), "off")
                yasara.SwitchObj(str(TotalWT + 2) + " - " + str(TotalMutant), "off")
                yasara.ShowMessage(NameMutation + " Subunit " + ListSubunits[CurrentSubunit -1])


        #################################################################
        # The next ~90 lines of code are nearly exactly the same as above and define the "Go to..." button
        # ---------------------------------------------------------------
        elif button == "Goto":
            def change(item):
                if "_" in item:
                    if item.split("_")[0].endswith("ps"):
                        return item.split("_")[2] + "-" + item.split("_")[4] + ", SU " + item.split("_")[1] + "-" + item.split("_")[3] + ", at MD " + item.split("_")[0] + ", conf: " + item.split("_")[-1]
                    else:
                        return item.split("_")[1] + "-" + item.split("_")[3] + ", SU " + item.split("_")[0] + "-" + item.split("_")[2]  + ", conf: " + item.split("_")[-1]
                else:
                    return item

            optionlist =\
              yasara.ShowWin("List","Mutant Selection List",
                             "No",
                             "Select mutant: ",
                             [item for item in sorted(Fresco_yob_dict, key=fresco_grab)])[1:]

            try:
                i = sorted(Fresco_yob_dict, key=fresco_grab).index(optionlist[0])
            except:
                print "Nothing selected! Loading previous view"

            NameMutation = sorted(Fresco_yob_dict, key=fresco_grab)[i]
            MutantFile = Fresco_yob_dict[sorted(Fresco_yob_dict, key=fresco_grab)[i]]
            MutantPDB = re.sub(r"_N_\d.*",".pdb", MutantFile)

            if not "-" in NameMutation:
                WTFile = fresco_design_dir + Fresco_wildtype_file
                MutatedResidue = fresco_grab(NameMutation,"residue")
                MutatedSU = "A"
            else:
                WTList = [item for item in os.listdir(fresco_unique_dsb_dir + "Subdir_templates") if re.sub(".*MD ((?:\d)*ps).*",r"\1", NameMutation) in item and item.endswith(".pdb")]
                if len(WTList) == 0:
                    if "ps" in NameMutation:
                        print "Warning: Couldn't find the MD snapshot file upon which this disulfide bond was predicted. Expected file in: " + os.path.realpath(fresco_namedpdbs_dir + "Subdir_template") + os.sep
                    WTFile = fresco_design_dir + Fresco_wildtype_file
                else:
                    WTFile = os.path.realpath(fresco_unique_dsb_dir + "Subdir_templates") + os.sep + [item for item in os.listdir(fresco_unique_dsb_dir + "Subdir_templates") if re.sub(".*MD ((?:\d)*ps).*",r"\1", NameMutation) in item and item.endswith(".pdb")][0]
                MutatedResidue = " ".join(re.findall("[A-Z](\d+)[A-Z]", NameMutation))
                MutatedSU = re.findall(", SU ([A-Z])", NameMutation)[0]

            yasara.DelObj("All")
            yasara.DelView("2")
            yasara.DelView("2")

            print "Loading wildtype: " + WTFile
            yasara.LoadPDB(WTFile)
            yasara.NameObj("Obj 1", "WT_static")
            if not os.path.isfile(fresco_namedpdbs_dir + "Subdir_template" + os.sep + Fresco_wildtype_file[:-4] + "_N_20000fs_50000fs_05_LSOn_Avg.yob"):
                yasara.RaiseError("Template MD of " + Fresco_wildtype_file + " not found. Exiting. Expected path: " + fresco_namedpdbs_dir + os.sep + "Subdir_template" + os.sep + Fresco_wildtype_file[:-4] + "_N_20000fs_50000fs_05_LSOn_Avg.yob")
                yasara.plugin.end()

            for x in range(1, 6):
                if not "-" in NameMutation:
                    yasara.LoadYOb(fresco_namedpdbs_dir + "Subdir_template" + os.sep + Fresco_wildtype_file[:-4] + "_N_20000fs_50000fs_0" + str(x) + "_LSOn_Avg.yob")
                else:
                    yasara.LoadYOb(fresco_unique_dsb_dir + "Subdir_templates" + os.sep + Fresco_wildtype_file[:-4] + "_N_20000fs_50000fs_0" + str(x) + "_LSOn_Avg.yob")
                yasara.NameObj("Obj " + str(x + 1), "WT_" + str(x) + "_MD")

            TotalWT = yasara.CountObj("All")
            yasara.ColorAtom("Obj 1 - " + str(TotalWT) + " element C", first="057000")

            print "Loading mutant: " + MutantFile
            yasara.LoadPDB(MutantPDB)
            yasara.NameObj("Obj " + str(TotalWT + 1), "Mut_static")
            # No need to check if Yob files exist, the script wouldn't be at this point if they weren't all there
            for x in range(1, 6):
                yasara.LoadYOb(re.sub("_05_","_0" + str(x) + "_", MutantFile))
                yasara.NameObj("Obj " + str(x + TotalWT + 1), "Mutant_" + str(x) + "_MD")

            TotalMutant = str(yasara.CountObj("All"))
            yasara.ColorAtom("Obj " + str(1 + TotalWT) + " - " + str(TotalMutant) + " element C", first="0002350")

            yasara.SupObj("2-12", "1", match="yes")
            yasara.Style("tube")
            yasara.ShowRes("All with distance < 8 from protein res " + str(MutatedResidue))
            yasara.ColorAtom("Obj " + str(1 + TotalWT) + " - " + str(TotalMutant) + " protein res " + str(MutatedResidue) + " element C", "magenta")
            yasara.ColorAtom("Obj " + str(1 + TotalWT) + " - " + str(TotalMutant) + " protein res " + str(MutatedResidue) + " element C", "magenta")
            yasara.HideAtom("element H with bond to atom element C")
            ListSubunits = yasara.ListMol("Obj 1 Res Protein res " + str(MutatedResidue), format="MOLNAME")
            CurrentSubunit = 1
            NumberOfSubunits = len(ListSubunits)
            yasara.CenterAtom("Obj 1 Mol " + MutatedSU + " protein res " + str(MutatedResidue if not "-" in NameMutation else MutatedResidue.split(" ")[0]) + " atom CA", coordsys="global")
            yasara.ZoomRes("Obj 1 Mol " + MutatedSU + " protein res " + str(MutatedResidue), steps=5)
            yasara.LabelPar("Arial")
            if not "-" in NameMutation:
                yasara.LabelRes("Obj " + str(1 + TotalWT) + " protein res " + str(MutatedResidue), format=NameMutation, height=0.8, color="black")
            else:
                yasara.LabelRes("Obj " + str(1 + TotalWT) + " protein res " + str(MutatedResidue.split(' ')[0]) + " Mol " + MutatedSU, format=NameMutation.split('-')[0], height=0.8, color="black")
                yasara.LabelRes("Obj " + str(1 + TotalWT) + " protein res " + str(MutatedResidue.split(' ')[1]) + " Mol " + MutatedSU, format=NameMutation.split('-')[1].split(",")[0], height=0.8, color="black")
            yasara.ShowHBoRes("visible")
            yasara.BallStickRes("!protein and !HOH")
            yasara.ColorBG("white")

            yasara.DuplicateView("Main", "WT_MD")
            yasara.DuplicateView("Main", "Mutant_MD")
            yasara.ShowView("2")
            yasara.SwitchObj("1 " + str(TotalWT + 1) + " - " + str(TotalMutant), "off")
            yasara.ShowView("3")
            yasara.SwitchObj("1 - " + str(TotalWT + 1), "off")
            yasara.ShowView("1")
            yasara.SwitchObj("2 - " + str(TotalWT), "off")
            yasara.SwitchObj(str(TotalWT + 2) + " - " + str(TotalMutant), "off")
            yasara.ShowMessage(NameMutation + " Subunit " + ListSubunits[CurrentSubunit -1])


        #################################################################
        # The other buttons
        # ---------------------------------------------------------------
        elif button == "HUD":
            if HudStatus is True:
                yasara.HUD(show="Off")
                HudStatus = False
            else:
                yasara.HUD(show="Obj")
                HudStatus = True
        elif button == "NextSubunit":
            if NumberOfSubunits == 1:
                yasara.ShowMessage("There is only 1 subunit detected")
                yasara.Wait(50)
                yasara.ShowMessage(NameMutation + " Subunit " + ListSubunits[CurrentSubunit -1])
            else:
                CurrentSubunit += 1
                if CurrentSubunit > NumberOfSubunits:
                    CurrentSubunit = 1
                yasara.CenterAtom("OBJ 1 MOL " + ListSubunits[CurrentSubunit -1] + " protein res " + str(MutatedResidue if not "-" in NameMutation else MutatedResidue.split(" ")[0]) + " atom CA", coordsys="global")
                yasara.ZoomRes("OBJ 1 MOL " + ListSubunits[CurrentSubunit -1] + " protein res " + str(MutatedResidue), steps=5)

                # do some reporting
                yasara.ShowMessage("Zoomed to subunit " + ListSubunits[CurrentSubunit -1])
                yasara.Wait(50)
                yasara.ShowMessage(NameMutation + " Subunit " + ListSubunits[CurrentSubunit -1])

        elif button == "WaterVisible":
            if WaterVisibility == True:
              yasara.HideRes("CIP CIM HOH")
              WaterVisibility = False
            else:
              yasara.ShowRes("CIP CIM HOH with distance < 8 from protein res " + str(MutatedResidue))
              yasara.ShowHBoRes("visible")
              WaterVisibility = True

        elif button == "Sidechains":
            if BackboneStatus == 0:
              yasara.HideRes("protein res")
              yasara.HideHBoRes("protein res")
              BackboneStatus = True
            else:
              yasara.ShowRes("All with distance < 8 from protein res " + str(MutatedResidue))
              yasara.ShowRes("Res " + str(MutatedResidue))
              yasara.ShowHBoRes("visible")
              yasara.HideAtom("element H with bond to atom element C")
              BackboneStatus = False

        elif button == "Status":
            if HideStatus == True:
                yasara.ShowMessage(NameMutation + " Subunit " + ListSubunits[CurrentSubunit -1])
                HideStatus = False
            else:
                yasara.HideMessage()
                HideStatus = True

        elif button == "FogONOFF":
            if CurrentFogStatus is True:
                yasara.ColorBG("white")
                yasara.Fog(density="0%")
                CurrentFogStatus = False
            else:
                yasara.ColorBG("white")
                yasara.ColorFog("white")
                yasara.Fog(density="100%",Range="Static",dismin=str(0.5) + "%", dismax= str(1) + "%")
                CurrentFogStatus = True

        elif button == "FogFur":
            if FogStart == 20:
                yasara.ShowMessage("Already maximum fog distance")
                yasara.Wait(50)
                yasara.ShowMessage(NameMutation + " Subunit " + ListSubunits[CurrentSubunit -1])
            else:
                FogStart = FogStart + 0.25
                yasara.Fog(density="100%",Range="Static",dismin=str(FogStart) + "%",dismax=(str(FogStart + 0.5) + "%"))
                yasara.ShowMessage(NameMutation + " Subunit " + ListSubunits[CurrentSubunit -1])

        elif button == "FogClo":
            if FogStart == 0:
                yasara.ShowMessage("Already minimum fog distance")
                yasara.Wait(50)
                yasara.ShowMessage(NameMutation + " Subunit " + ListSubunits[CurrentSubunit -1])
            else:
                FogStart = FogStart - 0.25
                yasara.Fog(density="100%",Range="Static",dismin=str(FogStart) +"%",dismax=(str(FogStart + 0.5) +"%"))
                yasara.ShowMessage(NameMutation + " Subunit " + ListSubunits[CurrentSubunit -1])

        elif button == "Raytrace":
            # determine the screensize for later.
            try:
                w,h = yasara.ScreenSize()
            except:
                w,h,n = yasara.ScreenSize()
            if not os.path.exists(fresco_design_dir + os.sep + "pictures"):
                os.makedirs(fresco_design_dir + os.sep + "pictures")
            FileAlreadyExists = True
            i = 0000
            while FileAlreadyExists:
              i += 1
              FileAlreadyExists = yasara.FileSize("pictures/" + NameMutation + "_Subunit_" + ListSubunits[CurrentSubunit -1] + "_" + str(i) + ".png")
            yasara.ShowMessage("A ray-traced screenshot of " + str(w) + " by " + str(h) + "pixels will be saved as pictures/" + NameMutation + "_Subunit_" + ListSubunits[CurrentSubunit -1] + "_" + str(i) + ".png")
            yasara.RayTrace("pictures/" + NameMutation + "_Subunit_" + ListSubunits[CurrentSubunit -1] + "_" + str(i) + ".png",x=w,y=h,zoom=1,atoms="Balls",labelshadow="No",secalpha=100,display="No",outline=0,background="Yes")
            yasara.ShowMessage(NameMutation + " Subunit " + ListSubunits[CurrentSubunit -1])
        elif button == "x":
            yasara.DelImage("All")
            yasara.DelView("2")
            yasara.DelView("2")
            yasara.HideMessage()
            yasara.UnlabelRes("protein")
            yasara.Fog(density="0%")
            yasara.plugin.end()
    yasara.Console("ON")

#################################################################

#################################################################
# Preparation of the Excel File
# ------------------------------------------------------

elif yasara.request == "Excel":
    Fresco_wildtype = Fresco_wildtype_file[:-12]
    if len(yasara.ListObj("all")) > 1:
        yasara.RaiseError("You have more than one object loaded. Try again, when only one object - that is, the wildtype file - is loaded.")
        yasara.plugin.end()
    yasara.ShowWin("Custom", "Warning",580,145,
                   "Text",               80, 48, "If the list of mutations is large, this can take a while.",
                   "Text",               20, 75, "Yasara may appear to be frozen. Wait until you see a window pop up. ",
                   "Button",           285, 99, "_O_K")
    excel_out = open(fresco_design_dir + os.sep + Fresco_wildtype + "MutantsForInspections_XL.tab", "w+")
    excel_out.write("#\tMD available\tMut\tddG (kJ)\tStd Dev (kJ)\tHydrophob +/-\tKept or dicarded for reason\n")
    excel_out_verbose = open(fresco_design_dir + os.sep + Fresco_wildtype + "MutantsForInspections_verbose.tab", "w+")
    excel_out_verbose.write("#\tWT Res\tMut Res\tMD availablet\tMut\tddG (kJ)\tStd Dev (kJ)\tWT Hbonds\tMut Hbonds\tHBonds +/-\tWT Surf Exp\tMut Surf Exp\tSurf Exp +/-\tWT Hydrophob\tMut Hydrophob\tHydrophob +/-\tKept or dicarded for reason\n")
    excel_in = open(yasara.selection[0].filename[0], "r")
    next(excel_in)
    yasara.Console("OFF")
    # Go through all the lines in the input .tab file
    for line in excel_in:
        num = re.search("^(([A-Z])(\d*)([A-Z]))",line)

        # try to find the corresponding files and folders
        fileloc = fresco_design_dir + "NamedPdbFiles" + os.sep + "Subdir_" + num.group(1) + os.sep + Fresco_wildtype + "_cleaned_" + num.group(1) + "_WOW.pdb"
        folderloc = fresco_design_dir + "NamedPdbFiles" + os.sep + "Subdir_" + num.group(1)

        if os.path.exists(folderloc) and len([item for item in os.listdir(folderloc) if ".yob" in item]) == 5:
            Mdexists = "Yes"
        else:
            Mdexists = "No"

        if os.path.isfile(os.path.realpath(fileloc)):
            # calculate the reported values
            objectnr = yasara.LoadPDB(fileloc)[0]
            yasara.AddHydAll("All")
            yasara.ShowHBoRes("All") #Show all the H Bonds
            yasara.AddEnvAll()
            hbonds_mut = len(yasara.ListHBoAtom("Obj 2 Res " + num.group(3),"Obj 2",Min=6.25,results=1))
            hbonds_wt = len(yasara.ListHBoAtom("Obj 1 Res " + num.group(3),"Obj 1",Min=6.25,results=1))
            hbonds = str("{0:+.03f}".format(hbonds_mut - hbonds_wt))
            surfexp_mut = yasara.SurfAtom("Obj 2 Res " + num.group(3), Type="accessible")[0]
            surfexp_wt = yasara.SurfAtom("Obj 1 Res " + num.group(3), Type="accessible")[0]
            surfexp = str("{0:+.03f}".format(round(surfexp_mut - surfexp_wt,2)))
            hydrophob_mut = float(Amino_acid_dict[num.group(4)][1])
            hydrophob_wt = float(Amino_acid_dict[num.group(2)][1])
            hydrophob = str("{0:+.03f}".format(hydrophob_mut - hydrophob_wt))

            out_string = num.group(3) + "\t" + Mdexists + "\t" + re.sub("\s+","\t",line) + hydrophob + "\n"
            excel_out.write(out_string)

            out_string_verbose = num.group(3) + "\t" + Amino_acid_dict[num.group(2)][0] + "\t" + Amino_acid_dict[num.group(4)][0] + "\t" + Mdexists + "\t" + re.sub("\s+","	",line) + str(hbonds_wt) + "\t" + str(hbonds_mut) + "\t" + hbonds + "\t" + str(surfexp_wt) + "\t" + str(surfexp_mut) + "\t" + surfexp + "\t" + str(hydrophob_wt) + "\t" + str(hydrophob_mut) + "\t" + hydrophob + "\n"
            excel_out_verbose.write(out_string_verbose)

            yasara.DelObj("Obj " + str(objectnr))

        else:
            yasara.write("Error! Couldn't find file: " + fileloc)
            file_error = True
    excel_out.close()
    excel_out_verbose.close()
    if "file_error" in locals():
        yasara.ShowWin("Custom", "Job done!",542,225,
               "Text",               162, 48, "Succesfully created file:",
               "Text",               20, 85, "../fresco" + Fresco_wildtype + "/designsMD/" + Fresco_wildtype + "MutantsForInspections_XL.tab",
               "Text",              20, 125, "Some mutants specified in the input list couldn't be found in the",
               "Text",              20, 155, "corresponding folders, however. Check the console for details.",
               "Button",           267, 180, "_O_K")
    else:
        yasara.ShowWin("Custom", "Job done!",580,150,
               "Text",               180, 48, "Succesfully created file:",
               "Text",               20, 75, "../designsMD/" + Fresco_wildtype + "MutantsForInspections_XL.tab",
               "Button",           285, 100, "_O_K")
    os.system(("start " if os.sep is "\\" else "open ") + fresco_design_dir + os.sep + Fresco_wildtype + "MutantsForInspections_XL.tab")
    yasara.Console("ON")

yasara.plugin.end()
