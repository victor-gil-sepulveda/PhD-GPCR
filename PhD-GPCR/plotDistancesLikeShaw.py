"""
Created on 13/03/2014

@author: victor
"""


import prody
import os.path
from pyRMSD.RMSDCalculator import RMSDCalculator
import matplotlib
import pylab
import os
from tools import distance, create_dir


data = [
    {
       'dir':'2rh1_refined',
       'pdb_traj':'2rh1_refined.pdb',
        'selection':{
            'backbone':"backbone and chain A",
            'motif_all': "noh resnum 322:327", # Todos heavy
            'motif_backbone': "backbone resnum 322:327", # Solo backbone
            'arg131_leu272': "ca resnum 131 272", # No hace falta alignment
            'ser207':"backbone resnum 207", # Alignment backbone / RMSD backbone
            'ser207_OG':"name OG resnum 207",
            'lig_OAE':"name OAE", # Unique name in the ligand
            'ile121': "backbone resnum 121",
            'phe282': "backbone resnum 282"
        }
    }
]

cwd = os.getcwd()
for datum in data:
    prot_name = datum['dir']
    print "========================\nWorking with %s\n========================"%(prot_name)
    # Look for the directory and enter it
    base_dir = os.path.join(cwd, prot_name)
    os.chdir(base_dir)
    # Create dirs to put results
    create_dir("plots")
    create_dir("selections")

    # Generate all selections
    pdb = prody.parsePDB("%s"%datum['pdb_traj'])
    selections = {}
    for selection in datum['selection']:
        print selection
        selections[selection] = pdb.select(datum['selection'][selection])
        prody.writePDB(os.path.join("selections",selection), selections[selection])

    #############################
    # Motif VS Helix Distance
    #############################
    calculator = RMSDCalculator( calculatorType = "QCP_OMP_CALCULATOR",
                 fittingCoordsets = selections['motif_all'].getCoordsets(),
                 calculationCoordsets = selections['motif_backbone'].getCoordsets())

    motif_rmsd = calculator.oneVsTheOthers( conformation_number = 0,
                                            get_superposed_coordinates = False)

    residue_distances = []
    for conf in selections['arg131_leu272'].getCoordsets():
        arg131 = conf[0]
        leu272 = conf[1]
        residue_distances.append(distance(arg131, leu272))
    exp_motif_rmsd = [0]+list(motif_rmsd)
    matplotlib.pyplot.scatter(residue_distances, exp_motif_rmsd)
    matplotlib.pyplot.savefig(os.path.join("plots", "motif_vs_helix_dist.svg"))
    matplotlib.pyplot.close()

    ###########################################
    # Backbone of ser 207 rmsd Vs ser to ligand distance
    ###########################################

    calculator = RMSDCalculator( calculatorType = "QCP_OMP_CALCULATOR",
                 fittingCoordsets = selections['backbone'].getCoordsets(),
                 calculationCoordsets = selections['ser207'].getCoordsets())

    ser207_rmsd = calculator.oneVsTheOthers( conformation_number = 0,
                                             get_superposed_coordinates = False)

    # Distance from ligand to serine
    distances = []
    lig_coords = selections["lig_OAE"].getCoordsets()
    ser_coords = selections["ser207_OG"].getCoordsets()
    for i, ser_OG in enumerate(ser_coords):
        lig_OAE = lig_coords[i]
        distances.append(distance(ser_OG[0], lig_OAE[0]))
    exp_ser_207_rmsd = [0]+list(ser207_rmsd)
    matplotlib.pyplot.scatter(distances, exp_ser_207_rmsd)
    matplotlib.pyplot.savefig(os.path.join("plots", "ser207_dist_lig.svg"))
    matplotlib.pyplot.close()

    ###########################################
    # Ile 121 rmsd vs Phe 282 rmsd
    ###########################################
    ile121_coordsets = pdb.select(datum['selection']['ile121'])
    phe282_coordsets = pdb.select(datum['selection']['phe282'])

    calculator = RMSDCalculator( calculatorType = "QCP_OMP_CALCULATOR",
                 fittingCoordsets = selections['backbone'].getCoordsets(),
                 calculationCoordsets = selections['ile121'].getCoordsets())

    ile121_rmsd = calculator.oneVsTheOthers( conformation_number = 0,
                                             get_superposed_coordinates = False)

    calculator = RMSDCalculator( calculatorType = "QCP_OMP_CALCULATOR",
                 fittingCoordsets = selections['backbone'].getCoordsets(),
                 calculationCoordsets = selections['phe282'].getCoordsets())

    phe282_rmsd = calculator.oneVsTheOthers( conformation_number = 0,
                                         get_superposed_coordinates = False)

    matplotlib.pyplot.scatter(phe282_rmsd, ile121_rmsd)
    matplotlib.pyplot.savefig(os.path.join("plots", "ile121_vs_phe282.svg"))
    matplotlib.pyplot.close()
