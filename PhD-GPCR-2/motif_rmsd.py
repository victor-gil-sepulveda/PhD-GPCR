'''
Created on Dec 21, 2015

@author: victor
'''
from optparse import OptionParser
import seaborn as sns
from histogram import parse_motifs
import os
from prody.proteins.pdbfile import parsePDB
from pyRMSD.RMSDCalculator import RMSDCalculator
from prepare_clusters_for_visualization import create_directory
import cPickle as pickle
import matplotlib.pyplot as plt
import pandas as pd

if __name__ == '__main__':
    sns.set_style("whitegrid")
    
    parser = OptionParser()
    parser.add_option("--motifs", dest="motifs")
    parser.add_option("--results", dest="results")
    parser.add_option("-d", dest="data")
    
    
    (options, args) = parser.parse_args()
    
    if not options.results:  
        parser.error('Results folder not given')
    else:
        create_directory(options.results)
        
    if options.motifs is None:
        parser.error('You must specify the motifs definition file.')
    else:
        motifs = parse_motifs(options.motifs)
        ordered_motifs = ['F/I-II', 'A/IV', 'B/V', 'C/VI', 'D', 'E/VII', 'Priming Loop']
    
    if options.data is None:
        rmsd_results = {"Drug":[], "Protein":[],"RMSD":[],"Motif":[]} 
        for drug in ["CMA","CMC","DMA"]:
            for protein in ["BVDV", "HCV", "JEV", "Polio", "TBEV" ]:
                path = os.path.join(drug, "ca_%s_%s.pdb"%(protein,drug))
                print "Working with", path
                pdb = parsePDB(path, subset = "ca", csets=range(1000))
                print "Loaded"
                for motif in ordered_motifs:
                    if motif in motifs[protein]:
                        cas = pdb.select("name CA")
                        motif_cas = pdb.select("resid %d to %d"%(motifs[protein][motif][0],motifs[protein][motif][1]))
                        calculator = RMSDCalculator(calculatorType = "QCP_OMP_CALCULATOR",
                                                    fittingCoordsets = cas.getCoordsets(),
                                                    calculationCoordsets = motif_cas.getCoordsets())
                        rmsds = calculator.oneVsFollowing(0)
                        rmsd_results["RMSD"].extend(rmsds)
                        rmsd_results["Protein"].extend([protein]*len(rmsds))
                        rmsd_results["Drug"].extend([drug]*len(rmsds))
                        rmsd_results["Motif"].extend([motif]*len(rmsds))
                        del cas
                        del motif_cas
                del pdb
                    
        pickle.dump(rmsd_results, open(os.path.join(options.results,"rmsd_results.pickle"),"w"))
    else:
        rmsd_results = pickle.load(open(options.data))
    
    db = pd.DataFrame.from_dict(rmsd_results)
    sns.FacetGrid(data = db, row = "Drug", column = "Protein", hue = "Motif")
    sns.barplot("Protein", "RMSD", hue = "Motif")
    plt.show()
            
                