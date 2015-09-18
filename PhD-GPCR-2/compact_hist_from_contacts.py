"""
Created on 14/9/2015

@author: victor
"""

import sys
import histogram
import tools
from _collections import defaultdict
import matplotlib.pyplot as plt
import numpy
import os

if __name__ == '__main__':
    contacts_file_with_paths = sys.argv[1]
    results_dir = sys.argv[2]
    tools.create_dir(results_dir)
    
    handler = open(contacts_file_with_paths,"r")
    
    data = defaultdict(list)
    for line in handler:
        drug, protein, contacts_file_path = line.split()
        _, contacts_per_residue = histogram.parse_contacts(contacts_file_path)
        data[protein].append((drug, contacts_per_residue))
    handler.close()
    
    colors = tools.define_tableau20_cm()
    for protein in data:
        fig, ax = plt.subplots()
        ax.set_title(protein)
        
        THRESHOLD = 1000
        # filter residues, leave only those that have at least one value higher than a threshold
        all_residue_labels = []
        for drug, values in data[protein]:
#             all_residue_labels.extend( values.keys())
            for residue in values:
                if values[residue]>=THRESHOLD:
                    all_residue_labels.append(residue)
        all_residue_labels = sorted(list(set(all_residue_labels)))
            
        
        n_groups = len(all_residue_labels)
        bar_width = 0.30
        index = numpy.arange(n_groups) + bar_width/2
        legend_label = []
        for i, (drug, values) in enumerate(data[protein]):
            x = numpy.array([values[res] for res in all_residue_labels])
            print x
            rects = ax.bar(index, x, bar_width, color = colors[i])
            index = index + bar_width
            legend_label.append(drug)
        
        plt.subplots_adjust(bottom=0.16) 
        #ax.autoscale(tight=True)
        ax.legend( legend_label)
        ax.set_xticks(numpy.arange(n_groups)+bar_width/2)
        ax.set_xticklabels(all_residue_labels, rotation= 45)
        plt.savefig("%s.svg"%(os.path.join(results_dir, protein)))
        #plt.show()