'''
Created on 18/9/2015

@author: victor
'''

import sys
import histogram
import tools
import os
from collections import defaultdict
from histogram import get_labels, filter_less_contacts_than,\
    plot_histogram, parse_drug_info, parse_frames_info, get_max_val

if __name__ == '__main__':
    contacts_file_with_paths = sys.argv[1]
    results_dir = sys.argv[2]
    tools.create_dir(results_dir)
    
    num_atoms_per_drug = parse_drug_info(sys.argv[3])
    frames_per_prot_drug = parse_frames_info(sys.argv[4])
    
    handler = open(contacts_file_with_paths,"r")
    data = []
    for line in handler:
        drug, protein, contacts_file_path = line.split()
        contacts_per_cluster, contacts_per_residue = histogram.parse_contacts(contacts_file_path)
        data.append((protein, drug, contacts_per_residue, contacts_per_cluster))
    handler.close()
    
    colors = tools.define_tableau20_cm()
    
    max_val = 0
    for protein, drug, contacts_per_residue, contacts_per_cluster in data:
        contact_residue_labels =  get_labels (contacts_per_cluster)
        weight = num_atoms_per_drug[drug]*frames_per_prot_drug[protein][drug]
        max_val = max(max_val, get_max_val(contacts_per_cluster, contact_residue_labels, weight))

    for protein, drug, contacts_per_residue, contacts_per_cluster in data:
        contact_residue_labels =  get_labels (contacts_per_cluster)
        weight = num_atoms_per_drug[drug]*frames_per_prot_drug[protein][drug]
        print protein, drug, contacts_per_residue
        
        total = 0
        for res_label in contacts_per_residue:
            total+=contacts_per_residue[res_label]
        print total, total/weight
        for cluster in  contacts_per_cluster:
            print cluster, len(contacts_per_cluster[cluster])
        print "******"
        for residue in contacts_per_residue:
            print residue, contacts_per_residue[residue]
        print "******"
#         # A plot filtering
        filtered_contact_residue_labels =  filter_less_contacts_than(1000, contact_residue_labels, contacts_per_residue)
        total_filtered = 0
        for residue in filtered_contact_residue_labels:
            print residue, contacts_per_residue[residue]
            total_filtered += contacts_per_residue[residue]
        print total_filtered, float(total_filtered)/total
        exit()
        target = os.path.join(results_dir, "%s_%s_histogram_f_1k.svg"%(drug,protein))
        plot_histogram(protein, drug, max_val, contacts_per_cluster, filtered_contact_residue_labels, target, weight, False)
         

#         # A normal plot
#         target = os.path.join(results_dir, "%s_%s_histogram.svg"%(drug,protein))
#         plot_histogram(contacts_per_cluster, contact_residue_labels, target, weight, False)
#         
#         # A plot averaging
#         target = os.path.join(results_dir, "%s_%s_histogram_a.svg"%(drug,protein))
#         plot_histogram(contacts_per_cluster, contact_residue_labels, target, weight, True)
    
#         # A plot filtering + averaging
#         filtered_contact_residue_labels =  filter_less_contacts_than(1000, contact_residue_labels, contacts_per_residue)
#         target = os.path.join(results_dir, "%s_%s_histogram_fa_1k.svg"%(drug,protein))
#         plot_histogram(contacts_per_cluster, filtered_contact_residue_labels, target, weight, True)
    
#         # A plot filtering
#         filtered_contact_residue_labels =  filter_less_contacts_than(2000, contact_residue_labels,  contacts_per_residue)
#         target = os.path.join(results_dir, "%s_%s_histogram_f_2k.svg"%(drug,protein))
#         plot_histogram(contacts_per_cluster, filtered_contact_residue_labels, target, weight, False)
#     
#         # A plot filtering + averaging
#         filtered_contact_residue_labels =  filter_less_contacts_than(2000, contact_residue_labels,  contacts_per_residue)
#         target = os.path.join(results_dir, "%s_%s_histogram_fa_2k.svg"%(drug,protein))
#         plot_histogram(contacts_per_cluster, filtered_contact_residue_labels, target, weight, True)
    