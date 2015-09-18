'''
Created on 18/9/2015

@author: victor
'''

import sys
import histogram
import tools
from _collections import defaultdict
import matplotlib.pyplot as plt
import numpy
import os
from histogram import get_labels, filter_less_contacts_than,\
    plot_histogram, get_num_contacts_per_residue

if __name__ == '__main__':
    contacts_file_with_paths = sys.argv[1]
    results_dir = sys.argv[2]
    tools.create_dir(results_dir)
    
    handler = open(contacts_file_with_paths,"r")
    
    data = []
    for line in handler:
        drug, protein, contacts_file_path = line.split()
        contacts_per_cluster, contacts_per_residue = histogram.parse_contacts(contacts_file_path)
        data.append((protein, drug, contacts_per_residue, contacts_per_cluster))
    handler.close()
    
    colors = tools.define_tableau20_cm()

    for protein, drug, contacts_per_residue, contacts_per_cluster in data:
        
        contact_residue_labels =  get_labels (contacts_per_cluster)
        
        # A normal plot
        target = os.path.join(results_dir, "%s_%s_histogram.svg"%(drug,protein))
        plot_histogram(contacts_per_cluster, contact_residue_labels, target, False)
        
        # A plot averaging
        target = os.path.join(results_dir, "%s_%s_histogram_a.svg"%(drug,protein))
        plot_histogram(contacts_per_cluster, contact_residue_labels, target, True)
    
        # A plot filtering
        filtered_contact_residue_labels =  filter_less_contacts_than(1000, contact_residue_labels, contacts_per_residue)
        target = os.path.join(results_dir, "%s_%s_histogram_f_1k.svg"%(drug,protein))
        plot_histogram(contacts_per_cluster, filtered_contact_residue_labels, target, False)
    
        # A plot filtering + averaging
        filtered_contact_residue_labels =  filter_less_contacts_than(1000, contact_residue_labels, contacts_per_residue)
        target = os.path.join(results_dir, "%s_%s_histogram_fa_1k.svg"%(drug,protein))
    
        # A plot filtering
        filtered_contact_residue_labels =  filter_less_contacts_than(2000, contact_residue_labels, contacts_per_residue)
        target = os.path.join(results_dir, "%s_%s_histogram_f_2k.svg"%(drug,protein))
        plot_histogram(contacts_per_cluster, filtered_contact_residue_labels, target, False)
    
        # A plot filtering + averaging
        filtered_contact_residue_labels =  filter_less_contacts_than(2000, contact_residue_labels, contacts_per_residue)
        target = os.path.join(results_dir, "%s_%s_histogram_fa_2k.svg"%(drug,protein))