'''
Created on 18/9/2015

@author: victor
'''

import histogram
import tools
import os
from collections import defaultdict
from histogram import get_labels, filter_less_contacts_than,\
    plot_histogram, parse_drug_info, parse_frames_info, get_max_val,\
    plot_histogram_to_axis, order_labels, calc_contacts_per_motif
from optparse import OptionParser
import seaborn as sns
import matplotlib.pyplot as plt
import numpy
from prody import proteins

if __name__ == '__main__':
    sns.set_style("whitegrid")
    
    parser = OptionParser()
    parser.add_option("--contacts", dest="contacts")
    parser.add_option("--motifs", dest="motifs")
    parser.add_option("--results", dest="results")
    parser.add_option("--drug-atoms", dest="drug_atoms")
    parser.add_option("--frames-per-drug", dest="frames_per_drug")
    
    (options, args) = parser.parse_args()
    
    if options.contacts is None:
        parser.error('Contacts file is needed."')
    
    if options.motifs is None:
        parser.error('You must specify the motifs definition file.')
    
    if options.results is None:
        parser.error('You must specify the results folder')
        
    if options.drug_atoms is None:
        parser.error('You must specify the drug_atoms specification file.')
                     
    if options.frames_per_drug is None:
        parser.error('You must specify the frames_per_drug specification file.')

    tools.create_dir(options.results)
    
    THRESHOLD = 1000
    num_atoms_per_drug = parse_drug_info(options.drug_atoms)
    frames_per_prot_drug = parse_frames_info(options.frames_per_drug)
    
    # Parse contact files
    handler = open(options.contacts,"r")
    proteins = []
    drugs = []

    data = {
            "contacts_per_cluster":defaultdict(lambda:defaultdict()),
            "contacts_per_residue":defaultdict(lambda:defaultdict())
            }
    for line in handler:
        if line[0]!= "#":
            drug, protein, contacts_file_path = line.split()
            if not protein in proteins: proteins.append(protein)
            if not drug in drugs: drugs.append(drug)
            contacts_per_cluster, contacts_per_residue = histogram.parse_contacts(contacts_file_path)
            data["contacts_per_cluster"][protein][drug] =  contacts_per_cluster
            data["contacts_per_residue"][protein][drug] =  contacts_per_residue
    handler.close()
    
    # Parse motifs
    motifs = defaultdict(dict)
    for line in open(options.motifs):
        if line[0]!= "#":
            parts = [ s.strip() for s in line.split(",")]
            sel_parts = parts[2].split()
            motifs[parts[0]][parts[1]] = (int(sel_parts[1]),int(sel_parts[3]))
    
    # overide motif order (can be checked with commented code below)
    ordered_motifs = ['F/I-II', 'A/IV', 'B/V', 'C/VI', 'D', 'E/VII', 'Priming Loop']
    # Sort drug names
    drugs.sort() 
    # philogenic ordering
    proteins = ["JEV", "WNV", "TBEV", "BVDV", "HCV", "Polio"]
#     # Order motifs per protein
#     ordered_motifs = {} 
#     for protein in motifs:
#         ordered_motifs[protein] = [c for s,c in sorted([(motifs[protein][k],k) for k in motifs[protein]])]
#         print protein, ordered_motifs[protein]
#     
    
    
    # PLOT CONTACTS PER RESIDUE
    f, axes = plt.subplots(len(proteins), len(drugs), sharey='row')
    f.subplots_adjust(hspace=0.3, wspace=0.1 )
    f.set_size_inches(9, 18, forward=True)
    colors = sns.color_palette("hls", 20)
    for i, protein in enumerate(proteins):
        for j, drug in enumerate(drugs):
            contacts_per_cluster = data["contacts_per_cluster"][protein][drug]
            contact_residue_labels =  get_labels(contacts_per_cluster)
            weight = num_atoms_per_drug[drug]*frames_per_prot_drug[protein][drug]
            filtered_contact_residue_labels =  filter_less_contacts_than(THRESHOLD, 
                                                                         contact_residue_labels, 
                                                                         data["contacts_per_residue"][protein][drug])
            plot_histogram_to_axis(axes[i,j],
                                   protein, drug, contacts_per_cluster, 
                                   filtered_contact_residue_labels, 
                                   colors, weight, False)
            axes[i,j].xaxis.grid()
            axes[i,j].set_aspect("auto")
            if i == 0:
                axes[i,j].set_title(drug)
            if j == 0:
                axes[i,j].set_ylabel(protein)
            axes[i,j].autoscale(tight=True)
    max_val = 0
    for i, protein in enumerate(proteins):
        for j, drug in enumerate(drugs):
            max_val = max(max_val, axes[i,j].get_ylim()[1])
    for i, protein in enumerate(proteins):
        for j, drug in enumerate(drugs):
            axes[i,j].set_ylim(0,max_val)  
#     plt.show()
    plt.savefig(os.path.join(options.results, "contacts_per_residue_and_cluster_histogram.svg"))
    plt.close()
        
    
    # PLOT CONTACTS PER MOTIF
    contacts_per_motif = calc_contacts_per_motif(proteins, drugs, data, num_atoms_per_drug, frames_per_prot_drug, motifs, ordered_motifs)
    
    f, axes = plt.subplots(len(proteins), len(drugs), sharex='col', sharey='row')
    f.subplots_adjust(hspace=0.1, wspace=0.1 )
    f.set_size_inches(9, 18, forward=True)
    for i, protein in enumerate(proteins):
        for j, drug in enumerate(drugs):
            weight = num_atoms_per_drug[drug]*frames_per_prot_drug[protein][drug]
            colors = sns.color_palette("hls", 20)
            base = numpy.zeros(len(ordered_motifs)) 
            for c,cluster_id in enumerate(contacts_per_motif[protein][drug]):
                row = numpy.zeros(len(ordered_motifs)) 
                for mn, motif_name in enumerate(ordered_motifs):
                    row[mn] = contacts_per_motif[protein][drug][cluster_id][motif_name]
                axes[i,j].bar(range(len(ordered_motifs)), 
                          row, 
                          label = cluster_id,
                          bottom = base,
                          color = colors[c])
                base += row
            if i == 0:
                axes[i,j].set_title(drug)
            if j == 0:
                axes[i,j].set_ylabel(protein)
            axes[i,j].xaxis.grid()
            axes[i,j].autoscale(tight=True)
            axes[i,j].set_autoscaley_on(False)
            axes[i,j].set_xticklabels(ordered_motifs, rotation= 45)
    
    max_val = 0
    for i, protein in enumerate(proteins):
        for j, drug in enumerate(drugs):
            max_val = max(max_val, axes[i,j].get_ylim()[1])
    for i, protein in enumerate(proteins):
        for j, drug in enumerate(drugs):
            axes[i,j].set_ylim(0,max_val)  
    
#     plt.show()
    plt.savefig(os.path.join(options.results, "contacts_per_motif_and_cluster_histogram.svg"))
    plt.close()
            
            
    # COMBO HISTOGRAMS
    common_residue_labels_per_protein = {}
    for protein in proteins:
        all_residue_labels = []
        for drug in drugs: 
            for residue in data["contacts_per_residue"][protein][drug]:
                if data["contacts_per_residue"][protein][drug][residue] >= THRESHOLD:
                    all_residue_labels.append(residue)
        common_residue_labels_per_protein[protein] = order_labels(all_residue_labels)
    
    colors = tools.define_tableau20_cm()
    f, axes = plt.subplots(3, 2, sharey='row')
    f.subplots_adjust(hspace=0.4, wspace=0.1 )
    f.set_size_inches(9, 18, forward=True)
    for i, protein in enumerate(proteins):
        n_groups = len(common_residue_labels_per_protein[protein])
        bar_width = 0.30
        index = numpy.arange(n_groups) + bar_width/2
        legend_labels = []
        for j, drug in enumerate(["CMC","DMA","CMA"]):#drugs):
            print protein, drug
            values = data["contacts_per_residue"][protein][drug]
            weight = num_atoms_per_drug[drug]*frames_per_prot_drug[protein][drug]
            x = numpy.array([ float(values[res]) for res in common_residue_labels_per_protein[protein]]) / weight
            rects = axes[i/2,i%2].bar(index, x, bar_width, color = colors[j])
            index = index + bar_width
            legend_labels.append(drug)
        if i==3:
            lgd = axes[i/2,i%2].legend(legend_labels, loc='right center', bbox_to_anchor=(0, 0))
        axes[i/2,i%2].set_title(protein)
        axes[i/2,i%2].set_xticklabels(common_residue_labels_per_protein[protein], rotation= 45)
        axes[i/2,i%2].set_xticks(numpy.arange(n_groups)+bar_width/2)
        axes[i/2,i%2].autoscale(tight=True)
    max_val = 0
    for i, protein in enumerate(proteins):
        max_val = max(max_val, axes[i/2,i%2].get_ylim()[1])
    print "MAX_VAL", max_val
    for i, protein in enumerate(proteins):
        axes[i/2,i%2].set_ylim(0,max_val) 
        
    #plt.show()
    plt.savefig(os.path.join(options.results, "contacts_per_residue_histogram.svg"))
    
    # COMBO HISTOGRAMS PER MOTIF
    f, axes = plt.subplots(3, 2, sharey='row')
    f.subplots_adjust(hspace=0.4, wspace=0.1 )
    f.set_size_inches(9, 18, forward=True)
    for i, protein in enumerate(proteins):
        n_groups = len(ordered_motifs)
        bar_width = 0.30
        index = numpy.arange(n_groups) + bar_width/2
        legend_labels = []
        # remove "cluster_id" info
        
        contacts_per_motif_wo_c = defaultdict(lambda:defaultdict(int))
        clusters =  contacts_per_motif[protein][drug].keys() 
        for j, drug in enumerate(drugs):
            for cluster in clusters:
                for motif_name in ordered_motifs:
                    contacts_per_motif_wo_c[drug][motif_name] += contacts_per_motif[protein][drug][cluster_id][motif_name]
        
        for j, drug in enumerate(drugs):
            
            values = data["contacts_per_residue"][protein][drug]
            weight = num_atoms_per_drug[drug]*frames_per_prot_drug[protein][drug]
            x = numpy.array([ float(contacts_per_motif_wo_c[drug][motif]) for motif in contacts_per_motif_wo_c[drug]])
            rects = axes[i/2,i%2].bar(index, x, bar_width, color = colors[j])
            index = index + bar_width
            legend_labels.append(drug)
        if i==3:
            lgd = axes[i/2,i%2].legend(legend_labels, loc='right center', bbox_to_anchor=(0, 0))
        axes[i/2,i%2].set_title(protein)
        axes[i/2,i%2].set_xticklabels(ordered_motifs, rotation= 45)
        axes[i/2,i%2].set_xticks(numpy.arange(n_groups)+bar_width/2)
        axes[i/2,i%2].autoscale(tight=True)
    max_val = 0
    for i, protein in enumerate(proteins):
        max_val = max(max_val, axes[i/2,i%2].get_ylim()[1])
    print "MAX_VAL", max_val
    for i, protein in enumerate(proteins):
        axes[i/2,i%2].set_ylim(0,max_val) 
        
#     plt.show()
    plt.savefig(os.path.join(options.results, "contacts_per_motif_histogram.svg"))