'''
Created on 18/9/2015

@author: victor
'''

import os
import tools
from histogram import get_labels, filter_less_contacts_than,\
    parse_drug_info, parse_frames_info,\
    plot_histogram_to_axis, order_labels,\
    calc_contacts_per_cluster_per_motif_for_prots_and_drugs, parse_motifs,\
    calc_contacts_per_motif, parse_contacts_file, scale_axes_to_max_val,\
    get_total_contacts_per_protein_and_drug
from optparse import OptionParser
import seaborn as sns
import matplotlib.pyplot as plt
import numpy
from prody import proteins
import pandas as pd

def check_options(options):
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

if __name__ == '__main__':
    sns.set_style("whitegrid")
    
    parser = OptionParser()
    parser.add_option("--contacts", dest="contacts")
    parser.add_option("--motifs", dest="motifs")
    parser.add_option("--results", dest="results")
    parser.add_option("--drug-atoms", dest="drug_atoms")
    parser.add_option("--frames-per-drug", dest="frames_per_drug")
    
    (options, args) = parser.parse_args()
    check_options(options)
    
    tools.create_dir(options.results)
    
    THRESHOLD = 1000
    
    # Parse all files
    num_atoms_per_drug = parse_drug_info(options.drug_atoms)
    frames_per_prot_drug = parse_frames_info(options.frames_per_drug)
    data, proteins, drugs = parse_contacts_file(options.contacts)
    motifs = parse_motifs(options.motifs)
    
    # Override order for stuff
    ordered_motifs = ['F/I-II', 'A/IV', 'B/V', 'C/VI', 'D', 'E/VII', 'Priming Loop']
    drugs.sort() 
    proteins = ["JEV", "WNV", "TBEV", "BVDV", "HCV", "Polio"]
    
    # PLOT TOTAL NUMBER OF CONTACTS PER PROTEIN
    t_c_p_d = get_total_contacts_per_protein_and_drug(data["contacts_per_residue"], 
                                                      proteins, drugs,
                                                      num_atoms_per_drug, frames_per_prot_drug)
    db = pd.DataFrame.from_dict(t_c_p_d)
    sns.barplot(x="Protein", y="Contacts", hue="Drug", data = db)
    plt.savefig(os.path.join(options.results, "total_contacts.svg"))
    
    # PLOT CONTACTS PER RESIDUE
    f, axes = plt.subplots(len(proteins), len(drugs), sharey='row')
    f.subplots_adjust(hspace=0.3, wspace=0.1 )
    f.set_size_inches(9, 18, forward=True)
    colors = sns.color_palette("hls", 15)
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
    
    scale_axes_to_max_val(axes, len(proteins), len(drugs))
    plt.savefig(os.path.join(options.results, "contacts_per_residue_and_cluster_histogram.svg"))
    plt.close()
        
    
    # PLOT CONTACTS PER MOTIF
    contacts_per_motif_and_cluster = calc_contacts_per_cluster_per_motif_for_prots_and_drugs(proteins, drugs, data, num_atoms_per_drug, frames_per_prot_drug, motifs, ordered_motifs)
    
    f, axes = plt.subplots(len(proteins), len(drugs), sharex='col', sharey='row')
    f.subplots_adjust(hspace=0.1, wspace=0.1 )
    f.set_size_inches(9, 18, forward=True)
    colors = sns.color_palette("hls", 15)
    for i, protein in enumerate(proteins):
        for j, drug in enumerate(drugs):
            weight = num_atoms_per_drug[drug]*frames_per_prot_drug[protein][drug]
            base = numpy.zeros(len(ordered_motifs)) 
            for c,cluster_id in enumerate(contacts_per_motif_and_cluster[protein][drug]):
                row = numpy.zeros(len(ordered_motifs)) 
                for mn, motif_name in enumerate(ordered_motifs):
                    row[mn] = contacts_per_motif_and_cluster[protein][drug][cluster_id][motif_name]
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
    
    scale_axes_to_max_val(axes, len(proteins), len(drugs))
    
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
    
    colors = sns.color_palette("deep")#tools.define_tableau20_cm()
    f, axes = plt.subplots(3, 2, sharey='row')
    f.subplots_adjust(hspace=0.4, wspace=0.1 )
    f.set_size_inches(10, 18, forward=True)
    for i, protein in enumerate(proteins):
        n_groups = len(common_residue_labels_per_protein[protein])
        bar_width = 0.30
        legend_labels = []
        index = numpy.arange(n_groups) + bar_width/2
        for j, drug in enumerate(drugs):
            values = data["contacts_per_residue"][protein][drug]
            weight = num_atoms_per_drug[drug]*frames_per_prot_drug[protein][drug]
            x = numpy.array([ float(values[res]) for res in common_residue_labels_per_protein[protein]]) / weight
            axes[i/2,i%2].bar(index, x, bar_width, color = colors[j])
            index = index + bar_width
            legend_labels.append(drug)
        if i==3:
            lgd = axes[i/2,i%2].legend(legend_labels, loc='center right', bbox_to_anchor=(1.3, 0.5))
        axes[i/2,i%2].xaxis.grid()
        axes[i/2,i%2].set_title(protein)
        axes[i/2,i%2].set_xticklabels(common_residue_labels_per_protein[protein], rotation= 45)
        axes[i/2,i%2].set_xticks(numpy.arange(n_groups)+bar_width/2)
        axes[i/2,i%2].autoscale(tight=True)
    scale_axes_to_max_val(axes, 3, 2)
    plt.savefig(os.path.join(options.results, "contacts_per_residue_histogram.svg"),
                bbox_extra_artists=(lgd,))
    
    # COMBO HISTOGRAMS PER MOTIF
    f, axes = plt.subplots(3, 2, sharey='row')
    f.subplots_adjust(hspace=0.4, wspace=0.1 )
    f.set_size_inches(10, 18, forward=True)
    for i, protein in enumerate(proteins):
        n_groups = len(ordered_motifs)
        bar_width = 0.30
        index = numpy.arange(n_groups) + bar_width/2
        legend_labels = []
        # remove "cluster_id" info
        
        for j, drug in enumerate(drugs):
            contacts_per_motif = calc_contacts_per_motif(contacts_per_motif_and_cluster[protein][drug], ordered_motifs)
            values = data["contacts_per_residue"][protein][drug]
            weight = num_atoms_per_drug[drug]*frames_per_prot_drug[protein][drug]
            x = numpy.array([ float(contacts_per_motif[motif]) for motif in ordered_motifs])
            rects = axes[i/2,i%2].bar(index, x, bar_width, color = colors[j]) 
            index = index + bar_width
            legend_labels.append(drug)
        if i==3:
            lgd = axes[i/2,i%2].legend(legend_labels, loc='center right', bbox_to_anchor=(1.3, 0.5))
        axes[i/2,i%2].xaxis.grid()
        axes[i/2,i%2].set_title(protein)
        axes[i/2,i%2].set_xticklabels(ordered_motifs, rotation= 45)
        axes[i/2,i%2].set_xticks(numpy.arange(n_groups)+bar_width/2)
#         axes[i/2,i%2].autoscale(tight=True)
    scale_axes_to_max_val(axes, 3, 2)
    plt.savefig(os.path.join(options.results, "contacts_per_motif_histogram.svg"),
                bbox_extra_artists=(lgd,))