"""
Created on 14/8/2015

@author: victor
"""
from _collections import defaultdict
import numpy
import matplotlib.pyplot as pl
import itertools
from tools import define_tableau20_cm
from matplotlib.font_manager import FontProperties

def parse_contacts(contacts_file_path):
    with open(contacts_file_path,"r") as contacts_file:
        contacts_per_cluster = {}
        contacts_per_residue = defaultdict(int)
        for line in contacts_file:
            parts = line.split()
            contacts_per_cluster[parts[0]]= parts[1:-1]
            for res in parts[1:-1]:
                contacts_per_residue[res] += 1
    return contacts_per_cluster, contacts_per_residue

def count_contacts_per_cluster_and_residue(contacts_per_cluster, contact_residue_labels, normalize = False):
    cluster_ids = contacts_per_cluster.keys()
    
    contact_to_data_id = dict(zip(contact_residue_labels, range(len(contact_residue_labels))))
    
    # Generate the data counts, each row is a cluster id
    data = numpy.zeros(( len(cluster_ids),len(contact_residue_labels)))
    
    for i, cluster_id in enumerate(cluster_ids):
        for contact_id in contacts_per_cluster[cluster_id]:
            try:
                j = contact_to_data_id[contact_id]
                data[i][j] += 1
            except KeyError:
                pass # This means we are using a subset of the labels
    
    if normalize:
        area = data.sum().sum()
        # Then average the values
        data /= area
    
    return data

def get_labels (contacts_per_cluster):
    cluster_ids = contacts_per_cluster.keys()
    all_contacts_tags = []
    for cid in cluster_ids:
        all_contacts_tags.extend(contacts_per_cluster[cid])
    all_tags = list(set(all_contacts_tags)); all_tags.sort()
    return all_tags
    
def filter_less_contacts_than(num_contacts, contact_labels, contacts_per_residue):
    filtered_contacts_residue_labels = []
    for contact_label in contact_labels:
        if contacts_per_residue[contact_label] >= num_contacts:
            filtered_contacts_residue_labels.append(contact_label)
    
    
    return filtered_contacts_residue_labels

def get_num_contacts_per_residue(contacts_per_cluster):
    contacts_per_residue = defaultdict(int)
    for cluster_id in contacts_per_cluster:
        for contact_id in contacts_per_cluster[cluster_id]:
            contacts_per_residue[contact_id] += 1
    return contacts_per_residue
    
def plot_histogram(contacts_per_cluster, contact_residue_labels, histogram_path, normalize = False ):
    cluster_ids = contacts_per_cluster.keys()
    
    data = count_contacts_per_cluster_and_residue(contacts_per_cluster, contact_residue_labels, normalize)
    
    # plot data    
    _, ax = pl.subplots()
    tab20cm = itertools.cycle(define_tableau20_cm())
    num_res_labels = len(contact_residue_labels)
    base = numpy.zeros(num_res_labels) 
    for i, row in enumerate(data):
        if row.sum() > 0:
            pl.bar(range(num_res_labels), 
                      row, 
                      label = cluster_ids[i],
                      bottom = base,
                      color = next(tab20cm))
            base += row
    ax.set_xticklabels(contact_residue_labels, rotation= 45)
    pl.xticks(range(len(contact_residue_labels)))
    fontP = FontProperties()
    fontP.set_size('small')
    pl.legend( ncol = 4, prop = fontP,fancybox=True, shadow=True, bbox_to_anchor=(0., 1.02, 1., .102), loc=3, mode="expand",borderaxespad=0.)
    pl.subplots_adjust(bottom=0.15, top=0.8) 
#     pl.show()
    pl.savefig(histogram_path)
    

if __name__ == "__main__":
    contacts_per_cluster, contacts_per_residue = parse_contacts("/media/victor/d9c4538f-7e46-492d-812b-8110c78c0af9/James/James_3/Trajectories_ParaVictor/Test/BVDV_CMA/residue_contacts.txt")
    
    print contacts_per_residue
    print get_num_contacts_per_residue(contacts_per_cluster)
    
    contact_residue_labels =  get_labels (contacts_per_cluster)
    
    filtered_contact_residue_labels =  filter_less_contacts_than(2000, contact_residue_labels, contacts_per_residue)
    
    cpcr =  count_contacts_per_cluster_and_residue(contacts_per_cluster, contact_residue_labels, 100)
    
#     print contacts_per_residue
#     print zip(contact_residue_labels, cpcr.sum(0))
#     
#     plot_histogram(contacts_per_cluster, filtered_contact_residue_labels, "", True)
    
