"""
Created on Dec 8, 2015

@author: victor
"""
import unittest
import os
import histogram

import Test
from histogram import get_labels, filter_less_contacts_than,\
    count_contacts_per_cluster_and_residue, calc_prot_drug_contacts_per_motif
import numpy
from _collections import defaultdict

class TestDataManaging(unittest.TestCase):
    contacts_file_name = os.path.join(Test.__path__[0], "data", "DMA_WNV.txt")
    small_contacts_file_name = os.path.join(Test.__path__[0], "data", "contacts_small.txt")

    def test_load_data(self):
        contacts_per_cluster, contacts_per_residue = histogram.parse_contacts(TestDataManaging.small_contacts_file_name)
        number_of_items = {
                            "cluster_4": 16, 
                            "cluster_2": 16,
                            "cluster_7": 16,
                            "cluster_5": 6,
                            "cluster_6": 16,
                            "cluster_0": 16,
                            "cluster_1": 16,
                            "cluster_3": 16,
                            "cluster_8": 6
                            }
        for cluster_id in contacts_per_cluster:
            self.assertEqual(number_of_items[cluster_id], len(contacts_per_cluster[cluster_id]),"%s is suspicious"%cluster_id)
        correct_contacts_per_residue = {'181:TYR:A': 4,  '395:ASP:A': 17, '421:LYS:A': 7, '268:ASP:A': 5, 
         '442:SER:A': 15, '461:ARG:A': 1, '337:TYR:A': 6, '333:GLN:A': 8, 
         '528:SER:A': 18, '340:ASN:A': 13, '266:GLY:A': 22, '527:TRP:A': 3, '329:ARG:A': 5}
        
        self.assertDictEqual(correct_contacts_per_residue, contacts_per_residue)
        
    def test_get_labels(self):
        contacts_per_cluster, contacts_per_residue = histogram.parse_contacts(TestDataManaging.small_contacts_file_name)
        labels = get_labels(contacts_per_cluster)
        self.assertItemsEqual(labels, contacts_per_residue.keys())
        ordered_labels = ['181:TYR:A', '266:GLY:A', '268:ASP:A', '329:ARG:A', '333:GLN:A', '337:TYR:A', '340:ASN:A', '395:ASP:A', '421:LYS:A', '442:SER:A', '461:ARG:A', '527:TRP:A', '528:SER:A']
        self.assertSequenceEqual(ordered_labels, labels)
        
    def test_filter(self):
        _, contacts_per_residue = histogram.parse_contacts(TestDataManaging.small_contacts_file_name)
        ordered_labels = ['181:TYR:A', '266:GLY:A', '268:ASP:A', '329:ARG:A', '333:GLN:A', '337:TYR:A', '340:ASN:A', '395:ASP:A', '421:LYS:A', '442:SER:A', '461:ARG:A', '527:TRP:A', '528:SER:A']
        filtered_labels =  filter_less_contacts_than( 15, ordered_labels, contacts_per_residue)
        self.assertItemsEqual(['266:GLY:A', '395:ASP:A', '442:SER:A', '528:SER:A'], filtered_labels)
    
    def test_count_contacts_per_cluster_and_residue(self):
        contacts_per_cluster, _ = histogram.parse_contacts(TestDataManaging.small_contacts_file_name)
        filtered_labels = ['266:GLY:A', '395:ASP:A', '442:SER:A', '528:SER:A']
        ordered_clusters = ["cluster_0",
                            "cluster_1",
                            "cluster_2",
                            "cluster_3",
                            "cluster_4", 
                            "cluster_5",
                            "cluster_6",
                            "cluster_7",
                            "cluster_8"]
        data = count_contacts_per_cluster_and_residue(
                                                     contacts_per_cluster, 
                                                     ordered_clusters,
                                                     filtered_labels, 
                                                     1., 
                                                     False)
        correct_data = [[  0., 7., 7., 0.],
                        [ 10., 0., 2., 0.],
                        [  0., 4., 0., 3.],
                        [  0., 0., 0., 6.],
                        [  0., 6., 0., 4.],
                        [  0., 0., 6., 0.],
                        [  1., 0., 0., 4.],
                        [ 11., 0., 0., 1.],
                        [  0., 0., 0., 0.]]
        numpy.testing.assert_equal(correct_data, data)
    
#     def test_parse_motifs
    
    def test_calc_prot_drug_contacts_per_motif(self):
        contacts_per_cluster, _ = histogram.parse_contacts(TestDataManaging.small_contacts_file_name)

        motifs = defaultdict(dict)
        for line in open(os.path.join(Test.__path__[0], "data", "motifs.txt")):
            if line[0]!= "#":
                parts = [ s.strip() for s in line.split(",")]
                sel_parts = parts[2].split()
                motifs[parts[0]][parts[1]] = (int(sel_parts[1]),int(sel_parts[3]))
        print motifs["WNV"]
        
        ordered_motifs = ['F/I-II', 'A/IV', 'B/V', 'C/VI', 'D', 'E/VII', 'Priming Loop']
   
        print calc_prot_drug_contacts_per_motif(contacts_per_cluster, 
                                                motifs["WNV"], 
                                                ordered_motifs, 
                                                1.)
        
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()