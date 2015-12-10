'''
Created on Dec 3, 2015

@author: victor
'''
from _collections import defaultdict
from optparse import OptionParser
import glob
import os

def gen_motifs_code(motifs):
    motif_code_strings = [] 
    motif_code_template = """mol selection "%(motif_selection)s"
#mol representation Licorice 0.300000 12.000000 12.000000
mol representation NewCartoon 0.300000 12.000000 4.500000 0
mol addrep $current
draw text [measure center [atomselect $current "%(motif_selection)s"]] "%(motif_name)s"  size 1 thickness 2
"""
    for motif_name in motifs:
        motif_code_strings.append(motif_code_template%{"motif_selection":motifs[motif_name],
                                                       "motif_name": motif_name})
    
    return "".join(motif_code_strings)


if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-i", dest="input")
    parser.add_option("--arrows", dest="arrows")
    parser.add_option("--motifs", dest="motifs")
    (options, args) = parser.parse_args()
    
    if options.input is None:
        parser.error('Input file is needed. Line format: "protein_name drug_name clusters_folder"')
    
    if options.arrows is None:
        parser.error('You must specify the arrow definition file.')
    
    if options.motifs is None:
        parser.error('You must specify the motifs file')

    # arrow selections
    arrows = defaultdict(dict)
    for line in open(options.arrows):
        parts = line.split(",")
        arrows[parts[0]]["begin"] = parts[1]
        arrows[parts[0]]["end"] = parts[2]
    
    # motifs
    motifs = defaultdict(dict)
    for line in open(options.motifs):
        if line[0]!= "#":
            parts = [ s.strip() for s in line.split(",")]
            motifs[parts[0]][parts[1]] = parts[2]
    
    # load colors per cluster
    #colors = [(1.,0.,0.), (0.,1.,0.), (0.,0.,1.)]*10
    import seaborn as sns
    colors = sns.hls_palette(15, l=.3, s=.8)
    
    # VMD execution template
    template = open("/home/victor/git/PhD-GPCR/PhD-GPCR-2/data/load_script_template.tcl").read()
    
    for line in open(options.input):
        protein, drug, folder = line.strip().split()
        print folder
        # sorted clusters and same color generation always make the same cluster_id, color pair
        cluster_files = sorted(glob.glob(os.path.join(folder, "ligand_cluster_*.pdb")))
        clusters_filename = os.path.join(folder, "clusters.lst")
        clusters_file = open(clusters_filename, "w")
        for i, cf in enumerate(sorted(cluster_files)):
            clusters_file.write("%s %.2f %.2f %.2f%s"%(cf, 
                                                       colors[i][0],
                                                       colors[i][1],
                                                       colors[i][2],
                                ("\n" if i <(len(cluster_files)-1) else "")))
        clusters_file.close()
        
        vmd_file_contents = template%{"cluster_lst": clusters_filename, 
                                  "main_structure":os.path.join(folder,"%s.pdb"%protein),
                                  "arrow_begin_selection": arrows[drug]["begin"],
                                  "arrow_end_selection": arrows[drug]["end"],
                                  "motifs_code": gen_motifs_code(motifs[protein])}
        vmd_vis_filename = os.path.join(folder,"vmd_vis")
        open(vmd_vis_filename,"w").write(vmd_file_contents)
        print vmd_vis_filename
