'''
Created on Dec 3, 2015

@author: victor
'''
from _collections import defaultdict
from optparse import OptionParser
import glob
import os
import numpy
from pyproct.tools.scriptTools import create_directory
from prody.proteins.pdbfile import parsePDB, writePDB

def gen_motifs_code(motifs, motif_colors, draw_text):
    motif_code_strings = [] 
    motif_code_template = """mol selection "%(motif_selection)s"
mol representation NewCartoon 0.350000 12.000000 4.500000 0
color change rgb $my_color %(motif_color)s %(motif_color)s %(motif_color)s 
mol color ColorId $my_color
mol material "my_material"
mol addrep $current
%(option_text)s draw text [measure center [atomselect $current "%(motif_selection)s"]] "%(motif_name)s"  size 1 thickness 2
set my_color [expr $my_color + 1]
"""
    for motif_name in motifs:
        motif_code_strings.append(motif_code_template%{
                                                       "motif_selection":motifs[motif_name],
                                                       "motif_name": motif_name,
                                                       "motif_color": motif_colors[motif_name],
                                                       "option_text": "" if draw_text else "#" 
                                                       }) 
    return "".join(motif_code_strings)

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-i", dest="input")
    parser.add_option("-o", dest="output_folder")
    
    parser.add_option("--motifs", dest="motifs")
    parser.add_option("--camera", dest="camera")
    parser.add_option("-z", action="store_true", default = False, dest="do_zoom")
    (options, args) = parser.parse_args()
    
    if options.input is None:
        parser.error('Input file is needed. Line format: "protein_name drug_name clusters_folder"')
    
    if options.motifs is None:
        parser.error('You must specify the motifs file')

    # motifs
    motifs = defaultdict(dict)
    for line in open(options.motifs):
        if line[0]!= "#":
            parts = [ s.strip() for s in line.split(",")]
            motifs[parts[0]][parts[1]] = parts[2]
    
    # prepare motif colors
    motif_colors = {}
    start = 0.1 ; end = 0.9
    for c, motif_id in zip(numpy.arange(start, end+0.01, (end - start) / len(motifs)),motifs[motifs.keys()[0]].keys()):
        motif_colors[motif_id] = c
    
    # camera setups
    if options.camera is not None:
        camera = {}
        for line in open(options.camera):
            if line[0] != "#":
                parts = [ s.strip() for s in line.split(":")]
                # general and zoomed
                camera[parts[0]]=(parts[1], parts[2])
    
    # load colors per cluster
    #colors = [(1.,0.,0.), (0.,1.,0.), (0.,0.,1.)]*10
    import seaborn as sns
    colors = sns.hls_palette(15, l=.3, s=.8)
    
    # VMD execution template
    template = open("/home/victor/git/PhD-GPCR/PhD-GPCR-2/data/load_script_representatives.tcl").read()
    
    for line in open(options.input):
        protein, drug, folder = line.strip().split()

        # sorted clusters and same color generation always make the same cluster_id, color pair
        representatives_file = os.path.join(folder, "representatives.pdb")
        
        output_folder = os.path.join(options.output_folder, drug, protein)
        create_directory(output_folder)
        
        pdb = parsePDB(representatives_file)
        writePDB(os.path.join(output_folder,"protein.pdb"), pdb.select("protein"), csets = [0])
        writePDB(os.path.join(output_folder,"ligands.pdb"), pdb.select("resname %s"%drug))
        
        num_clusters = pdb.numCoordsets()
        clusters_file = open(os.path.join(output_folder,"cluster_colors"), "w")
        for i in range(num_clusters):
            clusters_file.write("%.2f %.2f %.2f%s"%(   colors[i][0],
                                                       colors[i][1],
                                                       colors[i][2],
                                ("\n" if i <(num_clusters-1) else "")))
        clusters_file.close()
        
        camera_settings = ""; camera_settings_zoomed = ""; option_camera = "#"; pre_render_file = ""; rendered_file = ""; option_zoom = "#"
        if options.camera is not None:
            camera_settings = camera[protein][0]
            camera_settings_zoomed = camera[protein][1]
            option_camera = "" 
            pre_render_file = os.path.join(options.output_folder, drug, protein, "%s_%s_render.dat"%(protein, drug))
            rendered_file = os.path.join(options.output_folder, drug, protein, "%s_%s_render.psd"%(protein, drug))
            pre_render_zoom_file = os.path.join(options.output_folder, drug, protein, "%s_%s_zoom_render.dat"%(protein, drug))
            rendered_zoom_file = os.path.join(options.output_folder, drug, protein, "%s_%s_zoom_render.psd"%(protein, drug))
            option_zoom = ""
        
        
        if not options.do_zoom:
            option_zoom = "#"
        else:
            option_zoom = ""
        
        vmd_file_contents = template%{"cluster_lst": os.path.join(output_folder,"cluster_colors"), 
                                      "protein_structure":os.path.join(output_folder,"protein.pdb"),
                                      "ligand_structures":os.path.join(output_folder,"ligands.pdb"),
                                      "motifs_code": gen_motifs_code(motifs[protein], motif_colors, False),
                                      "viewpoint_values": "",
                                      "protein": protein,
                                      "drug": drug,
                                      "camera_settings":camera_settings,
                                      "camera_settings_zoomed": camera_settings_zoomed,
                                      "option_camera": option_camera,
                                      "pre_render_file": pre_render_file,
                                      "rendered_file": rendered_file,
                                      "pre_render_zoom_file": pre_render_zoom_file,
                                      "rendered_zoom_file": rendered_zoom_file,
                                      "option_zoom": option_zoom
        }
            
        vmd_vis_filename = os.path.join(output_folder,"vmd_vis")
        open(vmd_vis_filename,"w").write(vmd_file_contents)
        print vmd_vis_filename
