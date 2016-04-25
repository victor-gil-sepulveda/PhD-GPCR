'''
Created on Dec 3, 2015

@author: victor
'''
from _collections import defaultdict
from optparse import OptionParser
import glob
import os
import numpy

def gen_motifs_code(motifs, motif_colors, draw_text):
    motif_code_strings = [] 
    motif_code_template = """mol selection "%(motif_selection)s"
mol representation NewCartoon 0.300000 12.000000 4.500000 0
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
    parser.add_option("--arrows", dest="arrows")
    parser.add_option("--motifs", dest="motifs")
    parser.add_option("--camera", dest="camera")
    (options, args) = parser.parse_args()
    
    if options.input is None:
        parser.error('Input file is needed. Line format: "protein_name drug_name clusters_folder"')
    
    if options.motifs is None:
        parser.error('You must specify the motifs file')

    # arrow selections
    if options.arrows is not None:
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
    template = open("/home/victor/git/PhD-GPCR/PhD-GPCR-2/data/load_script_template.tcl").read()
    
    for line in open(options.input):
        protein, drug, folder = line.strip().split()

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
        
        if options.arrows is not None:
            option_arrow = ""
            option_drug = "#"
            arrow_begin = arrows[drug]["begin"]
            arrow_end = arrows[drug]["end"]
        else:
            option_arrow = "#"
            option_drug = ""
            arrow_begin = ""
            arrow_end = ""
        
        camera_settings = ""; camera_settings_zoomed = ""; option_camera = "#"; pre_render_file = ""; rendered_file = ""; option_zoom = "#"
        if options.camera is not None:
            camera_settings = camera[protein][0]
            camera_settings_zoomed = camera[protein][1]
            option_camera = "" 
            pre_render_file = os.path.join(folder,"%s_%s_render.dat"%(protein, drug))
            rendered_file = os.path.join(folder,"%s_%s_render.psd"%(protein, drug))
            pre_render_zoom_file = os.path.join(folder,"%s_%s_zoom_render.dat"%(protein, drug))
            rendered_zoom_file = os.path.join(folder,"%s_%s_zoom_render.psd"%(protein, drug))
            option_zoom = ""
        #OVERRIDE OPTION ZOOM
        option_zoom = "#"
        
        vmd_file_contents = template%{"cluster_lst": clusters_filename, 
                                      "main_structure":os.path.join(folder,"%s.pdb"%protein),
                                      "arrow_begin_selection": arrow_begin,
                                      "arrow_end_selection": arrow_end,
                                      "motifs_code": gen_motifs_code(motifs[protein], motif_colors, False),
                                      "viewpoint_values": "",
                                      "drug": drug,
                                      "option_arrow": option_arrow,
                                      "option_drug": option_drug,
                                      "camera_settings":camera_settings,
                                      "camera_settings_zoomed": camera_settings_zoomed,
                                      "option_camera": option_camera,
                                      "pre_render_file": pre_render_file,
                                      "rendered_file": rendered_file,
                                      "pre_render_zoom_file": pre_render_zoom_file,
                                      "rendered_zoom_file": rendered_zoom_file,
                                      "option_zoom": option_zoom
        }
            
        vmd_vis_filename = os.path.join(folder,"vmd_vis")
        open(vmd_vis_filename,"w").write(vmd_file_contents)
        print vmd_vis_filename
