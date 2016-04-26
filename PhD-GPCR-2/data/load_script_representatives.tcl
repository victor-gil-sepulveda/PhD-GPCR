display projection orthographic
display ambientoclusion on
display antialias on
display shadows on
color Display Background white
axes location off

material add "my_material"
material change opacity "my_material" 1.0
material change shininess "my_material" 1.0
material change specular "my_material" 0.0
material change outline "my_material" 0.9
material change outlinewidth "my_material" 0.9


# prepare protein
set protein [mol load pdb "%(protein_structure)s"]
mol rename top "%(protein)s"

mol color Name
mol representation NewCartoon 0.300000 12.000000 4.500000 0
mol selection all
mol material Opaque
mol addrep $protein

# Load ligands
set ligands [mol load pdb "%(ligand_structures)s"]
mol rename top "%(drug)s"

# To get the colors
set fp [open "%(cluster_lst)s" r]
set file_data [read $fp]
set data [split $file_data "\n"]

# For each frame (each color)
set i 0
mol delrep 0 $protein
mol delrep 0 $ligands
foreach line $data {
    split $line " "
    set file_color [split $line " "]
    	
	# Get a unique color from the file
	set r [lindex $file_color 0]
	set g [lindex $file_color 1]
	set b [lindex $file_color 2]
	
	#Set colors starting from 17 (current starts in 0)
	set my_color [expr $i + 18]
	puts "COLOR $r $g $b";
	color change rgb $my_color $r $g $b
	graphics $ligands color $my_color
	graphics $ligands materials on
	graphics $ligands material "my_material"
	
	mol selection "resname %(drug)s and noh"
	mol representation Licorice 0.300000 12.000000 12.000000
	mol color ColorId $my_color
	mol addrep $ligands
	#puts "REP ID $rep_id";
	mol drawframes $ligands $i $i

	incr i
}
close $fp

color change rgb white

mol top 0
set current 0

mol color ColorId 3

# show motifs
set my_color 11
%(motifs_code)s

# Save (http://www.ks.uiuc.edu/Research/vmd/mailing_list/vmd-l/19933.html)
#foreach mol [molinfo list] { 
#    # save orientation and zoom parameters 
#      set viewpoints($mol) [molinfo $mol get {center_matrix rotate_matrix scale_matrix global_matrix}] 
#    } 

# Center display port at selection 
set current 1
display resetview

# Load display matrices
%(option_camera)s foreach mol [molinfo list] { 
%(option_camera)s    molinfo $mol set {center_matrix rotate_matrix scale_matrix global_matrix} {%(camera_settings)s}
%(option_camera)s } 

# Render
%(option_camera)s render Tachyon %(pre_render_file)s "/usr/local/lib/vmd/tachyon_LINUX -aasamples 12 %(pre_render_file)s -add_skylight 1.5  -format PSD48 -res 1024 1024 -o %(rendered_file)s"

# Load zoomed display matrices
%(option_zoom)s foreach mol [molinfo list] { 
%(option_zoom)s    molinfo $mol set {center_matrix rotate_matrix scale_matrix global_matrix} {%(camera_settings_zoomed)s}
%(option_zoom)s } 

# Render

%(option_zoom)s render Tachyon %(pre_render_zoom_file)s "/usr/local/lib/vmd/tachyon_LINUX -aasamples 12 %(pre_render_zoom_file)s -add_skylight 1.5  -format PSD48 -res 1024 1024 -o %(rendered_zoom_file)s"

exit
