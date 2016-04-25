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
set current [mol load pdb "%(main_structure)s"]
mol rename top "%(main_structure)s"

set nframes [molinfo top get numframes] 

# For each frame
	# Get a unique color from the file
	set r [lindex $file_color 1]
	set g [lindex $file_color 2]
	set b [lindex $file_color 3]
	
	#Set colors starting from 17 (current starts in 1)
	set my_color [expr $current + 17]
	color change rgb $my_color $r $g $b
	graphics $current color $my_color
	graphics $current materials on
	graphics $current material "my_material"
	
	mol selection "resname %(drug)s and noh"
	mol representation Licorice 0.300000 12.000000 12.000000
	mol color ColorId $my_color
	set rep_id [mol addrep $current]

	mol drawframes top $rep_id $frame_index

color change rgb white

mol top 0
set current 0

# Hide full prot
mol showrep $current 0 off
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
