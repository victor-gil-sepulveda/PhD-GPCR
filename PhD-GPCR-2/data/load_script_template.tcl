#*******************
# http://www.ks.uiuc.edu/Research/vmd/current/ug/node127.html
#*******************
proc vmd_draw_arrow {mol start end} {
    # an arrow is made of a cylinder and a cone
    set middle [vecadd $start [vecscale 0.9 [vecsub $end $start]]]
    graphics $mol cylinder $start $middle radius 0.3
    graphics $mol cone $middle $end radius 0.5
}

proc label_atom {selection_string label_string} {
    set sel [atomselect top $selection_string]
    if {[$sel num] != 1} {
        error "label_atom: '$selection_string' must select 1 atom"
    }
    # get the coordinates of the atom
    lassign [$sel get {x y z}] coord
    # and draw the text
    draw text $coord $label_string
}

#*******************

set fp [open "%(cluster_lst)s" r]
set file_data [read $fp]
set data [split $file_data "\n"]

display projection orthographic
display ambientoclusion on
display antialias on
display shadows on
color Display Background white
axes location off

material add "my_material"
material change opacity "my_material" 1.0 
material change shininess "my_material" 1.0 
material change outline "my_material" 0.9
material change outlinewidth "my_material" 0.9 


# prepare protein
set current [mol load pdb "%(main_structure)s"]
mol rename top "%(main_structure)s"

# load ligand clusters
foreach line $data {
    	split $line " "
    	set file_color [split $line " "]
    	set filename [lindex $file_color 0]
	
	set current [mol load pdb $filename]
	mol rename top $filename
	mol showrep $current 0 off
	
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
	
	
	# Add balls in place of ligands
#	 set n [molinfo $current get numframes]
#	 for { set i 0 } { $i < $n } { incr i } {
#	     set ligand [atomselect $current "resname CMA and noh" frame $i ]
#	     set lig_center [measure center $ligand]
#	     graphics $current sphere $lig_center 
#	 }

%(option_arrow)s	# Add arrows in place of ligands 
%(option_arrow)s	set n [molinfo $current get numframes]
%(option_arrow)s	for { set i 0 } { $i < $n } { incr i } {
%(option_arrow)s		set atom_begin [atomselect $current "%(arrow_begin_selection)s" frame $i ]
%(option_arrow)s		set begin_coords [lindex [$atom_begin get {x y z}] 0]
%(option_arrow)s		set atom_end [atomselect $current "%(arrow_end_selection)s" frame $i ]
%(option_arrow)s		set end_coords [lindex [$atom_end get {x y z}] 0]
%(option_arrow)s		vmd_draw_arrow top $begin_coords $end_coords
%(option_arrow)s	}
	
%(option_drug)s	# In case we need to show the actual ligand...
%(option_drug)s	mol selection "resname %(drug)s and noh"
%(option_drug)s	mol representation Licorice 0.300000 12.000000 12.000000
%(option_drug)s	mol color ColorId $my_color
%(option_drug)s mol addrep $current

}
close $fp

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
%(option_camera)srender Tachyon %(pre_render_file)s "/usr/local/lib/vmd/tachyon_LINUX -fullshade -aasamples 12 %(pre_render_file)s -trans_vmd -rescale_lights 0.3 -add_skylight 1.0  -normalize -format PSD48 -res 1024 1024 -o %(rendered_file)s"

# Load zoomed display matrices
%(option_camera)s foreach mol [molinfo list] { 
%(option_camera)s    molinfo $mol set {center_matrix rotate_matrix scale_matrix global_matrix} {%(camera_settings_zoomed)s}
%(option_camera)s } 

# Render
%(option_camera)srender Tachyon %(pre_render_zoom_file)s "/usr/local/lib/vmd/tachyon_LINUX -fullshade -aasamples 12 %(pre_render_zoom_file)s -trans_vmd -rescale_lights 0.3 -add_skylight 1.0  -normalize -format PSD48 -res 1024 1024 -o %(rendered_zoom_file)s"

