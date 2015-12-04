#*******************
# http://www.ks.uiuc.edu/Research/vmd/current/ug/node127.html
#*******************
proc vmd_draw_arrow {mol start end} {
    # an arrow is made of a cylinder and a cone
    set middle [vecadd $start [vecscale 0.9 [vecsub $end $start]]]
    graphics $mol cylinder $start $middle radius 0.15
    graphics $mol cone $middle $end radius 0.25
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
axes location off

material add "balls"
material change opacity "balls" 0.9 
material change shininess "balls" 0.3 

# prepare protein
set current [mol load pdb "%(main_structure)s"]
mol rename top "%(main_structure)s"
mol showrep $current 0 off
#mol color ColorID 4
#mol selection "protein"
#mol representation NewCartoon 0.300000 12.000000 4.500000 0
#mol addrep $current

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
	
	color change rgb $current $r $g $b
	graphics $current color $current
	graphics $current materials on
	graphics $current material "balls"
	
	
	# Add balls in place of ligands
#	 set n [molinfo $current get numframes]
#	 for { set i 0 } { $i < $n } { incr i } {
#	     set ligand [atomselect $current "resname CMA and noh" frame $i ]
#	     set lig_center [measure center $ligand]
#	     graphics $current sphere $lig_center 
#	 }

	# Add arrows in place of ligands 
	set n [molinfo $current get numframes]
	for { set i 0 } { $i < $n } { incr i } {
		set atom_begin [atomselect $current "%(arrow_begin_selection)s" frame $i ]
		set begin_coords [lindex [$atom_begin get {x y z}] 0]
		set atom_end [atomselect $current "%(arrow_end_selection)s" frame $i ]
		set end_coords [lindex [$atom_end get {x y z}] 0]
		vmd_draw_arrow top $begin_coords $end_coords
	}
	
	# In case we need to show the actual ligand...
	#mol selection "resname CMA and noh"
	#mol representation Licorice 0.300000 12.000000 12.000000
    #mol addrep $current

    mol showrep $current 0 off
	
}
close $fp

mol top 0
set current 0

# show motifs
%(motifs_code)s

# Center display port at selection 
display resetview
