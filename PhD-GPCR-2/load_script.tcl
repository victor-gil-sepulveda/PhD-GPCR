set fp [open "clusters.lst" r]
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
	set n [molinfo $current get numframes]
	for { set i 0 } { $i < $n } { incr i } {
	    set ligand [atomselect $current "resname CMA and noh" frame $i ]
	    set lig_center [measure center $ligand]
	    graphics $current sphere $lig_center 
	}
	
	# In case we need to show the actual ligand...
	#mol selection "resname CMA and noh"
	#mol representation Licorice 0.300000 12.000000 12.000000
    #mol addrep $current

    # Show only the first guy
    if {
        $current == 0
    } then {
	    mol color ColorID 4
	    mol selection "protein"
	    mol representation NewCartoon 0.300000 12.000000 4.500000 0
	    mol addrep $current
	} else {
	    mol showrep $current 0 off
	}
}
#8934031100001545261
close $fp

mol top 0

