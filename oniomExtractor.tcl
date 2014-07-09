
#### PROCEDURES

proc rect_to_spherical {x y z} {
    list [set rad [expr {hypot($x, hypot($y, $z))}]] [expr {atan2($y,$x)}] [expr {acos($z/($rad+1.0e-20))}]
}

proc spherical_to_rect {rad phi theta} {
    list [expr {$rad * cos($phi) * sin($theta)}] [expr {$rad * sin($phi) * sin($theta)}] [expr {$rad * cos($theta)}]
}


proc progress_init {tot} {
    set ::progress_start     [clock seconds]
    set ::progress_last      0
    set ::progress_last_time 0
    set ::progress_tot       $tot
}

proc progress_tick {file cur} {
   set now [clock seconds]
   set tot $::progress_tot

   if {$cur > $tot} {
       set cur $tot
   }
   if {($cur >= $tot && $::progress_last < $cur) ||
       ($cur - $::progress_last) >= (0.05 * $tot) ||
       ($now - $::progress_last_time) >= 5} {
       set ::progress_last $cur
       set ::progress_last_time $now
       set percentage [expr round($cur*100/$tot)]
       set ticks [expr $percentage/2]
       if {$cur == 0} {
           set eta   ETA:[format %7s Unknown]
       } elseif {$cur >= $tot} {
           set eta   TOT:[format %7d [expr int($now - $::progress_start)]]s
       } else {
           set eta   ETA:[format %7d [expr int(($tot - $cur) * ($now - $::progress_start)/$cur)]]s
       }
       set lticks [expr 50 - $ticks]
       set str "Reading File: $file |[string repeat = $ticks]"
       append str "[string repeat . $lticks]| [format %3d $percentage]% |$eta\r"
       puts -nonewline stdout $str
       if {$cur >= $tot} {
           puts ""
       }
       flush stdout
   }
}


proc readONIOM {readFile} {
# read coordinates, layer info, and redundant and return them
    
        set x "" 
        set y ""
        set z ""
        set layer ""
        set redundant 0
        set x [lindex $readFile 2]
        set y [lindex $readFile 3]
        set z [lindex $readFile 4]
        set redundant [lindex $readFile 1]
        set pos [expr [string last $z $readFile] +[string length $z]]
        set layer1 [string range $readFile $pos [string length $readFile]]
        set pos [string last [lindex $layer1 0] $layer1]
        set layer "[string range $layer1 $pos [string length $readFile]]"

	# read Atoms Data
        set pdbName "H"
        set resName "ABC"
        set resid "9999"
        set atomType "H"
        set charge "0.000"        
        
        
        if {[string first "(" $readFile]!=-1} {
            # 1- remove () characters
            set readFile [lindex $readFile 0]
            set readFileMod [string map {( " "}  $readFile]
            set readFileMod [string map {) " "}  $readFileMod]

            # 2- remove item:  , and =
            set readFileMod [string map {, " "}  $readFileMod]
            set readFileMod [string map {= " "}  $readFileMod]

            set pdbName [lindex $readFileMod 2]
            set resName [lindex $readFileMod 4]
            set resid   [lindex $readFileMod 6]

            # 3- get charge & atomtype
            if {[string match "*--*" $readFileMod]==1} {
        	    set charge [string range [lindex $readFileMod 0] [expr [string length [lindex $readFileMod 0]]-9] [string length [lindex $readFileMod 0]]]
                set atomType [string range [lindex $readFileMod 0] 0 [expr [string length [lindex $readFileMod 0]]-11] ]

            } else {

                set charge [string range [lindex $readFileMod 0] [expr [string length [lindex $readFileMod 0]]-8] [string length [lindex $readFileMod 0]]]
                set atomType [string range [lindex $readFileMod 0] 0 [expr [string length [lindex $readFileMod 0]]-10] ]
            }

        } else {
            #There is no PDB information
            # 1- get charge & atomtype
            set readFileMod [lindex $readFile 0]
            if {[string match "*--*" $readFileMod]==1} {
        	    set charge [string range [lindex $readFileMod 0] [expr [string length [lindex $readFileMod 0]]-9] [string length [lindex $readFileMod 0]]]
                set atomType [string range [lindex $readFileMod 0] 0 [expr [string length [lindex $readFileMod 0]]-11] ]

            } else {
                set charge [string range [lindex $readFileMod 0] [expr [string length [lindex $readFileMod 0]]-8] [string length [lindex $readFileMod 0]]]
                set atomType [string range [lindex $readFileMod 0] 0 [expr [string length [lindex $readFileMod 0]]-10] ]
            }

        }

	return "$charge $pdbName $atomType $resName $resid $x $y $z $redundant {$layer}"
}


proc readFile {File} {
    global loadFile charges input keywords

    ### Open file
    set loadFile [open $File r]

    ### ProgressBar
    progress_init [file size $File]

    ###a read File
    while {![eof $loadFile]} {
	    set readFile [gets $loadFile]

        ## READ INPUT KEYWORDS
        if {[string first "Gaussian 09:" $readFile]!=-1} {
            readKeywords $readFile
        }


        ## READ INPUTFILE COORDINATES
        if {$readFile==" Symbolic Z-matrix:"} {
            readInput
        }


        ## READING Mulliken Charges
        if {[string first "Mulliken charges and spin densities:" $readFile]!=-1 || \
            [string first "Mulliken atomic charges:" $readFile]!=-1 || \
            [string first "Mulliken charges:" $readFile]!=-1} {

            readMullikenCharges 
        }


        ## READING GEOMETRIES
        if {[string first "Z-Matrix orientation:" $readFile]!=-1 || \
            [string first "Input orientation:" $readFile]!=-1  || \
            [string first "Standard orientation:" $readFile]!=-1  } {

            readGeometries 
        }

    }

    ### ProgressBar
    progress_tick  $File [file size $File]
    close $loadFile
}


proc readKeywords {readFile} {
    global keywords loadFile File

    set keywords ""

    # ProgressBar
    progress_tick $File [tell $loadFile]

    set readFile [gets $loadFile]
    set end 0

    while {$end==0} {

        if {[string first "%" $readFile]!=-1} {
            set readFile [string range $readFile 1 [string length $readFile]]
            set keywords "$keywords$readFile \n"
        }

        if {[string first "#" $readFile]!=-1} {
            while {[string first "--" $readFile]==-1} {
                set readFile [string range $readFile 1 [string length $readFile]]
                set keywords "$keywords$readFile"
                set readFile [gets $loadFile]
            }
            set end 1
        }
        set readFile [gets $loadFile]
    }
}

proc readInput {} {
    global loadFile File charges input keywords 
    
    # ProgressBar
    progress_tick $File [tell $loadFile]

    # Obtain charges
    set charges ""
    set readFile [gets $loadFile]

    while {[lindex $readFile 0]=="Charge"} {
        set readFile [gets $loadFile]
        regsub -all "=" $readFile " " ReadFile
        set charges "[lappend charges [lindex $readFile 1] [lindex $readFile 3]] "
    }

    # Obtain atoms Data and coordinates
    set atomCount 1 
    set readFile [gets $loadFile]
    while {$readFile!=" "} {
        set data [readONIOM $readFile]
        dict set input $atomCount charge [lindex $data 0]
        dict set input $atomCount pdbName [lindex $data 1]
        dict set input $atomCount atomType [lindex $data 2]
        dict set input $atomCount resName [lindex $data 3]
        dict set input $atomCount resid [lindex $data 4]
        dict set input $atomCount x [lindex $data 5]
        dict set input $atomCount y [lindex $data 6]
        dict set input $atomCount z [lindex $data 7]
        dict set input $atomCount redundant [lindex $data 8]
        dict set input $atomCount layer [lindex $data 9]
        incr atomCount
        set readFile [gets $loadFile]
    }

}

proc readMullikenCharges {} {
    global loadFile input Mcharge
    set readFile [gets $loadFile]
    for {set i 1} {$i<=[dict size $input]} {incr i} {
        set readFile [gets $loadFile]
        dict set Mcharge $i charge [lindex $readFile 2]
    }
}

proc readGeometries {} {
    global loadFile File standardCount xyz input xyzCount

    # Progress Bar
    progress_tick  $File [tell $loadFile]

    # increment geometry number
    incr standardCount

    # read header
    for {set a 0} {$a<=4} {incr a} {set readFile [gets $loadFile]}
    
    # read Geometry
    set xyzCount 1

    dict for {id info} $input {
        dict with info {
                 dict set xyz $standardCount-$xyzCount x [lindex $readFile 3]
                 dict set xyz $standardCount-$xyzCount y [lindex $readFile 4]
                 dict set xyz $standardCount-$xyzCount z [lindex $readFile 5]
                 dict set xyz $standardCount-$xyzCount charge $charge
                 dict set xyz $standardCount-$xyzCount pdbName $pdbName
                 dict set xyz $standardCount-$xyzCount atomType $atomType
                 dict set xyz $standardCount-$xyzCount resName $resName
                 dict set xyz $standardCount-$xyzCount resid $resid
                 dict set xyz $standardCount-$xyzCount redundant $redundant
                 dict set xyz $standardCount-$xyzCount layer $layer
        }
    incr xyzCount
    set readFile [gets $loadFile]
    }
}

proc savePDBGeometries {fileName} {

    global xyz xyzCount count

    set savePDB [open $fileName-all-Geometries.pdb w]
    #puts "All Geometries in PDB format                : $fileName-all-Geometries.pdb  |"; set countC 1

    # Calculate the number of Steps
    set optSteps  [expr 1 + ([dict size $xyz]/$xyzCount)]

    # Print geometries
    set count 1
    #set countC 0
    puts $savePDB "HEADER"
    dict for {id1 info} $xyz {

     dict with info {
        set idnew [string map {- " "} $id1]

        # write text for counting
        #if {$count!=[lindex $idnew 0]} {incr countC; if {$countC==10} {puts -nonewline "#";set countC 0} else { puts -nonewline $countC; flush stdout}} 
        if {$count!=[lindex $idnew 0]} {puts -nonewline stdout "All Geometries in PDB format                : $fileName-all-Geometries.pdb  [format %-8s ($count/[expr $optSteps-1])]\r";  flush stdout}

        ## Create PDB format
        if {$count!=[lindex $idnew 0]} {puts $savePDB "END"; puts $savePDB "HEADER"; puts $savePDB "Standard Geometry [lindex $idnew 0]"; set count [lindex $idnew 0]}
        set chain "X"
        if {[lindex $layer 0]=="L"} {set chain L} elseif {[lindex $layer 0]=="H"} {set chain H} else {set chain "M"}
        set atomType [lindex [string map {- " "}  $atomType] 0]
        puts $savePDB "[format %-4s "ATOM"] [format %6s [lindex $idnew 1]] [format %-4s $pdbName][format %4s $resName] [format %-1s $chain] [format %-7s $resid] [format %7.3f $x] [format %7.3f $y] [format %7.3f $z] [format %5s "1.00"] [format %-8s "00.00"] [format %8s $atomType]"

     }
    }

    puts $savePDB "END"
    close $savePDB
}


proc infoCharges {} {

    global xyz xyzCount 

    # Calculate the number of Steps
    set optSteps  [expr 1 + ([dict size $xyz]/$xyzCount)]

    # Gather residues information
    set resIdList ""
    set resNameList ""
    set resNameListHL ""
    set resNameListLL ""
    set chargePositiveResidues "ARG LYS HIP"
    set chargeNegativeResidues "GLU ASP"
    set neutralResidues "PRO GLN SER TRP HOH CYS TYR ILE ASN VAL GLY HIS PHE ALA LEU THR MET"
    set unknownResidues ""
    set chargeUnknown 0
    set chargeNegative 0
    set chargePositive 0
    set chargeUnknownHL 0
    set chargeNegativeHL 0
    set chargePositiveHL 0

    set resUnknown 0
    set resNegative 0
    set resPositive 0
    set resUnknownHL 0
    set resNegativeHL 0
    set resPositiveHL 0

    for {set a 1} {$a<=[expr $xyzCount -1]} {incr a} {

        set pdbName [dict get $xyz $optSteps-$a pdbName]
        set resName [dict get $xyz $optSteps-$a resName]
        set resid [dict get $xyz $optSteps-$a resid]
        set layer [dict get $xyz $optSteps-$a layer]
        set charge [dict get $xyz $optSteps-$a charge]

        # make residues list
        if {[lsearch -exact $resIdList $resid]==-1} {
            set resNameList [lappend resNameList $resName]
            set resIdList   [lappend resIdList $resid]

            if {$layer=="L" || $layer=="l"} {set resNameListLL [lappend resNameListLL $resName]}
            if {$layer=="H" || $layer=="h"} {set resNameListHL [lappend resNameListHL $resName]}

            # Calculate charges
            if       {[lsearch $chargePositiveResidues $resName]!=-1} {set chargePositive [expr $chargePositive + $charge]; incr resPositive
            } elseif {[lsearch $chargeNegativeResidues $resName]!=-1} {set chargeNegative [expr $chargeNegative + $charge]; incr resNegative
            } elseif {[lsearch $neutralResidues $resName]!=-1       } {
            } else   {set unknownResidues "$unknownResidues $resName"; set chargeUnknown [expr $chargeUnknown + $charge]; incr resUnknown}


            # charges HL
            if {[lindex $layer 0]=="H"} {
                puts "ok"
                if {[lsearch $chargePositiveResidues $resName]!=-1} {set chargePositiveHL [expr $chargePositiveHL + $charge]; incr resPositiveHL
                } elseif {[lsearch $chargeNegativeResidues $resName]!=-1} {set chargeNegativeHL [expr $chargeNegativeHL + $charge]; incr resNegativeHL
                } else   {set chargeUnknownHL [expr $chargeUnknownHL + $charge]; incr resUnknownHL  }

            }




        }


    }

    puts "Number of Residues : [expr $resPositive +  $resNegative + $resUnknown ] ([expr $chargePositive +  $chargeNegative + $chargeUnknown ])" 
    puts "Positive Residues  : $resPositive ($chargePositive).  High Level : $resPositiveHL ($chargePositiveHL) ($chargePositiveResidues)"
    puts "Negative Residues  : $resNegative ($chargeNegative).  High Level : $resNegativeHL ($chargeNegativeHL) ($chargeNegativeResidues)"
    puts "Unknown  Residues  : $resUnknown  ($chargeUnknown).   High Level : $resUnknownHL  ($chargeUnknownHL) ($unknownResidues)"

    puts "\n"

}

proc saveInputLastGeometryAll {fileName} {

    global xyz xyzCount charges count Mcharge keywords

    set newCOM [open $fileName-input-final-all.com w]
    puts -nonewline "\nInput File with final geometry (all)        : $fileName-input-final-all.com"

    # Puts initial Data
     puts $newCOM $keywords
     puts $newCOM ""
     puts $newCOM "Title"
     puts $newCOM ""
     puts $newCOM  "$charges"

    # Calculate the number of Steps
    set optSteps  [expr 1 + ([dict size $xyz]/$xyzCount)]

    # Print geometries
    for {set a 1} {$a<=[expr $xyzCount -1]} {incr a} {

        set x [dict get $xyz $optSteps-$a x]
        set y [dict get $xyz $optSteps-$a y]
        set z [dict get $xyz $optSteps-$a z]
        set pdbName [dict get $xyz $optSteps-$a pdbName]
        set resName [dict get $xyz $optSteps-$a resName]
        set resid [dict get $xyz $optSteps-$a resid]
        set atomType [dict get $xyz $optSteps-$a atomType]
        set layer [dict get $xyz $optSteps-$a layer]
        set charge [dict get $xyz $optSteps-$a charge]
        set redundant [dict get $xyz $optSteps-$a redundant]

        ## Create PDB format
        set chain "X"
        if {[lindex $layer 0]=="L"} {set chain "L"} else {set chain "H"}

        ## Only changes the charge of the high level
        if {$chain =="H"} {set charge [dict get $Mcharge $a charge]}

        if {$count==$optSteps} {
            set atomData "$atomType-$charge\(PDBName=$pdbName,ResName=$resName,ResNum=$resid\)"
            set coord "[format "%10s"  [format "%5.5f" $x]] [format "%10s" [format "%5.5f" $y]] [format "%10s"  [format "%5.5f" $z]]"
            puts $newCOM " [format "%-60s" $atomData] [format "%-4s" $redundant] $coord $layer"
        }

    }
    close $newCOM
}



proc saveInputLastGeometryHL {fileName} {

    global xyz xyzCount charges count Mcharge keywords

    
    set newHCOM [open $fileName-input-final-highlevel.com w]
    set count 0
    puts -nonewline "\nInput File with final geometry (high-level) : $fileName-input-final-highlevel.com"

    # Puts initial Data
    puts $newHCOM "$keywords"
    puts $newHCOM ""
    puts $newHCOM "Title"
    puts $newHCOM ""
    puts $newHCOM "[lindex $charges [expr [llength $charges] -2]] [lindex $charges [expr [llength $charges] -1]]"

    # Calculate the number of Steps
    set optSteps  [expr 1 + ([dict size $xyz]/$xyzCount)]

    # Print geometries

    for {set a 1} {$a<=[expr $xyzCount-1]} {incr a} {
        set x [dict get $xyz $optSteps-$a x]
        set y [dict get $xyz $optSteps-$a y]
        set z [dict get $xyz $optSteps-$a z]
        set pdbName [dict get $xyz $optSteps-$a pdbName]
        set resName [dict get $xyz $optSteps-$a resName]
        set resid [dict get $xyz $optSteps-$a resid]
        set atomType [dict get $xyz $optSteps-$a atomType]
        set layer [dict get $xyz $optSteps-$a layer]

        if {[lindex $layer 0]=="H" || ([lindex $layer 0]=="L" && [llength $layer]!=1 )} {

            #change coordinates and atomtype of link atoms
            if { [llength $layer]!=1} {
                set pdbName [lindex $layer 1]; set atomType $pdbName

                # search for the atom that is linked to it
                set linkatom [lindex $layer 2]
                set xlink0 [dict get $xyz $optSteps-$linkatom x]
                set ylink0 [dict get $xyz $optSteps-$linkatom y]
                set zlink0 [dict get $xyz $optSteps-$linkatom z]

                set xatomnew [expr $x - $xlink0]
                set yatomnew [expr $y - $ylink0]
                set zatomnew [expr $z - $zlink0]

                # Convert points0 to Polar coordinates
                set coord [rect_to_spherical $xatomnew $yatomnew $zatomnew]
                set rad   [lindex $coord 0]
                set phi   [lindex $coord 1]
                set theta [lindex $coord 2]

                #new rad
                set rad 1.012

                set newpoint [spherical_to_rect $rad $phi $theta]
                set x [expr [lindex $newpoint 0] + $xlink0]
                set y [expr [lindex $newpoint 1] + $ylink0]
                set z [expr [lindex $newpoint 2] + $zlink0]

            }

            set atomType [lindex [string map {- " "}  $atomType] 0]
            puts $newHCOM "[format %-6s $atomType] [format %7.3f $x] [format %7.3f $y] [format %7.3f $z]"
        }

    }
    puts $newHCOM ""
    close $newHCOM

}

proc saveInputFile {fileName} {

    global input keywords charges

    puts -nonewline "\nOriginal Input File                         : $fileName-input-inital-all.com"

    set saveCOM [open $fileName-input-initial-all.com w]

    # Puts initial Data
    puts $saveCOM $keywords
    puts $saveCOM ""
    puts $saveCOM "Title"
    puts $saveCOM ""
    puts $saveCOM $charges

    dict for {id info} $input {
        dict with info {
         set atomData "$atomType-$charge\(PDBName=$pdbName,ResName=$resName,ResNum=$resid\)"
         set coord "[format "%10s"  [format "%5.5f" $x]] [format "%10s" [format "%5.5f" $y]] [format "%10s"  [format "%5.5f" $z]]"
         puts $saveCOM " [format "%-60s" $atomData] [format "%-4s" $redundant] $coord $layer"
        }
    }
    puts $saveCOM ""
    close $saveCOM
}






################ START


set version "1.1 (2/7(2014)"

#### HEADER
puts "   ____  _   _ _____ ____  __  __             _                  _             "
puts "  / __ \\| \\ | |_   _/ __ \\|  \\/  |           | |                | |            "
puts " | |  | |  \\| | | || |  | | \\  / |   _____  _| |_ _ __ __ _  ___| |_ ___  _ __ "
puts " | |  | | . ` | | || |  | | |\\/| |  / _ \\ \\/ / __| '__/ _` |/ __| __/ _ \\| '__|"
puts " | |__| | |\\  |_| || |__| | |  | | |  __/>  <| |_| | | (_| | (__| || (_) | |   "
puts "  \\____/|_| \\_|_____\\____/|_|  |_|  \\___/_/\\_\\ __|_|  \\__,_|\\___|\\__\\___/|_|    Version: $version  \n"
       

#### READ Data from Shell 

set File     [lindex $argv 0]
if {[lindex $argv 0]== "--help" ||  $File==""} {
        puts  " Usage :  tclsh extract.tcl  \[a\] \n"
        puts  "           \[a\] : Gaussian Output File (*.log)"
        puts "\n"
        puts " Developers: Nuno Sousa Cerqueira (nscerque@fc.up.pt)"
        puts "             Erica Cristina Morena (ericamoreno@unb.br)"
	    puts "\n"
        exit
}


    set standardCount 0 ;# number of standard orientations


#### READFILE
readFile $File


#### INFO FILES
#puts "\nInfo:\n"
#infoCharges

#### MAKE Files
puts "Creating Files:\n"
set fileName [file rootname $File]


#### PRINT GEOMETRIES IN PDB FORMAT
savePDBGeometries $fileName


#### PRINT INPUT WITH LAST GEOMETRY - ALL
saveInputLastGeometryAll $fileName


#### Print INPUT LAST GEOMETRY - HIGH LEVEL
saveInputLastGeometryHL $fileName


#### PRINT INPUT FILE
saveInputFile $fileName


puts "\n\nDone!"



