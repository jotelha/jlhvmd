#!/usr/bin/tclsh
# JlhVmd, a VMD package to manipulate interfacial systems or other
# topology related properties in VMD with the help of TopoTools and PbcTools
#
# Copyright (c) 2018,2019
#               by Johannes Hoermann <johannes.hoermann@imtek.uni-freiburg.de>
#
# $Id: jlhvmd.tcl,v 0.2 2019/05/16 $
#
# Sample usage for indenter insertion:
#
#   vmd> package requie jlhvmd
#   vmd> jlh set distance 30.0 indenterInfile indenter.lammps interfaceInfile interface.lammps outputPrefix system
#   vmd> jlh use sds
#   vmd> jlh read bb bb.yaml
#   vmd> jlh insert
#
# Sample usage for wrapping periodic images and joining split residues:
#
#   vmd> package require jlhvmd
#   vmd> jlh set interfaceInfile initial_config.lammps outPrefix default
#   vmd> jlh use sds
#   vmd> jlh read bb bb.yaml
#   vmd> jlh init
#   vmd> jlh show surfactant
#   vmd> jlh wrap atom
#   vmd> jlh join residue
#   vmd> jlh write
#
# will use parameters for SDS, merge an interfacial system's 'interface.lammps'
# data file with an indenter's 'indenter.pdb', remove any overlap and write
# system.lammps, system.psf, system.pdb data files as well as a tga snapshot
# system.tga of the resulting system
#
# No routine for displacing charged ions implemented yed. THEY ARE REMOVED.
# Double-check for section
#   Info) Identify ovelap...
#   Info) #atoms in overlapping SOD:                0
#   Info) #atoms in overlapping TIP3:            9354
#   Info) #atoms in overlapping SDS:                0
# in output to make sure only water has been removed.
#
# CHANGELOG
#
# ## [0.2] 2019-05-14
# ### Changed
# - Suppose to have substrate terminated at z = 0
# - Package structure copied from Axel Kohlmeyer's topotools
#
# ## [0.2.1] 2019-05-19
# ### Changed
# - jlh user interface
#
# ## [0.2.2] 2019-05-20
# ### Changed
# - random position changed from small volume between indenter and substrate to
#   whole bounding box except lower substrate part
#
# ## [0.2.3] 2019-08-10
# ### Changed
# - number of digits to fill in with type names now determined by counting
#   characters of largest type name after integer sorting instead of
#   using the actual number of types.
#
# ## [0.2.4] 2019-08-11
# ### Changed
# - commands case-insenitive, allow i.e. "jlh use SDS" as well as "jlh use sds"
# - removed obsolete pdb processing functionality and other code snippets
#
# ## [0.2.5] 2019-08-22
# ### Changed
# - added bb_center as global variable
# - split bb positioning and wrapping
# - added several wrapping and joining commands

namespace eval ::JlhVmd:: {
    variable version 0.2.5

    package require topotools
    package require pbctools
    package require yaml

    # default values:

    # shift tip as to z-center apex 12 nm above substrate
    set desired_distance 120.0

    # distance within which molecules are regarded as overlapping
    set overlap_distance 2.0
    # clip indenter that far away from cell boundary:
    set padding_distance 5.0

    set bounding_box { {0. 0. 0.} {150. 150. 150.} }
    set bb_center    { 75. 75. 75. }
    # io
    set interface_infile "interface.lammps"
    set indenter_infile  "indenter.lammps"
    set out_prefix "system"

    # number of distinct atomic layers within AU substrate building block cell
    set sb_cell_inner_lattice_points { 2.0 6.0 3.0}
    # multiples of cell building block in substrate
    set sb_cell_multiples {51.0 30.0 2.0}

    set system_rep 0
    set indenter_rep 0
}

# help/usage/error message and online documentation.
proc ::JlhVmd::usage {} {
    vmdcon -info ""
    vmdcon -info "JlhVmd, a VMD package to manipulate interfacial systems or other"
    vmdcon -info "topology related properties in VMD with the help of TopoTools and PbcTools."
    vmdcon -info ""
    vmdcon -info ""
    vmdcon -info "usage: jlh <command> \[args...\] <flags>"
    vmdcon -info ""
    vmdcon -info ""
    vmdcon -info "common flags (not implemented):"
    vmdcon -info ""
    vmdcon -info "  -molid     <num>|top    molecule id (default: 'top')"
    vmdcon -info "  -sel       <selection>  atom selection function or text (default: 'all')"
    vmdcon -info ""
    vmdcon -info ""
    vmdcon -info "commands:"
    vmdcon -info ""
    vmdcon -info "  help                                          prints this message"
    vmdcon -info ""
    vmdcon -info "  read <key> <file> \[ <key> <file> \[ ... \] \]    reads parameter <key> from <file>."
    vmdcon -info "  set <key> <value> \[ <key> <value> \[ ... \] \]   sets parameter <key> to <value>."
    vmdcon -info "  use <key>                                     uses standard settings for surfactant <key>."
    vmdcon -info ""
    vmdcon -info "  info                                          display system information."
    vmdcon -info "  init                                          initializes system without manipulation."
    vmdcon -info "  insert                                        inserts indenter into interfacial system."
    vmdcon -info "  join <key>                                    (re-)joins residues split across boundries."
    vmdcon -info "  render <key>                                  render image of scene to .tga file."
    vmdcon -info "  wrap <key>                                    wrap system into one periodic image."
    vmdcon -info "  write                                         write .lammps, .psf and .pdb output files."
    vmdcon -info ""
    vmdcon -info "key - value pairs"
    vmdcon -info ""
    vmdcon -info "  join:   residue           "
    vmdcon -info "  read:   bb                bounding box, expects yaml file with keys 'xlo', 'xhi', 'ylo', 'yhi', 'zlo', 'zhi'."
    vmdcon -info "  render: nonsolvent        "
    vmdcon -info "          solvent           "
    vmdcon -info "          surfactant        "
    vmdcon -info "  set:    bb                bounding box ({ { float float float } { float float float } })"
    vmdcon -info "          distance          desired distance between surface and indenter (float)"
    vmdcon -info "          interfaceInfile   input LAMMPS data file of interface (str)"
    vmdcon -info "          indenterInfile    input LAMMPS data file of indenter (str)"
    vmdcon -info "          outputPrefix      output prefix prepended to all resulting files (str)"
    vmdcon -info "  show:   nonsolvent        "
    vmdcon -info "          solvent           "
    vmdcon -info "          surfactant        "
    vmdcon -info "  use:    sds               "
    vmdcon -info "          ctab              "
    vmdcon -info "  wrap:   atom              "
    vmdcon -info "          residue           "
    vmdcon -info ""
    vmdcon -info ""
    vmdcon -info "sample usage:"
    vmdcon -info ""
    vmdcon -info "  vmd> package requie jlhvmd"
    vmdcon -info "  vmd> jlh set distance 30.0 indenterInfile indenter.lammps interfaceInfile interface.lammps outputPrefix system"
    vmdcon -info "  vmd> jlh use sds"
    vmdcon -info "  vmd> jlh read bb bb.yaml"
    vmdcon -info "  vmd> jlh insert"
    vmdcon -info ""
    vmdcon -info "will use parameters for SDS, merge an interfacial system's 'interface.lammps'"
    vmdcon -info "data file with an indenter's 'indenter.pdb', remove any overlap and write"
    vmdcon -info "system.lammps, system.psf, system.pdb data files as well as a tga snapshot"
    vmdcon -info "system.tga of the resulting system."
    vmdcon -info ""
    vmdcon -info ""
    vmdcon -info "Copyright (c) 2018,2019"
    vmdcon -info "              by Johannes Hoermann <johannes.hoermann@imtek.uni-freiburg.de>"
    vmdcon -info ""
    return
}

# the main frontend command.
# this takes care of all sanity checks on arguments and
# then dispatches the subcommands to the corresponding
# subroutines.
proc ::JlhVmd::jlh { args } {
    variable version

    variable desired_distance
    variable overlap_distance
    variable padding_distance
    variable bounding_box

    variable interface_infile
    variable indenter_infile
    variable out_prefix

    set cmd {}
    # set sel {} ; # need to initialize it here for scoping

    # process generic arguments and remove them
    # from argument list.
    set newargs {}
    for {set i 0} {$i < [llength $args]} {incr i} {
        set arg [lindex $args $i]

        if {[string match -?* $arg]} {

            set val [lindex $args [expr $i+1]]

            switch -- $arg {
                -molid {
                    if {[catch {molinfo $val get name} res]} {
                        vmdcon -err "Invalid -molid argument '$val': $res"
                        return
                    }
                    set molid $val
                    if {[string equal $molid "top"]} {
                        set molid [molinfo top]
                    }
                    incr i
                }

                -sel {
                    # check if the argument to -sel is a valid atomselect command
                    if {([info commands $val] != "") && ([string equal -length 10 $val atomselect])} {
                        set localsel 0
                        set selmol [$val molid]
                        set sel $val
                    } else {
                        set localsel 1
                        set seltxt $val
                    }
                    incr i
                }

                -- break

                default {
                    vmdcon -info "default: $arg"
                }
            }
        } else {
            lappend newargs $arg
        }
    }

    set retval ""
    if {[llength $newargs] > 0} {
        set cmd [lindex $newargs 0]
        set newargs [lrange $newargs 1 end]
    } else {
        set newargs {}
        set cmd help
    }

    # check whether we have a valid command.
    set validcmd {
      "set" "read" "use" "init"
      "show" "render"
      "wrap" "join"
      "insert" "write"
      "help" "info" }
    if {[lsearch -exact $validcmd $cmd] < 0} {
        vmdcon -err "Unknown sub-command '$cmd'"
        usage
        return
    }

    # branch out to the various subcommands
    switch -nocase -- $cmd {
        "set" {
            while {[llength $newargs] > 1} {
                set key [lindex $newargs 0]
                set newargs [lrange $newargs 1 end]
                switch -nocase -- $key {
                    distance {
                        set desired_distance [lindex $newargs 0]
                        set newargs [lrange $newargs 1 end]
                    }

                    bb {
                        set bounding_box [lindex $newargs 0]
                        compute_bb_center
                        set newargs [lrange $newargs 1 end]
                    }

                    indenterInfile {
                        set indenter_infile [lindex $newargs 0]
                        set newargs [lrange $newargs 1 end]
                    }

                    interfaceInfile {
                        set interface_infile [lindex $newargs 0]
                        set newargs [lrange $newargs 1 end]
                    }

                    outputPrefix {
                        set out_prefix [lindex $newargs 0]
                        set newargs [lrange $newargs 1 end]
                    }

                    default {
                        vmdcon -warn "Unknown parameter: $key"
                    }
                }
            }
            set retval 0
        }

        "read" {
            set key [lindex $newargs 0]
            set newargs [lrange $newargs 1 end]
            switch -nocase -- $key {
                bb {
                    set bb_file [lindex $newargs 0]
                    set newargs [lrange $newargs 1 end]
                    set retval [read_bb_from_yaml $bb_file]
                }

                default {
                    vmdcon -err "Unknown parameter: $key"
                    set retval 1
                }
            }
        }

        "render" {
            set key [lindex $newargs 0]
            set newargs [lrange $newargs 1 end]
            switch -nocase -- $key {
                solvent {
                    render_solvent
                    set retval 0
                }

                surfactant {
                    render_surfactant
                    set retval 0
                }

                nonsolvent {
                    render_nonsolvent
                    set retval 0
                }

                default {
                    vmdcon -err "Unknown parameter: $key"
                    set retval 1
                }
            }
        }

        use {
            set key [lindex $newargs 0]
            set newargs [lrange $newargs 1 end]
            switch -nocase -- $key {
                ctab {
                    use_CTAB
                    set retval 0
                }

                sds {
                    use_SDS
                    set retval 0
                }

                default {
                    vmdcon -err "Unknown parameter: $key"
                    set retval 1
                }
            }
        }

        "show" {
            set key [lindex $newargs 0]
            set newargs [lrange $newargs 1 end]
            switch -nocase -- $key {
                solvent {
                    show_solvent_only
                    set retval 0
                }

                surfactant {
                    show_surfactant_only
                    set retval 0
                }

                nonsolvent {
                    show_nonsolvent
                    set retval 0
                }

                default {
                    vmdcon -err "Unknown parameter: $key"
                    set retval 1
                }
            }
        }

        "wrap" {
            set key [lindex $newargs 0]
            set newargs [lrange $newargs 1 end]
            switch -nocase -- $key {
                residue {
                    position_bb
                    wrap_residue_into_bb
                    set retval 0
                }

                atom {
                    position_bb
                    wrap_atom_into_bb
                    set retval 0
                }

                default {
                    vmdcon -err "Unknown parameter: $key"
                    set retval 1
                }
            }
        }

        "join" {
            set key [lindex $newargs 0]
            set newargs [lrange $newargs 1 end]
            switch -nocase -- $key {
                residue {
                    position_bb
                    join_residue
                    set retval 0
                }

                default {
                    vmdcon -err "Unknown parameter: $key"
                    set retval 1
                }
            }
        }

        "init" {
            initialize $interface_infile
            set retval 0
        }

        insert {
            batch_process_lmp_visual $interface_infile $indenter_infile $out_prefix
            set retval 0
        }

        "write" {
            write_top_all $out_prefix
            set retval 0
        }

        "info" {
            display_system_information
            set retval 0
        }

        "help" {
            usage
            set retval 0
        }

        default {
            vmdcon -err "Unknown sub-command: $cmd"
            set retval 1
        }
    }
    return $retval
}

# two small supportive functions to enhance vmd's matrix manipulation toolset
proc ::JlhVmd::scalar_times_matrix {s m} {
    set res {}
    foreach row $m {
        lappend res [vecscale $s $row ]
    }
    return $res
}

# returns a scaling transformation matrix
proc ::JlhVmd::transscale {s} {
  set res  [ scalar_times_matrix $s [transidentity] ]
  lset res {3 3} 1.0
  return $res
}

proc ::JlhVmd::use_CTAB {} {
  variable surfactant_resname CTAB
  variable counterion_name    BR
  variable counterion_resname BR
  variable solvent_resname    TIP3
  variable substrate_name     AU
  variable substrate_resname  AUM

  variable H2O_H_type 8
  variable H2O_O_type 9

  # suggestion from https://lammps.sandia.gov/threads/msg21297.html
  variable type_name_list { \
      1 HL  \
      2 HAL2  \
      3 HAL3  \
      4 CTL2  \
      5 CTL3  \
      6 CTL5  \
      7 NTL  \
      8 HT  \
      9 OT  \
     10 BR  \
     11 AU  \
  }
}

proc ::JlhVmd::use_SDS {} {
  variable surfactant_resname SDS
  variable counterion_name    SOD
  variable counterion_resname SOD
  variable solvent_resname    TIP3
  variable substrate_name     AU
  variable substrate_resname  AUM

  variable H2O_H_type 8
  variable H2O_O_type 9

  # from SDS-related data file
  variable type_name_list { \
      1 HAL2  \
      2 HAL3  \
      3 CTL2  \
      4 CTL3  \
      5 OSL  \
      6 O2L  \
      7 SL  \
      8 HT  \
      9 OT  \
     10 SOD  \
     11 AU  \
  }
}

# ##################
# au_cell_P1_111.pdb
# ##################
# extremes (~ measure minmax)
#   0.28837   0.49948   0.70637  nm
#   2.884     4.995     7.064    Angstrom
# ##############################################################################
# TITLE     Periodic slab: SURF, t= 0.0
# REMARK    THIS IS A SIMULATION BOX
# CRYST1    2.884    4.995    7.064  90.00  90.00  90.00 P 1           1
# MODEL        1
# ATOM	  1  Au  SURF    1	 0.000   0.000   0.000  1.00  0.00
# ATOM	  2  Au  SURF    1	 1.440   2.500   0.000  1.00  0.00
# ATOM	  3  Au  SURF    1	 0.000   1.660   2.350  1.00  0.00
# ATOM	  4  Au  SURF    1	 1.440   4.160   2.350  1.00  0.00
# ATOM	  5  Au  SURF    1	 1.440   0.830   4.710  1.00  0.00
# ATOM	  6  Au  SURF    1	 0.000   3.330   4.710  1.00  0.00
# TER
# ENDMDL
# ##############################################################################

proc ::JlhVmd::init_system { infile { psffile "" } } {
  variable system_id
  variable system
  variable type_name_list

  if { $psffile ne "" } {
    set system_id [mol new $psffile type psf waitfor all]
    #mol addfile XYZ.dcd type dcd first 999 last 999 waitfor all molid $mol
    topo readlammpsdata $infile full -molid $system_id
  } else {
    # no psf topology, use topotools to derive types
    set system_id [topo readlammpsdata $infile full]

    # https://sites.google.com/site/akohlmey/software/topotools/topotools-tutorial---various-tips-tricks
    topo guessatom element mass
    # topo guessatom name element
    topo guessatom radius element

    # suggestion from https://lammps.sandia.gov/threads/msg21297.html
    foreach {type name} $type_name_list {
      set sel [atomselect $system_id "type '$type'"]
      $sel set name $name
      $sel delete
    }
  }

  set system [atomselect $system_id all]
  $system global

  mol rename $system_id interface
}

proc ::JlhVmd::make_types_ascii_sortable {} {
  # preserve ordering of types when writing output, as TopoTools 1.7
  # sorts types alphabeticall, not numerically,
  # see topotools/topolammps.tcl::TopoTools::writelammpsmasses, line 900:
  #   set typemap  [lsort -unique -ascii [$sel get type]]

  # number of digits necessary to address all types with decimal numbers
  variable system
  variable H2O_H_type
  variable H2O_O_type

  set num_digits [
    string length [ lindex [ lsort -integer [ topo atomtypenames ] ] end ] ]

  vmdcon -info "Prepending zeros to fill ${num_digits} digits to types."

  proc map {lambda list} {
    #upvar num_digits
    set result {}
    foreach item $list {
        lappend result [apply $lambda $item]
    }
    return $result
  }
  # fill types with leading zeroes if necessary
  $system set type [
    map { x {
      upvar 2 num_digits num_digits
      return [format "%0${num_digits}d" $x]
      } } [ $system get type] ]

  # also set type-dependent variables
  set H2O_H_type [format "%0${num_digits}d" $H2O_H_type]
  set H2O_O_type [format "%0${num_digits}d" $H2O_O_type]
  # the following types reside within TopoTools, thus the retyping procedures
  # are placed within the according namespace
  ::TopoTools::make_bond_types_ascii_sortable $system
  ::TopoTools::make_angle_types_ascii_sortable $system
  ::TopoTools::make_dihedral_types_ascii_sortable $system
  ::TopoTools::make_improper_types_ascii_sortable $system
}

proc ::JlhVmd::display_system_information { {mol_id 0} } {
  vmdcon -info "Number of objects:"
  vmdcon -info "Number of atoms:           [format "% 12d" [topo numatoms -molid ${mol_id} ]]"
  vmdcon -info "Number of bonds:           [format "% 12d" [topo numbonds -molid ${mol_id} ]]"
  vmdcon -info "Number of angles:          [format "% 12d" [topo numangles -molid ${mol_id} ]]"
  vmdcon -info "Number of dihedrals:       [format "% 12d" [topo numdihedrals -molid ${mol_id} ]]"
  vmdcon -info "Number of impropers:       [format "% 12d" [topo numimpropers -molid ${mol_id} ]]"

  vmdcon -info "Number of object types:"
  vmdcon -info "Number of atom types:      [format "% 12d" [topo numatomtypes -molid ${mol_id} ]]"
  vmdcon -info "Number of bond types:      [format "% 12d" [topo numbondtypes -molid ${mol_id} ]]"
  vmdcon -info "Number of angle types:     [format "% 12d" [topo numangletypes -molid ${mol_id} ]]"
  vmdcon -info "Number of dihedral types:  [format "% 12d" [topo numdihedraltypes -molid ${mol_id} ]]"
  vmdcon -info "Number of improper types:  [format "% 12d" [topo numimpropertypes -molid ${mol_id} ]]"

  vmdcon -info "Object type names:"
  vmdcon -info "Atom type names:      [topo atomtypenames -molid ${mol_id} ]"
  vmdcon -info "Bond type names:      [topo bondtypenames -molid ${mol_id} ]"
  vmdcon -info "Angle type names:     [topo angletypenames -molid ${mol_id} ]"
  vmdcon -info "Dihedral type names:  [topo dihedraltypenames -molid ${mol_id} ]"
  vmdcon -info "Improper type names:  [topo impropertypenames -molid ${mol_id} ]"
}

proc ::JlhVmd::populate_selections {} {
  variable counterion_name
  variable counterion_resname
  variable substrate_name
  variable substrate_resname
  variable solvent_resname
  variable surfactant_resname
  variable H2O_H_type
  variable H2O_O_type

  variable system
  variable system_id
  variable counterion
  variable nonsolvent
  variable solvent
  variable substrate
  variable surfactant

  variable sb_cell_inner_lattice_points
  variable sb_cell_multiples

  vmdcon -info "Selecting substrate ..."
  set substrate [atomselect $system_id "name $substrate_name"]
  $substrate global
  $substrate set resname $substrate_resname
  vmdcon -info [format "%-30.30s %12d" "#atoms in $substrate_resname:" [$substrate num]]

  set counterion [atomselect $system_id "name $counterion_name"]
  $counterion global
  $counterion set resname $counterion_resname
  vmdcon -info [format "%-30.30s %12d" "#atoms in $counterion_resname:" [$counterion num]]

  # for types with leading zeroes: single quotation marks necessary, otherwise selection fails
  vmdcon -info "Solvent selection by 'type '$H2O_H_type' '$H2O_O_type''"
  set solvent [atomselect $system_id "type '$H2O_H_type' '$H2O_O_type'"]
  $solvent global
  $solvent set resname $solvent_resname
  vmdcon -info [format "%-30.30s %12d" "#atoms in $solvent_resname:" [$solvent num]]

  set surfactant [atomselect $system_id "not resname $substrate_resname \
    $counterion_resname $solvent_resname"]
  $surfactant global
  $surfactant set resname $surfactant_resname
  vmdcon -info [format "%-30.30s %12d" "#atoms in $surfactant_resname:" [$surfactant num]]

  set nonsolvent [atomselect $system_id "not resname $solvent_resname"]
  $nonsolvent global
  vmdcon -info [format "%-30.30s %12d" "#atoms in nonsolvent:" [$nonsolvent num]]

  # get substrate COM and measures
  set sb_center [measure center $substrate]
  vmdcon -info "substrate center:   [format "%8.4f %8.4f %8.4f" {*}$sb_center]"

  # measure inertia returns COM as a 3-vector in first list entry
  set sb_com    [lindex [measure inertia $substrate] 0]
  vmdcon -info "substrate COM:      [format "%8.4f %8.4f %8.4f" {*}$sb_com]"

  # low z cooridnates in 1st row 3rd entry of measure minmax
  set z_shift   [ expr -1.0 * [lindex [measure minmax $substrate] 0 2] ]
  vmdcon -info "substrate low z:    [format "%8.4f" $z_shift]"

  # calculate desired reference COM of substrate
  # (in case of substrate corner in origin 0 0 0)
  set sb_measures [ vecscale -1.0 [ vecsub {*}[measure minmax $substrate] ] ]
  vmdcon -info "substrate measures: [format "%8.4f %8.4f %8.4f" {*}$sb_measures]"

  set sb_com_reference [ vecscale 0.5 $sb_measures ]
  vmdcon -info "reference COM:      [format "%8.4f %8.4f %8.4f" {*}$sb_com_reference]"
  # pbc get returns list of box descriptions as 6-vectors:
  # 3 measures, 3 angles. Use first entry in list, first 3 indices
  # corresponding to 3-vector of length, width, height
  set cell [ lrange [ lindex [pbc get -molid $system_id] 0 ] 0 2 ]
  vmdcon -info "box:                [format "%8.4f %8.4f %8.4f" {*}$cell]"
  set cell_center [vecscale 0.5 $cell]
  vmdcon -info "box center:         [format "%8.4f %8.4f %8.4f" {*}$cell_center]"

  # discrepancy between box measures and substrate measures
  # should meet about one grid constant in x-y direction due to periodicity
  set sb_spacing [ vecsub $cell $sb_measures ]

  # estimate lattice constants in 2D-periodic directions
  # based on subsrate extremes and box measures
  set sb_lattice_constant {}
  set sb_lattice_constant_reference {}
  foreach i {0 1} {
      lappend sb_lattice_constant [ expr [ lindex $sb_measures $i] / ( [lindex $sb_cell_multiples $i] * [lindex $sb_cell_inner_lattice_points $i] - 1.0 ) ]
      lappend sb_lattice_constant_reference [expr [ lindex $cell $i] / ( [lindex $sb_cell_multiples $i] * [lindex $sb_cell_inner_lattice_points $i] )  ]
  }

  vmdcon -info "################################################################################"
  vmdcon -info "effective lattice constant (substrate extremes as reference): [format "%8.4f %8.4f" {*}$sb_lattice_constant]"
  vmdcon -info "reference lattice constant (periodic box as reference):       [format "%8.4f %8.4f" {*}$sb_lattice_constant_reference]"
  vmdcon -info "difference between substrate measures and box measures:       [format "%8.4f %8.4f" {*}$sb_spacing]"

  vmdcon -info "substrate center: [format "%8.4f %8.4f %8.4f" {*}$sb_center]"
  vmdcon -info "substrate COM:    [format "%8.4f %8.4f %8.4f" {*}$sb_com]"
  vmdcon -info "reference COM:    [format "%8.4f %8.4f %8.4f" {*}$sb_com_reference]"
  vmdcon -info "box center:       [format "%8.4f %8.4f %8.4f" {*}$cell_center]"
  vmdcon -info "################################################################################"
}


# adjust position of bounding box
proc ::JlhVmd::position_bb {} {
  variable bb_center
  pbc box -center origin -shiftcenter $bb_center -on
}

# wraps everything into one periodic image
proc ::JlhVmd::wrap_atom_into_bb {} {
    variable bb_center
    pbc wrap -center origin -shiftcenter $bb_center -nocompound -all -verbose
}

# wraps into one periodic image, but keeps residues connected
proc ::JlhVmd::wrap_residue_into_bb {} {
    variable bb_center
    pbc wrap -center origin -shiftcenter $bb_center -compound residue -all -verbose
}

# tries to join residues split across bb boundaries
proc ::JlhVmd::join_residue {} {
    variable bb_center
    pbc join residue -bondlist -all -verbose
}

proc ::JlhVmd::read_indenter_lmp { infile } {
  variable system_id
  variable system
  variable scale_factor
  variable indenter_id
  variable indenter

  set indenter_id [topo readlammpsdata $infile full]
  set indenter [atomselect $indenter_id all]
  $indenter global
  vmdcon -info [format "%-30.30s %12d" "#atoms in indenter:" [$indenter num]]
}

proc ::JlhVmd::init_indenter {} {
  variable indenter
  variable substrate
  # use same atom type as substrate
  $indenter set name [lindex [$substrate get name] 0]
  $indenter set resname [lindex [$substrate get resname] 0]
  $indenter set type [lindex [$substrate get type] 0]
  $indenter set element [lindex [$substrate get element] 0]
  $indenter set mass [lindex [$substrate get mass] 0]
}

proc ::JlhVmd::scale_indenter {} {
  variable indenter
  variable scale_factor
  $indenter move [transscale $scale_factor]
}

proc ::JlhVmd::position_indenter {} {
  variable system_id
  variable substrate
  variable indenter
  variable desired_distance

  # indenter height (lowest )
  set in_measures [ vecscale -1.0 [ vecsub {*}[measure minmax $indenter] ] ]
  set in_height [lindex $in_measures 2]
  set in_com    [lindex [measure inertia $indenter] 0]

  vmdcon -info "indenter measures:    [format "%8.4f %8.4f %8.4f" {*}$in_measures]"
  vmdcon -info "indenter height:      [format "%8.4f" $in_height]"
  vmdcon -info "indenter COM:         [format "%8.4f %8.4f %8.4f" {*}$in_com]"

  set cell [ lrange [ lindex [pbc get -molid $system_id] 0 ] 0 2 ]
  set cell_center [vecscale 0.5 $cell]

  vmdcon -info "cell measures:        [format "%8.4f %8.4f %8.4f" {*}$cell]"
  vmdcon -info "cell center:          [format "%8.4f %8.4f %8.4f" {*}$cell_center]"

  set sb_height [lindex [measure minmax $substrate] { 1 2 } ]
  # set sb_surface_to_cell_center [ expr [ lindex $cell_center 2] - $sb_height ]
  vmdcon -info "substrate height:          [format "%8.4f" $sb_height]"
  vmdcon -info "cell center:          [format "%8.4f %8.4f %8.4f" {*}$cell_center]"

  set indenter_extents  [ measure minmax $indenter ]
  vmdcon -info [ format "      indenter extents: %26s; %26s" \
    [ format "%8.4f %8.4f %8.4f" {*}[ lindex $indenter_extents 0 ] ] \
    [ format "%8.4f %8.4f %8.4f" {*}[ lindex $indenter_extents 1 ] ] ]
  set substrate_extents [ measure minmax $substrate ]
  vmdcon -info [ format "     substrate extents: %26s; %26s" \
    [ format "%8.4f %8.4f %8.4f" {*}[ lindex $substrate_extents 0 ] ] \
    [ format "%8.4f %8.4f %8.4f" {*}[ lindex $substrate_extents 1 ] ] ]

  set xy_shift [ list \
    [ expr [lindex $cell_center 0] - [lindex $in_com 0] ] \
    [ expr [lindex $cell_center 1] - [lindex $in_com 1] ] \
    0. ]

  # maximum substrate z coordinate minus minimum indenter z coordinate
  set initial_substrate_indenter_distance \
    [ expr [lindex $indenter_extents {0 2}] - [lindex $substrate_extents {1 2}] ]
  vmdcon -info "initial substrate-indenter distance: [format "%8.4f" $initial_substrate_indenter_distance]"
  vmdcon -info "desired substrate-indenter distance: [format "%8.4f" $desired_distance]"
  set z_shift [ list 0. 0. \
    [expr $desired_distance - $initial_substrate_indenter_distance] ]

  vmdcon -info "intended xy-shift:    [format "%8.4f %8.4f %8.4f" {*}$xy_shift]"
  vmdcon -info "intended z-shift:     [format "%8.4f %8.4f %8.4f" {*}$z_shift]"
  set in_offset [ vecadd $xy_shift $z_shift ]
  vmdcon -info "intended total shift: [format "%8.4f %8.4f %8.4f" {*}$in_offset]"

  $indenter moveby $in_offset

  # check
  set sb_in_distance [expr [lindex [measure minmax $indenter] {0 2}] - [lindex [measure minmax $substrate] {1 2}] ]
  vmdcon -info "substrate - indenter apex distance after positioning: $sb_in_distance"
}

proc ::JlhVmd::clip_indenter {} {
  variable system_id
  variable indenter_id
  variable indenter
  variable padding_distance

  set a [ molinfo $system_id get a ]
  set b [ molinfo $system_id get b ]
  set c [ molinfo $system_id get c ]

  vmdcon -info [format "%-30.30s %12d" \
    "# atoms in original indenter:" [$indenter num]]

  # user overlap_distance as padding to box boundary
  set indenter [ atomselect $indenter_id \
    "(x > $padding_distance) and (y > $padding_distance) and (z > $padding_distance) and (x < [expr $a - $padding_distance]) and (y < [expr $b - $padding_distance]) and (z <[expr $c - $padding_distance])"]
  $indenter global

  vmdcon -info [format "%-30.30s %12d" \
      "# atoms in clipped indenter:" [$indenter num]]
  # check
}


proc ::JlhVmd::merge {} {
  variable system_id
  variable indenter_id
  variable combined_id

  variable surfactant_resname
  variable counterion_resname
  variable solvent_resname
  variable substrate_resname

  variable system
  variable substrate

  variable indenter
  variable surfactant
  variable counterion
  variable solvent
  variable nonsolvent

  variable overlap_distance

  set combined_id [::TopoTools::selections2mol "$system $indenter"]
  # transfer box measures
  molinfo $combined_id set {a b c alpha beta gamma} \
    [molinfo $system_id get {a b c alpha beta gamma}]

  # assumes substrate annd identer of same material
  set indenter [atomselect $combined_id "index >= [$system num]"]
  $indenter global
  vmdcon -info [format "%30s" "#atoms in indenter $substrate_resname:"] \
    [format "%12d" [$indenter num]]

  set substrate [atomselect $combined_id \
    "resname $substrate_resname and index < [$system num]"]
  $substrate global
  vmdcon -info [format "%30s" "#atoms in substrate $substrate_resname:"] \
    [format "%12d" [$substrate num]]

  set counterion [atomselect $combined_id "resname $counterion_resname"]
  $counterion global
  vmdcon -info [format "%30s" "#atoms in $counterion_resname:"] \
    [format "%12d" [$counterion num]]

  set solvent [atomselect $combined_id "resname $solvent_resname"]
  $solvent global
  vmdcon -info [format "%30s" "#atoms in $solvent_resname:"] \
    [format "%12d" [$solvent num]]

  set surfactant [atomselect $combined_id "resname $surfactant_resname"]
  $surfactant global
  vmdcon -info [format "%30s" "#atoms in $surfactant_resname:"] \
    [format "%12d" [$surfactant num]]

  set nonsolvent [atomselect $combined_id "not resname $solvent_resname"]
  $nonsolvent global
  vmdcon -info [format "%30s" "#atoms not in $solvent_resname:"] \
   [format "%12d" [$nonsolvent num]]

  mol off $system_id
  mol off $indenter_id

  mol rename $combined_id merged
}

proc ::JlhVmd::identify_overlap {} {
  variable combined_id

  variable system
  variable overlap_distance
  variable overlapping
  variable nonoverlapping

  variable surfactant_resname
  variable counterion_resname
  variable solvent_resname
  variable substrate_resname

  atomselect macro immersed "index >= [$system num]"

  set overlapping [atomselect $combined_id \
    "same fragment as (exwithin $overlap_distance of immersed)"]
  $overlapping global
  set nonoverlapping [atomselect $combined_id \
    "not (same fragment as (exwithin $overlap_distance of immersed))"]
  $nonoverlapping global

  # report on overlapping molecules
  variable overlapping_counterion [atomselect \
    $combined_id "resname $counterion_resname and ([$overlapping text])"]
  $overlapping_counterion global
  vmdcon -info [format "%-30.30s" "#atoms in overlapping $counterion_resname:"] \
    [format "%12d" [$overlapping_counterion num]]

  variable overlapping_solvent [atomselect $combined_id \
    "resname $solvent_resname and ([$overlapping text])"]
  $overlapping_solvent global
  vmdcon -info [format "%-30.30s" "#atoms in overlapping $solvent_resname:"] \
    [format "%12d" [$overlapping_solvent num]]

  variable overlapping_surfactant [atomselect $combined_id \
    "resname $surfactant_resname and ([$overlapping text])"]
  $overlapping_surfactant global
  vmdcon -info [format "%-30.30s" "#atoms in overlapping $surfactant_resname:"] \
    [format "%12d" [$overlapping_surfactant num]]

  return $overlapping
}

proc ::JlhVmd::move_nonsolvent_overlap {} {
  variable system_id
  variable combined_id
  variable nonoverlapping
  variable overlapping

  variable substrate
  variable indenter

  variable overlapping_surfactant
  variable overlapping_counterion

  variable surfactant_resname
  variable counterion_resname
  variable solvent_resname
  variable substrate_resname

  variable overlap_distance
  variable bounding_box

  # minmax selection: Returns two vectors, the first containing the minimum
  # x, y, and z coordinates of all atoms in selection, and the second
  # containing the corresponding maxima.

  set indenter_extents  [ measure minmax $indenter ]
  vmdcon -info [ format "      indenter extents: %26s; %26s" \
    [ format "%8.4f %8.4f %8.4f" {*}[ lindex $indenter_extents 0 ] ] \
    [ format "%8.4f %8.4f %8.4f" {*}[ lindex $indenter_extents 1 ] ] ]
  set substrate_extents [ measure minmax $substrate ]
  vmdcon -info [ format "     substrate extents: %26s; %26s" \
    [ format "%8.4f %8.4f %8.4f" {*}[ lindex $substrate_extents 0 ] ] \
    [ format "%8.4f %8.4f %8.4f" {*}[ lindex $substrate_extents 1 ] ] ]

  # before 2019/05/20:
  # determine allowed extents of a new random position for overlapping
  # non-solvent residues. Here box between substrate surface and indenter
  # apex, laterally bounded by indenter extents
  # limits [ [ min_x, min_y, min_z ], [max_x, max_y, max_z] ]
  # set position_limits [ list  \
  #   [ lreplace [ lindex $indenter_extents 0 ] 2 2 [ \
  #     expr [lindex $substrate_extents 1 2] + $overlap_distance ] ] \
  #   [ lreplace [ lindex $indenter_extents 1 ] 2 2 [ \
  #     expr [lindex $indenter_extents 0 2 ] - $overlap_distance ] ] ]

  # since 2019/05/20:
  # determine allowed extents of a new random position for overlapping
  # non-solvent residues. Here system's the substrate's upper coordinate a
  # and otherwise the lateral bounding box coordinates are used:
  # limits [ [ min_x, min_y, min_z ], [max_x, max_y, max_z] ]
  set position_limits [ list \
    [ list \
      [ expr [ lindex $bounding_box 0 0 ] + $overlap_distance ] \
      [ expr [ lindex $bounding_box 0 1 ] + $overlap_distance ] \
      [ expr [ lindex $substrate_extents 1 2 ] + $overlap_distance ] \
    ] [ list \
      [ expr [ lindex $bounding_box 1 0 ] - $overlap_distance ] \
      [ expr [ lindex $bounding_box 1 1 ] - $overlap_distance ] \
      [ expr [ lindex $bounding_box 1 2 ] - $overlap_distance ] \
    ] \
  ]

  vmdcon -info [ format "random position limits: %26s; %26s" \
    [ format "%8.4f %8.4f %8.4f" {*}[ lindex $position_limits 0 ] ] \
    [ format "%8.4f %8.4f %8.4f" {*}[ lindex $position_limits 1 ] ] ]

  set position_origin [ lindex $position_limits 0 ]
  vmdcon -info [ format "random position origin: %26s" \
    [ format "%8.4f %8.4f %8.4f" {*}$position_origin ] ]
  set position_offset [ vecscale -1.0 [ vecsub {*}$position_limits ] ]
  vmdcon -info [ format "random position offset: %26s" \
    [ format "%8.4f %8.4f %8.4f" {*}$position_offset ] ]
  # vmd's "residue" and "resid" differ:
  # former is a zero-based consecutive index,
  # latter is identifier as in input file (starts at 1 for standard numbering)
  # but might not be unique

  # get unique residue ids in overlapping surfactant and counterions
  set overlapping_surfactant_residues [ lsort -unique -integer [ \
    $overlapping_surfactant get residue ] ]
  vmdcon -info [ format "% 8d surfactant residues to be moved." \
    [ llength $overlapping_surfactant_residues ] ]

  set overlapping_counterion_residues [ lsort -unique -integer [ \
    $overlapping_counterion get residue ] ]
  vmdcon -info [ format "% 8d counterion residues to be moved." \
    [ llength $overlapping_counterion_residues ] ]

  set overlapping_residues [ \
    concat $overlapping_surfactant_residues $overlapping_counterion_residues ]

  set nmoved 0
  foreach res $overlapping_residues {
    vmdcon -info [ format "Treating overlapping residue %8d..." $res ]
    set cur [ atomselect $combined_id "residue $res"]
    if { [ $cur num ] == 0 } {
      vmdcon -warn "Selection empty, already removed!"
      continue
    }

    vmdcon -info [format "%-30.30s" "  #atoms in current residue:"] \
      [format "%12d" [$cur num]]

    # only try 1000 times to avoid endless loop
    for {set i 0} {$i<1000} {incr i} {
      set cur_overlapping [atomselect $combined_id \
        "same fragment as (exwithin $overlap_distance of residue $res)"]
      vmdcon -info [format "%-30.30s" "    #atoms overlapping in total:"] \
        [format "%12d" [$cur_overlapping num]]

      # report on overlapping molecules
      set cur_overlapping_counterion [atomselect \
        $combined_id "resname $counterion_resname and ([$cur_overlapping text])"]
      vmdcon -info [format "%-30.30s" "    #atoms in overlapping $counterion_resname:"] \
        [format "%12d" [$cur_overlapping_counterion num]]

      set cur_overlapping_surfactant [atomselect $combined_id \
        "resname $surfactant_resname and ([$cur_overlapping text])"]
      vmdcon -info [format "%-30.30s" "    #atoms in overlapping $surfactant_resname:"] \
        [format "%12d" [$cur_overlapping_surfactant num]]

      set cur_overlapping_substrate [atomselect $combined_id \
        "resname $substrate_resname and ([$cur_overlapping text])"]
      vmdcon -info [format "%-30.30s" "    #atoms in overlapping $substrate_resname:"] \
        [format "%12d" [$cur_overlapping_substrate num]]

      set cur_overlapping_solvent [atomselect $combined_id \
        "resname $solvent_resname and ([$cur_overlapping text])"]
      vmdcon -info [format "%-30.30s" "    #atoms in overlapping $solvent_resname:"] \
        [format "%12d" [$cur_overlapping_solvent num]]

      # check whether position is alright:
      # we allow overlap with solvent, which is subsequently removed, but not
      # with anything else, i.e. other surfactant chains, counterions or
      # substrate. This is the case for overlapping solvent atoms = total ovelap
      if { [ $cur_overlapping_solvent num ] == [ $cur_overlapping num ]} {
        vmdcon -info [
          format "  Found suitable location in iteration %3d." $i ]
        # remove solvent by modifying overlapping and nonoverlapping selections
        set previous_noverlap [ $overlapping num ]
        set previous_nnonoverlap [ $nonoverlapping num ]
        set overlapping [ atomselect $combined_id \
          "(not residue $res and [$overlapping text]) or ([$cur_overlapping_solvent text])" ]
        $overlapping global
        vmdcon -info [ \
          format "    Modified overlap selection from %d to %d atoms." \
              $previous_noverlap [ $overlapping num ] ]

        set nonoverlapping [ atomselect $combined_id \
          "(residue $res or [$nonoverlapping text]) and not ([$cur_overlapping_solvent text])" ]
        $nonoverlapping global
        vmdcon -info [ \
          format "    Modified non-overlap selection from %d to %d atoms." \
              $previous_nnonoverlap [ $nonoverlapping num ] ]

        incr nmoved
        break
      }

      vmdcon -info [ format "  Iteration %3d: Make a move..." $i ]

      set random_3vec [ list [ expr rand() ] [ expr rand() ] [ expr rand() ] ]
      vmdcon -info [ format "   random 3-vector: %26s" \
        [ format "%8.4f %8.4f %8.4f" {*}$random_3vec ] ]

      set random_position [ vecadd $position_origin [ \
        vecmul $random_3vec $position_offset] ]
      vmdcon -info [ format "   random position: %26s" \
        [ format "%8.4f %8.4f %8.4f" {*}$random_position ] ]

      set cur_position [measure center $cur]
      vmdcon -info [ format "  current position: %26s" \
        [ format "%8.4f %8.4f %8.4f" {*}$cur_position ] ]

      set cur_offset [vecsub $random_position $cur_position]
      vmdcon -info [ format "           move by: %26s" \
        [ format "%8.4f %8.4f %8.4f" {*}$cur_offset ] ]
      #set cur [atomselect $combined_id "residue $res"]

      $cur moveby $cur_offset
    }
  }
  vmdcon -info [ format "Successfully moved %3d residues." $nmoved ]
}

proc ::JlhVmd::remove_overlap {} {
  variable system_id
  variable combined_id
  variable nonoverlapping
  variable overlapping

  variable surfactant_resname
  variable counterion_resname
  variable solvent_resname
  variable substrate_resname

  variable substrate
  variable indenter
  variable surfactant
  variable counterion
  variable solvent
  variable nonsolvent

  variable desired_distance

  variable indenter_immersed_id [::TopoTools::selections2mol $nonoverlapping]
  variable overlap_id [::TopoTools::selections2mol $overlapping]

  # copy periodic box
  molinfo $indenter_immersed_id set {a b c alpha beta gamma} \
    [molinfo $system_id get {a b c alpha beta gamma}]

  # assumes substrate annd identer of same material, distinguishes by z loc

  # rough selection based on distance, not pretty, not robust
  # does not account for substrate thickness
  set indenter [atomselect $indenter_immersed_id  \
    "resname $substrate_resname and z >= $desired_distance"]
  $indenter global
  vmdcon -info [format "%30s" "#atoms in indenter $substrate_resname:"] \
    [format "%12d" [$indenter num]]

  set substrate [atomselect $indenter_immersed_id \
    "resname $substrate_resname and z < $desired_distance"]
  $substrate global
  vmdcon -info [format "%30s" "#atoms in substrate $substrate_resname:"] \
    [format "%12d" [$substrate num]]

  set counterion [atomselect $indenter_immersed_id \
    "resname $counterion_resname"]
  $counterion global
  vmdcon -info [format "%30s" "#atoms in $counterion_resname:"] \
    [format "%12d" [$counterion num]]

  set solvent [atomselect $indenter_immersed_id "resname $solvent_resname"]
  $solvent global
  vmdcon -info [format "%30s" "#atoms in $solvent_resname:"] \
    [format "%12d" [$solvent num]]

  set surfactant [atomselect $indenter_immersed_id \
    "resname $surfactant_resname"]
  $surfactant global
  vmdcon -info [format "%30s" "#atoms in $surfactant_resname:"] \
    [format "%12d" [$surfactant num]]

  set nonsolvent [atomselect $indenter_immersed_id \
    "not resname $solvent_resname"]
  $nonsolvent global
  vmdcon -info [format "%30s" "#atoms not in $solvent_resname:"] \
    [format "%12d" [$nonsolvent num]]

  mol off $combined_id
  mol off $overlap_id

  mol rename $indenter_immersed_id indenter_immersed
  mol rename $overlap_id overlap

  return $indenter_immersed_id
}

# writes lammps data file, pdb and psf
proc ::JlhVmd::write_out_indenter_immersed { outname } {
  variable indenter_immersed_id
  set sel [atomselect $indenter_immersed_id all]
  topo -molid $indenter_immersed_id writelammpsdata $outname.lammps full
  vmdcon -info "Wrote $outname.lammps"
  vmdcon -warn "The data files created by TopoTools don't contain any \
    potential parameters or pair/bond/angle/dihedral style definitions. \
    Those have to be generated in addition, however, the generated data \
    files contain comments that match the symbolic type names with the \
    corresponding numeric definitions, which helps in writing those input \
     segment. In many cases, this can be easily scripted, too."
  $sel writepsf $outname.psf
  vmdcon -info "Wrote $outname.psf"
  $sel writepdb $outname.pdb
  vmdcon -info "Wrote $outname.pdb"
}

proc ::JlhVmd::write_top_all { outname } {
  set sel [atomselect top all]
  topo writelammpsdata $outname.lammps full
  vmdcon -info "Wrote $outname.lammps"
  vmdcon -warn "The data files created by TopoTools don't contain any \
    potential parameters or pair/bond/angle/dihedral style definitions. \
    Those have to be generated in addition, however, the generated data \
    files contain comments that match the symbolic type names with the \
    corresponding numeric definitions, which helps in writing those input \
    segment. In many cases, this can be easily scripted, too."
  $sel writepsf $outname.psf
  vmdcon -info "Wrote $outname.psf"
  $sel writepdb $outname.pdb
  vmdcon -info "Wrote $outname.pdb"
}


proc ::JlhVmd::show_nonsolvent { {mol_id 0} {rep_id 0} } {
  # atomselect keywords
  # name type backbonetype residuetype index serial atomicnumber element residue
  # resname altloc resid insertion chain segname segid all none fragment pfrag
  # nfrag numbonds backbone sidechain protein nucleic water waters
  # vmd_fast_hydrogen helix alpha_helix helix_3_10 pi_helix sheet betasheet
  # beta_sheet extended_beta bridge_beta turn coil structure pucker user user2
  # user3 user4 x y z vx vy vz ufx ufy ufz phi psi radius mass charge beta
  # occupancy sequence rasmol sqr sqrt abs floor ceil sin cos tan atan asin acos
  # sinh cosh tanh exp log log10 volindex0 volindex1 volindex2 volindex3 volindex4
  # volindex5 volindex6 volindex7 vol0 vol1 vol2 vol3 vol4 vol5 vol6 vol7
  # interpvol0 interpvol1 interpvol2 interpvol3 interpvol4 interpvol5 interpvol6
  # interpvol7 at acidic cyclic acyclic aliphatic alpha amino aromatic basic
  # bonded buried cg charged hetero hydrophobic small medium large neutral polar
  # purine pyrimidine surface lipid lipids ion ions sugar solvent glycan carbon
  # hydrogen nitrogen oxygen sulfur noh heme conformationall conformationA
  # conformationB conformationC conformationD conformationE conformationF drude
  # unparametrized addedmolefacture qwikmd_protein qwikmd_nucleic qwikmd_glycan
  # qwikmd_lipid qwikmd_hetero

  variable solvent_resname
  variable substrate
  variable indenter
  variable counterion
  variable bounding_box
  variable bb_center
  variable bb_measure

  # make solid atoms appear as thick beads
  $substrate  set radius 5.0
  if { [ info exists indenter ] > 0 } {
    $indenter   set radius 5.0
  }
  $counterion set radius 3.0

  mol selection not resname $solvent_resname
  mol representation CPK
  # or VDW
  mol color element
  mol material Opaque
  # color by element name

  mol modrep $rep_id $mol_id

  pbc box -on -center origin -shiftcenter $bb_center -molid $mol_id
}

proc ::JlhVmd::show_solvent_only { {mol_id 0} {rep_id 0} } {
  variable solvent_resname
  mol selection resname $solvent_resname
  mol representation lines
  mol color element
  mol material Glass3

  mol modrep $rep_id $mol_id
}

proc ::JlhVmd::show_surfactant_only { {mol_id 0} {rep_id 0} } {
  variable surfactant_resname
  mol selection resname $surfactant_resname

  mol representation CPK
  mol color element
  mol material Opaque

  mol modrep $rep_id $mol_id
}

proc ::JlhVmd::show_overlap { {mol_id 0} {rep_id 0} } {
  variable system
  variable overlap_distance

  mol representation Licorice
  mol color element
  mol material Transparent

  mol selection \
    "same fragment as (exwithin $overlap_distance of (index >= [$system num]))"

  mol modrep $rep_id $mol_id
}

# hides solvent
proc ::JlhVmd::set_visual {} {
  variable substrate
  variable indenter
  variable counterion

  # make solid atoms appear as thick beads
  $substrate  set radius 5.0
  if { [ info exists indenter ] > 0 } {
    $indenter   set radius 5.0
  }
  $counterion set radius 3.0

  display resetview

  color Display Background    gray
  color Display BackgroundTop white
  color Display BackgroundBot gray
  color Element Na            green
  display backgroundgradient on

  # after resetview usually centered top view
  # these should result in a centered side view
  rotate x by -90
  # values set empirically
  translate by 0 0.5 0
  scale by 0.4
}

proc ::JlhVmd::render_scene { outname } {
  render TachyonInternal $outname.tga
}

# initialization without manipulation
proc ::JlhVmd::initialize { system_infile } {
  vmdcon -info "-------------------------------------------------------------"
  vmdcon -info "Read system from LAMMPS data file $system_infile..."
  init_system $system_infile
  vmdcon -info "-------------------------------------------------------------"
  vmdcon -info "Objects in system read from $system_infile:"
  display_system_information
  vmdcon -info "-------------------------------------------------------------"
  vmdcon -info "Make types ascii-sortable to preserve original order..."
  make_types_ascii_sortable
  vmdcon -info "-------------------------------------------------------------"
  vmdcon -info "Objects in system after type renaming:"
  display_system_information
  vmdcon -info "-------------------------------------------------------------"
  vmdcon -info "Populating global selections and variables..."
  populate_selections
  vmdcon -info "-------------------------------------------------------------"
  vmdcon -info "Position bounding box..."
  position_bb
}

# read both system and indenter from lammps data file, merges them and
# removes overlap
proc ::JlhVmd::batch_merge_lmp { system_infile indenter_infile } {
  initialize $system_infile
  vmdcon -info "-------------------------------------------------------------"
  vmdcon -info "Wrap system into bounding box, conserving residues..."
  wrap_residue_into_bb
  vmdcon -info "-------------------------------------------------------------"
  vmdcon -info "Read indenter from LAMMPS data file $indenter_infile..."
  read_indenter_lmp $indenter_infile
  vmdcon -info "-------------------------------------------------------------"
  vmdcon -info "Init indenter..."
  init_indenter
  vmdcon -info "-------------------------------------------------------------"
  vmdcon -info "Position indenter..."
  position_indenter
  vmdcon -info "-------------------------------------------------------------"
  vmdcon -info "Merge systems..."
  merge
  vmdcon -info "-------------------------------------------------------------"
  vmdcon -info "Identify ovelap..."
  identify_overlap
  vmdcon -info "-------------------------------------------------------------"
  vmdcon -info "Reposition non-solvent overlap..."
  move_nonsolvent_overlap
  vmdcon -info "-------------------------------------------------------------"
  vmdcon -info "Remove overlap..."
  return [ remove_overlap ]
}

proc ::JlhVmd::batch_process_lmp { system_infile indenter_infile outname } {
  set out_id [ batch_merge_lmp $system_infile $indenter_infile ]
  vmdcon -info "-------------------------------------------------------------"
  vmdcon -info "Objects in output system:"
  display_system_information $out_id
  vmdcon -info "-------------------------------------------------------------"
  vmdcon -info "Write output..."
  write_out_indenter_immersed $outname
  return $out_id
}

proc ::JlhVmd::render_nonsolvent {} {
  variable out_prefix
  vmdcon -info "-------------------------------------------------------------"
  vmdcon -info "Set visualization properties..."
  set_visual
  vmdcon -info "-------------------------------------------------------------"
  vmdcon -info "Show everything except solvent..."
  show_nonsolvent
  vmdcon -info "-------------------------------------------------------------"
  vmdcon -info "Render snapshot..."
  render_scene $out_prefix
}

proc ::JlhVmd::render_solvent {} {
  variable out_prefix
  vmdcon -info "-------------------------------------------------------------"
  vmdcon -info "Set visualization properties..."
  set_visual
  vmdcon -info "-------------------------------------------------------------"
  vmdcon -info "Show only solvent..."
  show_solvent_only
  vmdcon -info "-------------------------------------------------------------"
  vmdcon -info "Render snapshot..."
  render_scene $out_prefix
}

proc ::JlhVmd::render_surfactant {} {
  variable out_prefix
  vmdcon -info "-------------------------------------------------------------"
  vmdcon -info "Set visualization properties..."
  set_visual
  vmdcon -info "-------------------------------------------------------------"
  vmdcon -info "Show only surfactant..."
  show_surfactant_only
  vmdcon -info "-------------------------------------------------------------"
  vmdcon -info "Render snapshot..."
  render_scene $out_prefix
}

proc ::JlhVmd::batch_process_lmp_visual { system_infile indenter_infile outname } {
  set out_id [ batch_process_lmp $system_infile $indenter_infile $outname ]
  vmdcon -info "-------------------------------------------------------------"
  vmdcon -info "Set visualization properties..."
  set_visual
  vmdcon -info "-------------------------------------------------------------"
  vmdcon -info "Show everything except solvent for output system..."
  show_nonsolvent $out_id
  vmdcon -info "-------------------------------------------------------------"
  vmdcon -info "Render snapshot of output system..."
  render_scene $outname
}

proc ::JlhVmd::read_bb_from_yaml { bb_file } {
    variable bounding_box
    # read bounding box from .yaml file
    set bb [::yaml::yaml2dict -file $bb_file]
    set bounding_box [ list \
      [ list [ dict get $bb xlo ] [ dict get $bb ylo ] [ dict get $bb zlo ] ] \
      [ list [ dict get $bb xhi ] [ dict get $bb yhi ] [ dict get $bb zhi ] ] ]
    vmdcon -info [ format "read bounding box from %s: %26s; %26s" \
        $bb_file \
        [ format "%8.4f %8.4f %8.4f" {*}[ lindex $bounding_box 0 ] ] \
        [ format "%8.4f %8.4f %8.4f" {*}[ lindex $bounding_box 1 ] ] ]
    compute_bb_center
    return $bounding_box
}

proc ::JlhVmd::compute_bb_center { } {
    variable bounding_box
    variable bb_center
    variable bb_measure

    vmdcon -info [ format "         bounding box: %26s; %26s" \
        [ format "%8.4f %8.4f %8.4f" {*}[ lindex $bounding_box 0 ] ] \
        [ format "%8.4f %8.4f %8.4f" {*}[ lindex $bounding_box 1 ] ] ]

    set bb_measures [ vecscale -1.0 [ vecsub {*}$bounding_box ] ]
    set bb_center   [ vecscale 0.5  [ vecadd {*}$bounding_box ] ]
    vmdcon -info [ format "bounding box measures: %26s" \
        [ format "%8.4f %8.4f %8.4f" {*}$bb_measures ] ]
    vmdcon -info [ format "  bounding box center: %26s" \
        [ format "%8.4f %8.4f %8.4f" {*}$bb_center ] ]
}

# namespace eval ::TopoTools::
# adapted from ::TopoTools::retypebonds
proc ::TopoTools::make_bond_types_ascii_sortable {sel} {
  set bondlist  [bondinfo getbondlist $sel type]

  set newbonds {}

  # sort type names as inetgers, select last element (largest) and
  # determine its field width:
  set num_digits [
    string length [ lindex [ lsort -integer [ topo bondtypenames ] ] end ] ]

  vmdcon -info "Prepending zeros to bond types filling ${num_digits} digits."

  foreach bond $bondlist {
      set type [format "%0${num_digits}d" [ lindex $bond 2 ]]
      lappend newbonds [list [lindex $bond 0] [lindex $bond 1] $type]
  }
  setbondlist $sel type $newbonds
}

# adapted from proc ::TopoTools::retypeangles
proc ::TopoTools::make_angle_types_ascii_sortable {sel} {
    set anglelist [angleinfo getanglelist $sel]
    set newanglelist {}

    # sort type names as inetgers, select last element (largest) and
    # determine its field width:
    set num_digits [
      string length [ lindex [ lsort -integer [ topo angletypenames ] ] end ] ]

    vmdcon -info "Prepending zeros to angle types filling ${num_digits} digits."
    foreach angle $anglelist {
        lassign $angle type i1 i2 i3
        set type [format "%0${num_digits}d" $type]
        lappend newanglelist [list $type $i1 $i2 $i3]
    }
    setanglelist $sel $newanglelist
}

# adapted from ::TopoTools::retypedihedrals
proc ::TopoTools::make_dihedral_types_ascii_sortable {sel} {
  set dihedrallist [dihedralinfo getdihedrallist $sel]
  set newdihedrallist {}

  # sort type names as inetgers, select last element (largest) and
  # determine its field width:
  set num_digits [
    string length [ lindex [ lsort -integer [ topo dihedraltypenames ] ] end ] ]

  vmdcon -info "Prepending zeros to angle types filling ${num_digits} digits."
  foreach dihedral $dihedrallist {
      lassign $dihedral type i1 i2 i3 i4
      set type [format "%0${num_digits}d" $type]
      lappend newdihedrallist [list $type $i1 $i2 $i3 $i4]
  }
  setdihedrallist $sel $newdihedrallist
}

# adapted from ::TopoTools::retypeimpropers
proc ::TopoTools::make_improper_types_ascii_sortable {sel} {
  set improperlist [improperinfo getimproperlist $sel]
  set newimproperlist {}
  set num_digits [
    string length [ lindex [ lsort -integer [ topo impropertypenames ] ] end ] ]

  vmdcon -info "Prepending zeros to improper types filling ${num_digits} digits."
  foreach improper $improperlist {
      lassign $improper type i1 i2 i3 i4
      set type [format "%0${num_digits}d" $type]
      lappend newimproperlist [list $type $i1 $i2 $i3 $i4]
  }
  setimproperlist $sel $newimproperlist
}

interp alias {} jlh {} ::JlhVmd::jlh
package provide jlhvmd $::JlhVmd::version
