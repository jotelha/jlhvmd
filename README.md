# JlhVmd

Combines current development versions of

* TopoTools (git@github.com:akohlmey/topotools.git)
* PbcTools  (git@github.com:frobnitzem/pbctools.git)

as submodules as well as some own routines (jlhvmd).
After cloning this repository, clone sub modules as well via

    git submodule init
    git submodule update

*jlhvmd* requires tcllib (for yaml file support).
Tested with release 1.19 obtained at
https://sourceforge.net/projects/tcllib/files/tcllib/1.19/tcllib-1.19.zip
Unpack, create empty installation path of choice, i.e.
`${HOME}/opt/tcllib-1.19/lib/tcllib1.19` and run

    ./configure --prefix ${HOME}/opt/tcllib-1.19/lib/tcllib1.19
    make install

from unpacked root.

Follow
https://gist.github.com/tonigi/a9cfaf7642a7fbc13293#file-install_vmd_plugin-md
to make these plugins available in VMD

One option to export this tree's location as well tcllib's path to
the environment variable `$TCLLIBPATH` before launching *vmd*, i.e.:

    export TCLLIBPATH="${HOME}/opt/tcllib-1.19/lib/tcllib1.19"
    export TCLLIBPATH="${HOME}/git/jlhvmd/lib/plugins/noarch/tcl ${TCLLIBPATH}"
    vmd

### Sample session:

```tcl
package require jlhvmd
package require topotools
jlh set interfaceInfile initial_config.lammps outputPrefix default
jlh use sds
jlh read bb bb.yaml
jlh init
jlh wrap atom
jlh join residue
jlh render nonsolvent
jlh write
```

Apparently, the absolute position of the systems's unit cell is not conserved.
The absolute atom positions, however, and unit cell measures are.

What is more, processing via *topotools * drops the image flags,
which can be a desired side-effect.

## PbcTools

PbcTools requires building `libpbc_core.so`. On Ubuntu 18.04
with `tcl8.5` and `tcl8.5-dev` installed via the apt package management,

These packages provide a configuration script


    source /usr/lib/tcl8.5/tclConfig.sh

that provides information sufficient to compile

    lib/plugins/noarch/tcl/pbctools3.0/src/pbc_core.c

with

    g++ -shared -fPIC -o libpbc_core.so ${TCL_INCLUDE_SPEC} pbc_core.c ${TCL_LIB_SPEC} --verbose

without providing the VMD-specific includes, other than described on
https://www.ks.uiuc.edu/Research/vmd/plugins/doxygen/compiling.html

The resulting `libpbc_core.so` is to be placed at
`lib/plugins/noarch/tcl/pbctools3.0`.

On CentOS 7.6, the above mentioned config file might be found at

  /usr/lib64/tcl8.5/tclConfig.sh

and if Tcl is available via a module system, i.e. LMOD arranged with EasyBuild,
then the according file might be found at

  $EBROOTTCL/lib/tclConfig.sh

after loading the according Tcl module.

## Notes

On Ubuntu, install `rlwrap` via `apt`. Afterwards, replace VMD's default *csh*
launcher (i.e. `/usr/local/bin/vmd`, file starts with shebang `#!/bin/csh`)
with the `bash` alternative, usually available within VMD's binary distribution
as `bin/vmd.sh`. Make sure to adapt the paths in this bash launcher's header.

During development,

    proc reloadPkg pkg {
        eval [package ifneeded $pkg [package require $pkg]]
    }

can be used to dynamically reload modified packages.

When updating the package version number, it has to be done a three places

* at the subtree name, i.e. from

    lib/plugins/noarch/tcl/jlhvmd0.2.4

  to

    lib/plugins/noarch/tcl/jlhvmd0.2.6

* within `pkgIndex.tcl`
* within `jlhvmd.tcl` at the namespace evaluation's head

    namespace eval ::JlhVmd:: {
        variable version 0.2.5
        ...

## LICENSE

The LICENSE applies to content of this repository, not necessarily to the
referenced submodules.
