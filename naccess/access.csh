set EXE_PATH = /home/wtaylor/util/naccess
#naccess_start
#
# The preceding line should be editted to point to
# the directory containing the execuable and default 
# vdw.radii file if things don't work properly.
#
# naccess - run the accessibility program  
#
set nargs = $#argv
if ( $nargs < 1 ) then
   echo "usage naccess pdb_file -p probe_size -r vdw_file -s stdfile -z zslice -[hwyfac]"
   exit
endif    
#
# options:
#      -p = probe size (next arg probe size)
#      -z = accuracy (next arg is accuracy)
#      -r = vdw radii file (next arg is filename)
#      -s = standard accessibilities file (next arg is filename)
#      -h = hetatoms ?
#      -w = waters ?
#      -y = hydrogens ?
#      -f = full asa output format ?
#      -a = produce asa file only, no rsa file ?
#      -c = output atomic contact areas in asa file instead of accessible areas
#
set PDBFILE = 0
set VDWFILE = 0
set STDFILE = 0
set probe = 1.40
set zslice = 0.05
set hets = 0
set wats = 0
set hyds = 0
set full = 0
set asao = 0
set cont = 0 

while ( $#argv )
  switch ($argv[1])
     case -[pP]:	
	shift 	
	set probe = $argv[1]
	breaksw
     case -[zZ]:
        shift
        set zslice = $argv[1]
	breaksw
     case -[hH]:
        set hets = 1
	breaksw
     case -[wW]:
        set wats = 1
        breaksw
     case -[yY]:
        set hyds = 1
        breaksw
     case -[rR]:
        shift
        set VDWFILE = $argv[1]
        breaksw
     case -[sS]:
        shift
        set STDFILE = $argv[1]
        breaksw
     case -[fF]:
        set full = 1
        breaksw
     case -[aA]:
        set asao = 1
        breaksw
     case -[cC]:
        set cont = 1
        breaksw
     default:
        if ( -e $argv[1] && $PDBFILE == 0 ) then
	   set PDBFILE = $argv[1]
        endif
	breaksw
  endsw
  shift
end
#
if ( $PDBFILE == 0 ) then
  echo "usage: you must supply a pdb format file"
  exit(1)
endif
#
# define the VDW radii file
# order of preference
# 1. supplied file with -r option
# 2. "vdw.radii" in current directory
# 3. "vdw.radii" in executable directory
# 4. give up !!
#
if ( $VDWFILE != 0 && ! -e $VDWFILE ) then
   echo "VDW FILE $VDWFILE not found - using default"
   set VDWFILE = 0
endif
if ( $VDWFILE == 0 && -e vdw.radii ) then
   set VDWFILE = vdw.radii
endif
if ( $VDWFILE == 0 && -e $EXE_PATH/vdw.radii ) then
   set VDWFILE = $EXE_PATH/vdw.radii
endif
if ( $VDWFILE == 0 ) then
   echo "FATAL ERROR: unable to assign a vdw radii file"
   exit(1)
endif
#
# is the standard data file present ?
# order of preference
# 1. supplied file with -r option
# 2. "standard.data" in current directory
# 3. "standard.data" in executable directory
# 4. use hard-coded values
#
if ( $STDFILE != 0 && ! -e $STDFILE ) then
   echo "STD FILE $STDFILE not found - using default"
   set STDFILE = 0
endif
if ( $STDFILE == 0 && -e standard.data ) then
  set STDFILE = standard.data
endif
if ( $STDFILE == 0 && -e $EXE_PATH/standard.data ) then
  set STDFILE = $EXE_PATH/standard.data
endif
#
# write input file
#
echo "PDBFILE $PDBFILE" >! accall.input
echo "VDWFILE $VDWFILE" >> accall.input
echo "STDFILE $STDFILE" >> accall.input
echo "PROBE $probe"     >> accall.input
echo "ZSLICE $zslice"   >> accall.input
if ( $hets ) then
  echo "HETATOMS"       >> accall.input
endif
if ( $wats ) then
  echo "WATERS"         >> accall.input
endif
if ( $hyds ) then
  echo "HYDROGENS"      >> accall.input
endif
if ( $full ) then
   echo "FULL"          >> accall.input
endif
if ( $asao ) then
   echo "ASAONLY"       >> accall.input
endif
if ( $cont ) then
   echo "CONTACT"       >> accall.input
endif
cat accall.input
#
# run accessibility calculations
#
$EXE_PATH/accall < accall.input
#
# delete temporary input file
#
#\rm accall.input 
#naccess_end
