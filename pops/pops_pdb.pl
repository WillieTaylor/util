#! /usr/bin/perl
#==============================================================================
# $Id$
# pops_pdb.pl : calculate POPS area of protein structure
#==============================================================================

#------------------------------------------------------------------------------
# command line parameters
$pdbfilename = $ARGV[0];
$outfilename = $ARGV[1];

if (!$ARGV[0] || !$ARGV[1])
{
    print ("Syntax: area_scop.pl <pdbfilename> <outfilename>\n");
    exit;
}

&prepare_pdb();
&calc_pops();

#------------------------------------------------------------------------------
# subroutines
#------------------------------------------------------------------------------
sub prepare_pdb
{
    open(PDB, "$pdbfilename") || die "Failed reading '$pdbfilename'\n";
    open(PDBPOPS, ">pops.pdb")  || die "Failed writing 'pops.pdb'\n";
    $modc = 0;
    $natom = 0;
    $rna = 0;
    $nres = 0;

    # extract/substitute PDB lines
    while (<PDB>)
    {
        if ($_ =~ /^ENDMDL/) { last; }
        if ($_ =~ /^END/) { last; }

        if ($_ =~ /  CA  /) { ++ $nres; }

        # grep out and substitute (if needed) the backbone atoms
        if ($_ =~ /^ATOM.{9}CA /)
        {
            s/ O1  .{3}/ OXT OXT/; # clean terminal oxygens
            s/ O2  .{3}/ OXT OXT/; # clean terminal oxygens
            s/ OXT .{3}/ OXT OXT/; # clean terminal oxygens
            s/ OT  .{3}/ OXT OXT/;
            s/ CD  ILE/ CD1 ILE/; # clean ILE CD1 atom

            printf PDBPOPS "$_";
            ++ $natom;
        }

        # detected RNA, switch to POPS-R-RNA.dat
        if ($_ =~ /^ATOM.{9}O2 /) { ++ $rna; }
    }

    # print the end mark
    printf PDBPOPS "END\n";

    close(PDB);
    close(PDBPOPS);

    if ($rna == 0)
    { system("cp POPS-R-DNA.dat sas.data"); }
    else
    { system("cp POPS-R-RNA.dat sas.data"); }
}

#------------------------------------------------------------------------------
# calculate POPS area
sub calc_pops
{
    system ("./pops_r pops.pdb > $outfilename");
    system ("rm sas.data");
}

