set prot = `echo $argv[1] | tr "_" " "`
echo $prot
# cp /mb/skye/home/databases/pdb/ftp.rcsb.org/pub/pdb/data/structures/divided/pdb/*/pdb$prot[1].ent.Z $prot[1].Z
cp /mb/skye/home/databases/pdb/*/pdb$prot[1].ent.Z $prot[1].Z
if ( -e $prot[1].Z ) then
	echo PDB file found
	if ( -e $prot[1] ) rm $prot[1]
	gunzip $prot[1].Z
else
	echo PDB file not found
	exit
endif
if ($#prot > 1) then
	if ( $prot[2] == '#' ) then
		echo getting single-chain
		~/util/getcas $prot[1] | grep -v  ' CA [B-Z]' | sed 's/ CA A/ CA  /' > pdb/$prot[1]
		ls -l pdb/$prot[1]
	else
		echo getting chain $prot[2]
		~/util/getcas $prot[1]$prot[2] | grep -v  ' CA [B-Z]' | sed 's/ CA A/ CA  /' > pdb/$prot[1]$prot[2]
		ls -l pdb/$prot[1]$prot[2]
	endif
else
	echo getting single chain
	~/util/getcas $prot[1] | grep -v  ' CA [B-Z]' | sed 's/ CA A/ CA  /' > pdb/$prot[1]
	ls -l pdb/$prot[1]
endif
rm $prot[1]
