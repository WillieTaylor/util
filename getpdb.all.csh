cp /mb/data/pdb/data/structures/all/pdb/pdb$argv[1].ent.Z $argv[1].Z
if ( -e $argv[1] ) rm $argv[1]
gunzip $argv[1].Z
if ($#argv > 1) then
	if ( $argv[2] == '#' ) then
		echo getting single-chain
		util/getcas $argv[1] > pdb/$argv[1]
	else
		echo getting chain $argv[2]
		util/getcas $argv[1]$argv[2] > pdb/$argv[1]$argv[2]
	endif
else
	echo getting single chain
	util/getcas $argv[1] > pdb/$argv[1]
endif
rm $argv[1]
