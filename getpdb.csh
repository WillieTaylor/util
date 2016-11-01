set prot = `echo $argv[1] | tr "_" " "`
echo $prot
# cp /mb/skye/home/databases/pdb/ftp.rcsb.org/pub/pdb/data/structures/divided/pdb/*/pdb$prot[1].ent.Z $prot[1].Z
# cp /mb/skye/home/databases/pdb/*/pdb$prot[1].ent.Z $prot[1].Z
cp /mb/databases/pdb/*/pdb$prot[1].ent.gz $prot[1].gz
if ( -e $prot[1].gz ) then
	echo PDB file found
	if ( -e $prot[1] ) rm $prot[1]
	gunzip $prot[1].gz
else
	echo PDB file not found
	exit
endif
cat $prot[1] | sed 's/^HETATM/HETAT /' | awk '{if($4=="MSE"){ print "X"$0}else{print $0}}' > temp1.pdb
cat temp1.pdb | sed 's/XHETAT /ATOM  /' | sed 's/ MSE / MET /' > temp2.pdb
cat temp2.pdb | grep -v  ' CA [B-Z]' | sed 's/ CA A/ CA  /' | grep -v  ' CA [2-9]' | sed 's/ CA 1/ CA  /' > $prot[1]
if ($#prot > 1) then
	set chain = `echo $prot[2] | tr "[a-z]" "[A-Z]"`
	if ( $prot[2] == '#' ) then
		echo getting single-chain
		~/util/getcas $prot[1] > temp.pdb
	else
		echo getting chain $chain
		~/util/getcas $prot[1]$chain > temp.pdb
	endif
	set name = `echo pdb/$prot[1]$chain`
else
	echo getting single chain
	~/util/getcas $prot[1] > temp.pdb
	set name = `echo pdb/$prot[1]`
endif
set len = `cat temp.pdb | grep '^ATOM ' | grep ' CA ' | wc -l`
mv temp.pdb $name
echo Protein = $name, length =  $len, size = `ls -l $name`
#rm $prot[1]
