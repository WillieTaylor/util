# 1 = file, 2 = chain ID (otherwise all)

cat $argv[1].pdb | grep -v REMARK | sed 's/HETATM/ATOM  /' | sed 's/MSE/MET/' | grep '^ATOM ' | grep ' CA ' | grep -v ' CA B' | sed 's/ CA A/ CA  /' > temp.cas
if ( $#argv > 1 ) then
	eval "grep ' $argv[2] ' temp.cas" > $argv[1].cas
	rm temp.cas
else
	echo no chain ID
	mv temp.cas $argv[1].cas
endif
