cat $argv[1].pdb | grep -v REMARK | sed 's/HETATM/ATOM  /' | sed 's/MSE/MET/' | grep '^ATOM ' | grep ' A ' | grep ' CA ' | grep -v ' CA B' | sed 's/ CA A/ CA  /' > $argv[1].cas

