cp /mb/data/rcsb/data/structures/divided/pdb/*/pdb$argv[1].ent.gz $argv[1].gz
gunzip $argv[1].gz
mv $argv[1] pdb/$argv[1].pdb
