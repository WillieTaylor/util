cat $argv[1] | sed 's/>/*@>/' | sed 's/|[a-z]*|/@/' | tr "@" "\n" > $argv[2]
