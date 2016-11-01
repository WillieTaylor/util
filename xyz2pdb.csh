awk '{n++; printf("ATOM %6d  CA  GLY A%4d     %7.3f %7.3f %7.3f  1.00  0.00\n", n,n,$1,$2,$3)}' $argv[1] > $argv[2]
