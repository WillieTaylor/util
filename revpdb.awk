grep ATOM pdb/1ntmC | sort -nr +1 | awk '{n++; printf "ATOM   %4d  CA  GLY A%4d     %7.3f %7.3f %7.3f  1.00  0.00\n",n,n,$7,$8,$9}'
