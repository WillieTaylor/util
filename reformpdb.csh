awk '/ATOM/ {printf("ATOM %6d  CA  %s A%4d     %7.3f %7.3f %7.3f %5.2f %5.2f\n", $2,$4,$6,$7,$8,$9,$10,$11)}' $argv[1] > $argv[2]
