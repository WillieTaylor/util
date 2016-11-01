cat $argv[1] | grep -v '^$' | awk '{s+=$1; ss+=$1*$1; n++; m=s/n; if(n>1) z=sqrt((ss-2*m*s+n*m*m)/(n-1)); print n,m,z}' | tail -1
