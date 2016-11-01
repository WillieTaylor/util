	rotate_all(ac,bc,mat);

rotate_all (ac, bc, mat) double *ac, *bc; Mat mat;
{
FILE	*si1, *si2, *su1, *su2;
float	x, y, z;
int	i, j;
char	line[222], lout[222];
	si1 = fopen("spin2.pdb","r");
	si2 = fopen("spin1.pdb","r");
	su1 = fopen("spun2.pdb","w");
	su2 = fopen("spun1.pdb","w");
	while (1) { int io;
		io = read_line(si1,line);
		if (io < 0 ) break;
		if (io==0) continue;
		strcpy(lout,line);
		if (!strncmp(line,"ATOM   ",7) || !strncmp(line,"HETATM ",7)) {
			sscanf(line+30,"%f %f %f", &x, &y, &z);
			x -= ac[0]; y -= ac[1]; z -= ac[2];
			sprintf(lout+30,"%8.3f%8.3f%8.3f", x,y,z);
			strcat(lout,line+54);
		}
		fprintf(su1,"%s\n", lout);
	}
	while (1) { int io; Vec a;
		io = read_line(si2,line);
		if (io < 0 ) break;
		if (io==0) continue;
		strcpy(lout,line);
		if (!strncmp(line,"ATOM   ",7) || !strncmp(line,"HETATM ",7)) {
			sscanf(line+30,"%f %f %f", &x, &y, &z);
			a.x = x-bc[0]; a.y = y-bc[1]; a.z = z-bc[2];
			MmulV(&mat,a,&a);
			sprintf(lout+30,"%8.3f%8.3f%8.3f", a.x,a.y,a.z);
			strcat(lout,line+54);
		}
		fprintf(su2,"%s\n", lout);
	}
}
