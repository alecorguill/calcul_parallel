/* small prog to create Mesh.Data from scotch informations */

#include <stdio.h>

main (int argc, char** argv)
{
FILE *InFile1,*InFile2,*OutFile;

int Nodes, Elements, I,J,K,L ;
float x,y;

  if( argc != 4 )
  {
	  fprintf( stderr, "usage: %s <original mesh> <partitionning file> <output file>\n", argv[0]);
	  return -1;
  }

InFile1=fopen(argv[1],"r");
InFile2=fopen(argv[2],"r");
OutFile=fopen(argv[3],"w");

	if( !InFile1 )
	{
		fprintf( stderr, "Can't open file %s for reading\n", argv[1] );
		return -1;
	}
	if( !InFile2 )
	{
		fprintf( stderr, "Can't open file %s for reading\n", argv[2] );
		return -1;
	}
	if( !OutFile )
	{
		fprintf( stderr, "Can't open file %s for writing\n", argv[3] );
		return -1;
	}
	
fscanf(InFile1,"%d\n",&Nodes);
fprintf(OutFile,"%d\n",Nodes);

for (I=0; I<Nodes; I++) {
   fscanf(InFile1,"%d %f %f\n", &K, &x,&y);
   fprintf(OutFile,"%d %f %f\n", K, x,y);
}

fscanf(InFile1,"%d\n",&Elements);
fprintf(OutFile,"%d\n",Elements);

for (I=0; I<Elements; I++) {
   fscanf(InFile1,"%d %d %d\n", &J, &K,&L);
   fprintf(OutFile,"%d %d %d\n", J,K,L);
}

fscanf(InFile2,"%d\n",&Elements);
for (I=0; I<Elements; I++) {
   fscanf(InFile2,"%d %d\n", &J, &K);
   fprintf(OutFile,"%d\n", K );
}




fclose(InFile1);
fclose(InFile2);
fclose(OutFile);


}
