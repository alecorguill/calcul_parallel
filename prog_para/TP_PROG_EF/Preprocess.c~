/* preprocess.c */
/* A program in C to partition a single mesh into a number of distributed
sub-meshes suitable for performing a parallel finite element solution*/

#include <stdio.h>
#include <stdlib.h>
/* A number of global data structures and variables relating to the mesh and its
partition are declared here */

typedef struct
{
	float Pos[2];
	int GlobalDOF;
	int LocalNode;
	int Type;
} NodeType;
typedef struct
{
	int Vertex[3];
} ElementType;


/* The main program reads in a mesh from a file followed by a
list of which elements should be assigned to whixh subdomain */

main(int argc, char** argv)
{
  FILE *InFile, *OutFile;
  int Nodes, Elements, I, J, K, V, *NodeCounter, *Common;
  int Increment, S1, S2, S3;
  char FileName[40];
  char *inputFileName;
  NodeType *Node;
  ElementType *Element;
  
	if( argc != 3 )
	{
		fprintf(stderr, "usage: %s <number of partitions> <input mesh file>", argv[0]);
		return -1;
	}

  int PROCNO  = atoi(argv[1]);
  NodeCounter =  (int*) malloc(PROCNO*sizeof(int));
  Common      =  (int*) malloc(PROCNO*sizeof(int));

  inputFileName = argv[2];
   
  InFile = fopen(inputFileName,"r");
  
  if( !InFile )
  {
	  fprintf(stderr,"Can't open %s for reading\n", inputFileName );
	  return -2;
  }

  printf("Generating the file for %d partitions\n", PROCNO );
  
  fscanf(InFile,"%d",&Nodes);
  printf("Reading %d nodes\n",Nodes);

  Node = (NodeType*) malloc(Nodes*sizeof(NodeType));
  for (I=0; I<Nodes; I++) {
    fscanf(InFile, "%d%f%f", &Node[I].GlobalDOF, &Node[I].Pos[0],&Node[I].Pos[1]);
    Node[I].Type = 0;
  }

  fscanf(InFile,"%d",&Elements);
  printf("Reading %d elements\n",Elements);

  Element = (ElementType*) malloc(Elements*sizeof(ElementType));
  for (I=0; I<Elements; I++) {
    fscanf(InFile, "%d%d%d",&S1,&S2,&S3); /*Modifications here due to the C prog which counts the nodes from 0 and not 1 */
    Element[I].Vertex[0] = S1-1; Element[I].Vertex[1] = S2 -1; Element[I].Vertex[2] = S3-1;
  }
  int PartSize[PROCNO];
  int *PartE[PROCNO];
  for (I=0; I<PROCNO; I++) {
    PartE[I] = (int*) malloc(Elements*sizeof(int));
    PartSize[I] = 0;
  }
  for (J=0; J<Elements; J++) {
    fscanf(InFile,"%d", &I);
    PartE[I][PartSize[I]++] = J;
  }
  fclose(InFile);
  int *PartN[PROCNO];
  for (I=0; I<PROCNO; I++)
  {
    PartN[I] = (int*) malloc(Nodes*sizeof(int));
  }

  for (I=0; I<PROCNO; I++) {
    NodeCounter[I] = 0;
    for (J=0; J<Nodes; J++) {
      Node[J].LocalNode = -1;
      PartN[I][J] = 0;
    }
    for (J=0; J<PartSize[I]; J++) {
      for (K=0; K<3; K++) {
        V = Element[PartE[I][J]].Vertex[K];
        if (Node[V].LocalNode == -1)
          Node[V].LocalNode = NodeCounter[I]++;
        if (Node[V].GlobalDOF > -1 && PartN[I][V] == 0) {
          Node[V].Type++;
          PartN[I][V] = 1;
        }
      }

    }
    printf("I = %d Node counter = %d\n",I,NodeCounter[I]);
  }

  /* For each processor a different file is created containing the sub-mesh
     that is to be stored on that processor. In addition to this sub-mesh each
     file contains lists of which nodes in that file are shared with each of the
     other processors. */

  for (I=0; I<PROCNO; I++) {
    if (I<10)
      sprintf(FileName,"Data0%d.In", I);
    else
      sprintf(FileName,"Data%d.In", I);
    OutFile = fopen(FileName, "w");
    fprintf(OutFile,"%d %d\n", NodeCounter[I], PartSize[I]);
    NodeCounter[I] = 0;
    for (J=0; J<Nodes; J++)
      Node[J].LocalNode = -1;

    for (J=0; J<PartSize[I]; J++) {
      for (K=0; K<3; K++) {
        V = Element[PartE[I][J]].Vertex[K];
        /*   printf("I %d J %d K %d PartE[I][J] %d \n",I,J,K,PartE[I][J]);}*/
      if (Node[V].LocalNode == -1) {
        Node[V].LocalNode = NodeCounter[I]++;
        fprintf(OutFile,"%f %f %d\n", Node[V].Pos[0], Node[V].Pos[1], Node[V].Type);
      }
    }
  }
  for (J=0; J<PartSize[I]; J++) {
    for (K=0; K<3; K++) {
      V = Element[PartE[I][J]].Vertex[K];
      fprintf(OutFile, " %d", Node[V].LocalNode);
    }
    fprintf(OutFile,"\n");
  }
  for (K=0; K<PROCNO; K++) {
    Common[K] = 0;
    if( I != K )
      for (J=0; J<Nodes; J++)
        if (PartN[I][J] > 0 && PartN[K][J] > 0)
          Common[K]++;
  }
  for (K=0; K<PROCNO; K++) {
    fprintf(OutFile," %d\n", Common[K]);
    if (Common[K] > 0) {
      for( J=0; J<Nodes; J++)
        if (PartN[I][J] > 0 && PartN[K][J] > 0)
          fprintf(OutFile," %d",Node[J].LocalNode);
      fprintf(OutFile,"\n");
    }
  }
  fclose(OutFile);
}

}
