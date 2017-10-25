/* FemPar.c */
/* A program using MPI (with C binding) computes a parallel finite element solution to a linear self-adjoint
elliptic problem using conjugate gradient iterations. */


/* C programme to solve in //:
   - div( a(x,y) grad(u) ) = f(x,y)
with a(x,y) = 1, f(x,y) = 0, Dirichlet BC g(x,y)=x+y

On a given triangulation the pb becomes the following: solve Kv=b

with Kij= int_{Triangle} a*Grad(Nj) . Grad(Ni) dx
here Grad(Nj) (resp Grad(Ni)) is the vector Gr0[j],Gr1[j] (resp Gr0[i],Gr1[i])

The system is stored in block matrix form

for p process

     A1             B1
       A2           B2
K =      .          .                   v = (v1,v2, ....,vp,vs)^t    vi: unknowns on partition i, vs = vertices lying on the partition boundary
          .         .
           .        .
             Ap     Bp
     C1 C2 . . . Cp As       with Ci the tranpose of matrix B1

Right Hand side , b = (b1,b2,...,bp,bs)^t

Each process : reads its own sub-mesh from the file and then assembles its own matrix blocks, Ap, Bp, Fp, As and Fs corresponding to
Ai, Bi, bi, Asi and bsi.


Data structures on each process:
Shared[i]      = number of processes which share the node on the sub-mesh boundary which has the local number i.
Common[i]      = number of nodes shared between this process and process i (this is set to zero when i = the rank of this process).
Neighbours[i,j] = local number of the node on the partition boundary which is the jth node shared with process i.

*/

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

/* A number of global data structures and variables relating to the mesh and its partition are declared
here */

typedef struct
{
  float Pos[2];
  int Local;
  int Type;  /* 
		0   => Condition au bord de dirichet
		1   => noeud interieur
		> 1 => noeud frontière partition
	     */

} NodeType;
typedef struct
{
   int Vertex[3];
} ElementType;
NodeType *Node;
ElementType *Element;
/*
  Shared[i] = nombre de process qui partage le noeud dont le numero local est i
  Common[i] = nb de noeud partagés entre notre process et le process i ( o si i == Me )
  Neighbours[i,j] = numero local du noeud sur la frontière de notre partition qui est le jeme noeud partagé avec le process i
*/
int *Shared, MaxCommon;// Common[PROCNO], *Neighbours[PROCNO];
int *Common, **Neighbours;
int ProcNo, ProcID, Nodes, Elements, IntNodes, IBNodes;
float *Buf;
int PROCNO;
/* MPI_ TAG */
int tag = 99;
/* The following function calculates the inner product of two vectors distributed in a prescribed manner. */

float InnerProduct(float *A1, float *A2, float *B1, float *B2)
{
  int I;
  float res=0.0;
  
  float gamma = 0.0;
  /* Points interieurs */
  for(int i=0; i<IntNodes; ++i){
    gamma += A1[i]*B1[i];
  }
  /* Points frontières */
  for(int i=0; i<IBNodes; ++i){
    gamma += A2[i]*B2[i]/Shared[i];
  }
  MPI_Allreduce(&gamma, &res, 1, MPI_FLOAT,  MPI_SUM, MPI_COMM_WORLD);
  return(res);
}

/* The following function updates entries in the array A which are shared between more than one processor. */

void Update(float *A)
{
   MPI_Status Status;
   int i,j;
   for(i=0; i<(ProcNo) ; ++i){
     for(j=0; j<Common[i]; ++j){
       Buf[j] = A[Node[Neighbours[i][j]].Local];
     }
     if(Common[i] > 0){
       MPI_Send(Buf, Common[i], MPI_FLOAT, i, tag, MPI_COMM_WORLD);
     }
   }
   for(i=0; i<(ProcNo); ++i){
     if(Common[i] > 0)
       MPI_Recv(Buf, Common[i], MPI_FLOAT, i, tag, MPI_COMM_WORLD, &Status);
     for(j=0; j<Common[i]; ++j){
       A[Node[Neighbours[i][j]].Local] += Buf[j];
     }
   }
}

/* The following function implements the conjugate gradient algorithm in a distrixibcuted manner. */

int CG(float **Ap, float **As, float **Bp, float *Fp, float *Fs, float *Vp, float *Vs)
{
   float *Rp, *Rs, *Pp, *Ps, *Qp, *Qs, GammaOld, GammaNew, Tau, Alpha, Beta, Tol=1.0e-4;
   int I, J, K, KMax=250;

   Pp = (float*) malloc(IntNodes*sizeof(float));
   Ps = (float*) malloc(IBNodes*sizeof(float));
   Qp = (float*) malloc(IntNodes*sizeof(float));
   Qs = (float*) malloc(IBNodes*sizeof(float));
   Rp = (float*) malloc(IntNodes*sizeof(float));
   Rs = (float*) malloc(IBNodes*sizeof(float));

   Update(Fs);
   for (I=0; I<IntNodes; I++)
       Pp[I] = Rp[I] = Fp[I];
   for (I=0; I<IBNodes; I++)
       Ps[I] = Rs[I] = Fs[I];

   GammaNew = InnerProduct(Rp,Rs,Rp,Rs);
   if (GammaNew < Tol*Tol)
     return(1);
   if (ProcID == 0)
      printf("Gamma= %f\n",GammaNew);
   for (K=0; K<KMax; K++)
   {
      for (I=0; I<IntNodes; I++) {
         for (J=0, Qp[I]=0.0; J<IntNodes; J++)
             Qp[I] += Ap[I][J]*Pp[J];
         for (J=0; J<IBNodes; J++)
             Qp[I] += Bp[I][J]*Ps[J];
      }
      for (I=0; I< IBNodes; I++) {
          for (J=0, Qs[I]=0.0; J<IntNodes; J++)
              Qs[I] += Bp[J][I]*Pp[J];
          for (J=0; J<IBNodes; J++)
              Qs[I] += As[I][J]*Ps[J];
      }
      Update(Qs);
      Tau = InnerProduct(Pp, Ps, Qp, Qs);
      Alpha = GammaNew/Tau;
      for (I=0; I< IntNodes; I++) {
          Vp[I] += Alpha*Pp[I];
          Rp[I] -= Alpha*Qp[I];
      }
      for (I=0; I< IBNodes; I++) {
          Vs[I] += Alpha*Ps[I];
          Rs[I] -= Alpha*Qs[I];
      }
      GammaOld = GammaNew;
      GammaNew = InnerProduct(Rp,Rs,Rp,Rs);
      if (ProcID == 0)
         printf( "Gamma = %f (K =%d)\n", GammaNew,K);
      if (GammaNew < Tol*Tol)
         return(1);
      Beta = GammaNew/GammaOld;
      for (I=0; I<IntNodes; I++)
          Pp[I] = Rp[I] + Beta*Pp[I];
      for (I=0; I<IBNodes; I++)
          Ps[I] = Rs[I] + Beta*Ps[I];
   }
   return(0);
}

/* The following function applies the Dirichlet boundary conditions. */

float BC( float X, float Y)
{
  return( X+Y );
}

/* The main function reads in the sub-mesh that is to be dealt with and then assembles the
contributions to the stiffness matrix and load vector which come from the elements in this
subdomain. The distributed conjugate gradient solver is then invoked. */

int main( int argc, char** argv)
{
   FILE *InFile, *OutFile;
   char FileName[40], OutName[40];
   int I, J, K, II, JJ;
   float **Ap, **Bp, **As, *Fp, *Fs, *Vp, *Vs;
   float SLoc[3][2], TwoA, Area, Gr0[3], Gr1[3], StiffJI;

   MPI_Init(&argc,&argv);
   double start = MPI_Wtime();
   MPI_Comm_rank(MPI_COMM_WORLD,&ProcID);
   MPI_Comm_size(MPI_COMM_WORLD,&ProcNo);
   PROCNO = ProcNo;
   Common= (int*)malloc(ProcNo*sizeof(int));
   Neighbours= (int**)malloc(ProcNo*sizeof(int *));
   if (ProcID < 10)
      sprintf(FileName, "Data0%d.In", ProcID);
   else
      sprintf(FileName, "Data%d.In", ProcID);

   InFile = fopen(FileName,"r");
   fscanf(InFile, "%d%d", &Nodes, &Elements);

   Node = (NodeType*) malloc(Nodes*sizeof(NodeType));
   for (I=0, IntNodes =0, IBNodes =0; I<Nodes; I++) {
       fscanf(InFile, "%f%f%d", &Node[I].Pos[0],&Node[I].Pos[1],&Node[I].Type);
       if (Node[I].Type == 1)
          Node[I].Local = IntNodes++;
       else if (Node[I].Type > 1)
          Node[I].Local = IBNodes++;
       else
          Node[I].Local = -1;
   }

   Shared = (int*)malloc(IBNodes*sizeof(int));

   for (I=0, IBNodes =0; I<Nodes; I++)
      if (Node[I].Type > 1)
         Shared[IBNodes++]=Node[I].Type;

   Element = (ElementType*) malloc(Elements*sizeof(ElementType));
   for (I=0; I<Elements; I++)
       fscanf(InFile,"%d%d%d",&Element[I].Vertex[0],&Element[I].Vertex[1],&Element[I].Vertex[2]);
   MaxCommon = 0;

   for (I=0; I<PROCNO; I++) {
       fscanf(InFile,"%d", &Common[I]);
       if (Common[I]>MaxCommon)
          MaxCommon = Common[I];
       if (Common[I]>0)
          Neighbours[I] = (int*)malloc(Common[I]*sizeof(int));
       for (J=0; J<Common[I]; J++)
           fscanf(InFile,"%d", &Neighbours[I][J]);
   }

   Buf = (float*) malloc(MaxCommon*sizeof(float));
   fclose(InFile);

   printf("ProcID = %d, IntNodes = %d, IBNodes %d\n", ProcID, IntNodes, IBNodes);

   Ap = (float**) calloc(IntNodes, sizeof(float *));  /* calloc: allocate + initialisation? */
   for (I=0; I<IntNodes; I++)
       Ap[I] = (float*) calloc(IntNodes, sizeof(float));

   Bp = (float**) calloc(IntNodes, sizeof(float *));
   for (I=0; I<IntNodes; I++)
       Bp[I] = (float*) calloc(IBNodes, sizeof(float));

   As = (float**) calloc(IBNodes, sizeof(float *));
   for (I=0; I<IBNodes; I++)
       As[I] = (float*) calloc(IBNodes, sizeof(float));

   Fp = (float*) calloc(IntNodes, sizeof(float));
   Fs = (float*) calloc(IBNodes, sizeof(float));
   Vp = (float*) calloc(IntNodes, sizeof(float));
   Vs = (float*) calloc(IBNodes, sizeof(float));


/*   for (I=0; I<IntNodes; I++) {
       Fp[I]=0.0;
       Vp[I]=0.0;
       for( J=0; J<IntNodes; J++)
          Ap[I][J] = 0.0;
       for( J=0; J<IBNodes; J++)
          Bp[I][J] = 0.0;
   }
   for (I=0; I<IBNodes; I++) {
       Fs[I]=0.0;
       Vs[I]=0.0;
       for( J=0; J<IBNodes; J++)
          As[I][J] = 0.0;
   } */ /* Done within the calloc */


   for (K=0; K<Elements; K++) {
       SLoc[0][0] = Node[Element[K].Vertex[0]].Pos[0];
       SLoc[1][0] = Node[Element[K].Vertex[1]].Pos[0];
       SLoc[2][0] = Node[Element[K].Vertex[2]].Pos[0];
       SLoc[0][1] = Node[Element[K].Vertex[0]].Pos[1];
       SLoc[1][1] = Node[Element[K].Vertex[1]].Pos[1];
       SLoc[2][1] = Node[Element[K].Vertex[2]].Pos[1];

       TwoA = SLoc[0][0]*( SLoc[1][1]-SLoc[2][1] ) +
              SLoc[1][0]*( SLoc[2][1]-SLoc[0][1] ) +
              SLoc[2][0]*( SLoc[0][1]-SLoc[1][1] );

       Area = 0.5*TwoA;
       Gr0[0] = (SLoc[1][1]-SLoc[2][1])/TwoA;
       Gr0[1] = (SLoc[2][1]-SLoc[0][1])/TwoA;
       Gr0[2] = (SLoc[0][1]-SLoc[1][1])/TwoA;
       Gr1[0] = -(SLoc[1][0]-SLoc[2][0])/TwoA;
       Gr1[1] = -(SLoc[2][0]-SLoc[0][0])/TwoA;
       Gr1[2] = -(SLoc[0][0]-SLoc[1][0])/TwoA;

      for (J=0; J<3; J++) {
          JJ = Element[K].Vertex[J];
          if (Node[JJ].Type == 1)
             for (I=0; I<3; I++) {
                 II = Element[K].Vertex[I];
                 StiffJI = Area*(Gr0[I]*Gr0[J]+Gr1[I]*Gr1[J]);
                 if (Node[II].Type == 1)                             /* Int Node */
                    Ap[Node[JJ].Local][Node[II].Local] += StiffJI;
                 else if (Node[II].Type > 1)                         /* Domain Boundary Node */
                    Bp[Node[JJ].Local][Node[II].Local] += StiffJI;
                 else if (Node[II].Type == 0)                        /* Dirichlet Node */
                    {Fp[Node[JJ].Local] -= BC(Node[II].Pos[0], Node[II].Pos[1])*StiffJI;
                    }
             }

          else if (Node[JJ].Type > 1)
             for (I=0; I<3; I++) {
                 II = Element[K].Vertex[I];
                 StiffJI = Area*(Gr0[I]*Gr0[J]+Gr1[I]*Gr1[J]);
                 if (Node[II].Type > 1)
                    As[Node[JJ].Local][Node[II].Local] += StiffJI;
                 else if (Node[II].Type == 0)
                    Fs[Node[JJ].Local] -= BC(Node[II].Pos[0],Node[II].Pos[1])*StiffJI;
             }
      }
   }

   CG(Ap, As, Bp, Fp, Fs, Vp, Vs);

   sprintf(OutName,"Sol0%d.plt",ProcID);
   OutFile = fopen(OutName,"w");
   fprintf(OutFile,"TITLE      =  \"2D FE Unstructured grid data\" \n");
   fprintf(OutFile,"VARIABLES = \"X\", \"Y\", \"U\", \"Dom\" \n");
   fprintf(OutFile,"ZONE T=\"P1\", F=FEPOINT, N= %d E= %d  ET=TRIANGLE \n",Nodes,Elements);
   for (I=0; I<Nodes;I++)
      if (Node[I].Type == 1)
         fprintf(OutFile," %f %f %f %d\n",Node[I].Pos[0],Node[I].Pos[1],Vp[Node[I].Local],ProcID);
       else if (Node[I].Type > 1)
         fprintf(OutFile," %f %f %f %d\n",Node[I].Pos[0],Node[I].Pos[1],Vs[Node[I].Local],ProcID);
       else
         fprintf(OutFile," %f %f %f %d\n",Node[I].Pos[0],Node[I].Pos[1],Node[I].Pos[0]+Node[I].Pos[1],ProcID);

   for (I=0; I<Elements; I++)
       fprintf(OutFile," %d %d %d \n", Element[I].Vertex[0] +1, Element[I].Vertex[1] +1, Element[I].Vertex[2] +1);
   fclose(OutFile);

   double end = MPI_Wtime();
   double max = end-start;
   double current = end-start;
   MPI_Reduce(&current,&max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
   if(ProcID == 0)
     printf("%f\n", max);
   MPI_Finalize();

}

