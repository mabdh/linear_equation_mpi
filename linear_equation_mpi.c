/*
 ============================================================================
 Name        : linearparallel.c
 Author      : abduh dong
 Version     :
 Copyright   : Your copyright notice
 Description : Calculate Pi in MPI
 ============================================================================
 */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MATRIX_DIMENSION 4
#define tagMatriks	1
#define tagVector	2
#define taghasilvector	3
#define taghasilmatriks 4

/*****************************************************/
void print_matriks(double* fSolution,double* fTriangular){
	FILE* fout;
	int i,n;
	int nDim = MATRIX_DIMENSION;
	n=1;
	fout=fopen("vektor.txt","wb");
	for(i=0;i<nDim;i++){
			fprintf(fout,"%d ",(int)fSolution[i]);
	}
	fclose(fout);
	fout=fopen("matriks.txt","wb");

		for(i=0;i<nDim*nDim;i++)
		{
			fprintf(fout,"%d ",(int)fTriangular[i]);
			if(i==n*nDim-1){
				fprintf(fout,"\n");
				n++;
			}
		}
	fclose(fout);
}
/***********************************************/

void print_vektor(double* vek){

	int i;
	int nDim = MATRIX_DIMENSION;
		for(i=0;i<nDim;i++)
		{
				printf("%d ",(int)vek[i]);
		}
}

void print_matriks2(double* fTriangular){

	int i;
	int n = 1;
	int nDim = MATRIX_DIMENSION;
		for(i=0;i<nDim*nDim;i++)
		{
				printf("%d ",(int)fTriangular[i]);

			if(i==n*nDim-1){
				printf("\n");
				n++;
			}
		}
}

void transpose(double* input,double* output){
	int n,p,q;
	int i,j,dim;

	dim=MATRIX_DIMENSION;
	
	n=1;p=0;q=0;
	for(j=0;j<dim;j++){
		for(i=j*dim;i<n*dim;i++){
			output[i]=input[p*dim+q];
			p++;
		}
	n++;
	p=0;
	q++;
	}
}

void init_matriks(double* matriks)
{
	int i;
	int nDim = MATRIX_DIMENSION;
	for(i=0;i<nDim*nDim;i++)
	{
		matriks[i]=0;
	}
}

void init_vektor(double* vektor)
{
	int i;
	int nDim = MATRIX_DIMENSION;
	for(i=0;i<nDim;i++)
	{
		vektor[i]=0;
	}
}

int main(int argc, char *argv[])
{
	int	p;		/* number of processes */
	int	my_rank;		/* rank of process */
	MPI_Status	status ;		/* return status for receive */

	double fMaxElem;
	double fAcc;
	int i , j, k, m,x;

	FILE *fin;

	  int nDim = MATRIX_DIMENSION;
	  double fMatr[MATRIX_DIMENSION*MATRIX_DIMENSION];
	  double fVec[MATRIX_DIMENSION];
	  double tempM[MATRIX_DIMENSION*MATRIX_DIMENSION];
	  // double tempV[MATRIX_DIMENSION];
	  double fSolution[MATRIX_DIMENSION];
	double fTriangular[MATRIX_DIMENSION*MATRIX_DIMENSION];
	double fVecTriangular[MATRIX_DIMENSION];

	
		//inisialisasi matriks
		init_matriks(fMatr);
		init_matriks(tempM);
		init_matriks(fTriangular);
		fin=fopen("matriks","r");
		if(fin==NULL){
			printf("\nFile Matriks tidak ada\n");
			exit(-1);
		}
		//baca semua data dari file
		for(i=0;i<(nDim*nDim);i++)
				fscanf(fin,"%lf",&fMatr[i]);
		fclose(fin);
		
		
		//inisialisasi vektor
		init_vektor(fVec);
		init_vektor(fVecTriangular);
		init_vektor(fSolution);

			fin=fopen("vektor","r");
			if(fin==NULL){
				printf("\nFile Matriks tidak ada\n");
				exit(-1);
			}
			//baca semua data dari file
			for(i=0;i<nDim;i++)
				fscanf(fin,"%lf",&fVec[i]);
			fclose(fin);


	/* start up MPI */
	MPI_Init(&argc, &argv);
	/* find out process rank */
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	/* find out number of processes */
	MPI_Comm_size(MPI_COMM_WORLD, &p);

		
			
			  if (my_rank == 0) {

				printf("rank0\n");
				
				// search of line with max element
				for(k=0; k<(nDim-1); k++) // base row of matrix
				{
				fMaxElem = fabs( fMatr[k*nDim + k] );
				m = k;
				for(i=k+1; i<nDim; i++)
				{
					if(fMaxElem < fabs(fMatr[i*nDim + k]) )
					{
						fMaxElem = fMatr[i*nDim + k];
						m = i;
					}
				}
				// permutation of base line (index k) and max element line(index m)
				if(m != k)
				{
					for(i=k; i<nDim; i++)
					{
						fAcc              = fMatr[k*nDim + i];
						fMatr[k*nDim + i] = fMatr[m*nDim + i];
						fMatr[m*nDim + i] = fAcc;
					}
					fAcc = fVec[k];
					fVec[k] = fVec[m];
					fVec[m] = fAcc;
				}
				if( fMatr[k*nDim + k] == 0.) return 1;
				}
				printf("Kirim\n\n");
				printf("Matriks:\n");
				print_matriks2(fMatr);
				printf("\n");
				printf("Vektor:\n");
				print_vektor(fVec);
				printf("\n\n");
				printf("Transpose:\n");
				transpose(fMatr,tempM);
				print_matriks2(fMatr);
				printf("\n");
				//kirim ke prosessor lain
				for (i=1; i<p; i++)
				{		
						MPI_Send(&(tempM[i * (nDim*nDim/p) + ((nDim*nDim)%p)]), (nDim*nDim)/p, MPI_DOUBLE, i, tagMatriks, MPI_COMM_WORLD);
						MPI_Send(&(fVec[i * nDim/p + nDim%p]), nDim/p, MPI_DOUBLE, i, tagVector, MPI_COMM_WORLD);
				}
				for (i = 1; i < p; i++) {
						MPI_Recv(&fMatr[i * (nDim*nDim/p) + (nDim*nDim%p)], (nDim*nDim)/p, MPI_DOUBLE, i, taghasilmatriks, MPI_COMM_WORLD, &status);
						MPI_Recv(&fVec[i * nDim/p + nDim%p], nDim/p, MPI_DOUBLE, i, taghasilvector, MPI_COMM_WORLD, &status);
				}

							/*	x=0;
								for(i=k*nDim;i<k*nDim+nDim;i++)
								{
									fTriangular[i]=fMatr[m+x];
									fVecTriangular[k]=fVec[m];
									x++;
								}
								 print_matriks2(fTriangular);
												  printf("\n");*/
				

			  	} else {
			  		printf("bukan rank nol\n");
	
			  		MPI_Recv(tempM, (nDim*nDim)/p, MPI_DOUBLE, 0, tagMatriks, MPI_COMM_WORLD, &status);
			  		MPI_Recv(fVec, nDim/p, MPI_DOUBLE, 0, tagVector, MPI_COMM_WORLD, &status);
					//print ke console
					printf("Diterima\n\n");
					printf("Matriks:\n");
					print_matriks2(tempM);
					printf("\n");
					printf("Vektor:\n");
					print_vektor(fVec);
					printf("\n\n");
					printf("Transpose:\n");
					transpose(tempM,fMatr);
					print_matriks2(fMatr);
					
					printf("\n");
				 for(k=0; k<(nDim-1); k++) // base row of matrix
				{
			  		//triangular
			  		for(j=(k+1); j<nDim/p; j++) // current row of matrix
			  		{
			  			      fAcc = - fMatr[j*nDim + k] / fMatr[k*nDim + k];
			  			      for(i=k; i<nDim; i++)
			  			      {
			  			    	 fMatr[j*nDim + i] = fMatr[j*nDim + i] + fAcc*fMatr[k*nDim + i];
			  			      }
			  			      	  fVec[j] = fVec[j] + fAcc*fVec[k]; // free member recalculation
			  		}
				}	
			  	
			  		// mengirim hasilnya ke proses 0
			  		for (i=1; i<p; i++){
					
			  			MPI_Send(&(fMatr[i * (nDim*nDim/p) + (nDim*nDim%p)]), (nDim*nDim)/p, MPI_DOUBLE, 0, taghasilmatriks, MPI_COMM_WORLD);
			  			MPI_Send(&(fVec[i * nDim/p + nDim%p]), nDim/p, MPI_DOUBLE, 0, taghasilvector, MPI_COMM_WORLD);
					}
			  	}
			
		  


		  for(k=(nDim-1); k>=0; k--)
		  {
		    fSolution[k] = fVecTriangular[k];
		    for(i=(k+1); i<nDim; i++)
		    {
		      fSolution[k] -= (fTriangular[k*nDim + i]*fSolution[i]);
		    }
		    fSolution[k] = fSolution[k] / fTriangular[k*nDim + k];
		  }
		  printf("selesai\n");
		  print_matriks(fVec,fMatr);
	/* shut down MPI */
	MPI_Finalize(); 

	return 0;
}

