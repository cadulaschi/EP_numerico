#include <stdio.h>
#include <stdlib.h> 
#include <math.h>

#define e 2.71828182846

/*IMPRIME UMA MATRIZ*/
void verMatriz (double **ver, int n);
/*IMPRIME UM VETOR*/
void verVetordouble (double *ver, int n);

void verVetorint (int *ver, int n);
/*FAZ A DECOMPOSICAO LU DE UMA MATRIZ nxn*/
void decomposicaoLU(double **temp, int *p, int n);
/*FAZ A SOLUCAO DE UM SISTEMA LINEAR temp*x = b*/
void solucao (double **temp, double *b, double *x, int *p, int n);
/*CONSTROI A MATRIZ JACOBIANA*/
void constroiJ (double** A, double* x, int n);
/*CONSTROI UM VETOR SOLUCAO b*/
void constroib (double* b, double *x, int n);
/*RESOLVE AS APROXIMACOES PELO METODO DE NEWTON*/
void metodoDeNewton (double** A, double* x, double* b, int n);
/*Calcula o vetor Pcalc de todas as barras j da rede*/
void calculaPcalc(int nBarras, double* tensao, double* vetorTeta, double** G, double** B, double* Pcalc); 
/*Calcula o vetor Qcalc de todas as barras j da rede*/
void calculaQcalc(int nBarras, double* tensao, double* vetorTeta, double** G, double** B, double* Qcalc);

void constroiDesvioDePotencia (int nPQ, int nPV, double* desvioP, double* PcalcOrdenado, 
                              double* QcalcOrdenado, double* campo4Ordenado);

void constroiJacobiana (int nBarras, int nPQ, int nPV, double* tensao, double* vetorTeta, 
                        int* tipoDaBarra, double* Pcalc, double* Qcalc, double** G, 
                        double** B, double** J);

void ordenaVetor (int nPQ, int nPV, int nBarras, int* tipoDaBarra, double* ordenado, 
                  double *original);

void constroiTroca (int nPQ, int nPV, int nBarras, int* tipoDaBarra, int* troca);