#include <stdio.h>
#include <stdlib.h> 
#include <math.h>

#define e 2.71828182846
#define tensaoNom1 132790.56
#define PI 3.14159265359

/*IMPRIME UMA MATRIZ*/
void verMatriz (double **ver, int n);
/*IMPRIME UM VETOR DOUBLE*/
void verVetordouble (double *ver, int n);
/*IMPRIME UM VETOR INT*/
void verVetorint (int *ver, int n);

void constroiTroca (int nPQ, int nPV, int nBarras, int* tipoDaBarra, int* troca);
/*FAZ A DECOMPOSICAO LU DE UMA MATRIZ nxn*/
void ordenaVetor (int nPQ, int nPV, int nBarras, int* tipoDaBarra, double* ordenado, 
                  double *original);

void decomposicaoLU(double **temp, int *p, int n);
/*FAZ A SOLUCAO DE UM SISTEMA LINEAR temp*x = b*/
void solucao (double **temp, double *b, double *x, int *p, int n);

/*Calcula o vetor Pcalc de todas as barras j da rede*/
void calculaPcalc(int nBarras, double* tensao, double* vetorTeta, double** G, double** B, double* Pcalc); 
/*Calcula o vetor Qcalc de todas as barras j da rede*/
void calculaQcalc(int nBarras, double* tensao, double* vetorTeta, double** G, double** B, double* Qcalc);

void constroiDesvioDePotencia (int nPQ, int nPV, double* desvioP, double* PcalcOrdenado, 
                              double* QcalcOrdenado, double* campo4Ordenado);
/*CONSTROI A MATRIZ JACOBIANA*/
void constroiJacobiana (int nBarras, int nPQ, int nPV, double* tensao, double* vetorTeta, 
                        int* tipoDaBarra, double* Pcalc, double* Qcalc, double** G, 
                        double** B, double** J, int* troca, double** fptheta, double** fpV,
                        double** fqtheta, double** fqV, double* fpthetaOrdenado, double* fpVOrdenado,
                        double* fqthetaOrdenado, double* fqVOrdenado, double* coluna, double* colunaOrdenado);
/*FAZ O METODO DE NEWTON*/
void metodoDeNewton(int nBarras, int nPQ, int nPV, int nS, double* tensao, double* vetorTeta, int* tipoDaBarra,
                    double* campo4, double** G, double** B, int* troca);

void calculoDa_PotenciaAtiva_e_PerdaAtiva (int nBarras, int nPQ, int nPV, double* tensao, double* vetorTeta, 
                               int* tipoDaBarra, double** G, double** B, double** Sreal, double** Simag,
                               double** perdaAtiva);

void calculoDaPotenciaAtivaDeCarga (int nBarras, double* tensao, double** G, double** B, 
                                    double* consumoDaCarga);