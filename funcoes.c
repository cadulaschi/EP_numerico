#include "funcoes.h"


void verMatriz (double **ver, int n) {
    printf ("\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf ("%lf \t", ver[i][j]);
        }
        printf ("\n");
    }
    printf("\n");
    //free(ver);
}

void verVetordouble (double *ver, int n) {
    //printf ("\n");
    for (int i = 0; i < n; i++) {
        printf ("%lf  ", ver[i]);
    }
    printf ("\n");
    //free(ver);
}

void verVetorint (int *ver, int n) {
    //printf ("\n");
    for (int i = 0; i < n; i++) {
        printf ("%d  ", ver[i]);
    }
    printf ("\n");
    //free(ver);
}

void decomposicaoLU(double **temp, int *p, int n) {
double somai = 0;
double somaj = 0;
double max, linhatemp;
int l;

    for (int k = 0; k < n; k++) {
        for (int i = k; i < n; i++) {
            for (int j = 0; j < k; j++) somaj = somaj + temp[i][j]*temp[j][k];
            temp [i][k] = temp [i][k] - somaj;
            somaj = 0;
        }

        max = fabs(temp[k][k]);
        l = k;
        for (int i = k; i < n; i++) { 
            if ((fabs(temp [i][k]) > max) && i > k) {
                max = fabs(temp[i][k]);
                l = i; 
            }
        }
        
        p[k] = l;

        //troca entre as linha k e p(k) = l
        if (p[k] != k) { 
            for (int j = 0; j < n; j++) {
                linhatemp = temp[k][j];
                temp[k][j] = temp[p[k]][j];
                temp[p[k]][j] = linhatemp;
            }
        }
        //Parte final do algoritmo        
        for (int j = k + 1; j < n; j++) {
            for (int i = 0; i < k; i++) somai = somai + temp[k][i]*temp[i][j];
            temp[k][j] = temp[k][j] - somai;
            temp[j][k] = temp[j][k]/temp[k][k];
            somai = 0;
        }
    }
}

void solucao (double **temp, double *b, double *x, int *p, int n) {
double btemp;
double soma = 0;
    //Mudando a ordem do vetor b
    for (int i = 0; i < n; i++) {
        if (p[i] != i) {
            btemp = b[i];
            b[i] = b[p[i]];
            b[p[i]] = btemp;
        }
    }

    //Alterando o vetor b a partir dos multiplicadores da matriz A (temp)
    for (int i = 0; i < n; i++) {
        for (int k = i + 1; k < n; k++) {
            b[k] = b[k] - temp[k][i]*b[i];
        }
    }

    //printf("Vetor b: ");
    //for (int i = 0; i < n; i++) printf("%f ", b[i]);
    //printf ("\n\n");

    //Resolvendo o sistema A*x = b
    x[n - 1] = b[n - 1]/temp[n - 1][n - 1];
    for (int i = n - 2; i >= 0; i--) {
        for (int k = i + 1; k < n; k++) {
            soma = soma + temp[i][k]*x[k];
        }
        x[i] = (b[i] - soma)/temp[i][i];
        soma = 0;
    }
}

void constroiJ (double** A, double* x, int n) {
    
}

void constroib (double* b, double *x, int n) {
    b[0] = -(2*x[0] - x[1] - (pow(e, x[0]))/((n+1)*(n+1))); 
    for (int i = 1; i < n - 1; i++) {
      b[i] = -(-x[i-1] + 2*x[i] - x[i+1] -((pow(e, x[i]))/((n+1)*(n+1))));
    }
    b[n-1] = -(-x[n-2] + 2*x[n-1] - (pow(e, x[n-1]))/((n+1)*(n+1)));   
}

void metodoDeNewton (double** A, double* x, double* b, int n) {
//double** A = (double **)calloc(n, sizeof(double*));
double** temp = (double **)calloc(n, sizeof(double*));
for (int i = 0; i < n; i++) {
    //A[i] = (double *)calloc(n, sizeof(double*));
    temp[i] = (double *)calloc(n, sizeof(double*));
}

double* c = (double *)calloc(n, sizeof(double));
int* p = (int *)calloc(n, sizeof(int));

    printf ("\n\n\n\nValores de iniciais utilizados: ");
    for (int i = 0; i < n; i++) printf ("%.1lf | ", x[i]);
        printf ("\n\n");
    
    int iter = 0;
    while (iter < 10) {
        iter++;
        constroiJ (A, x, n);
        //ver(A, n);
        constroib (b, x, n);
        /*for (int i = 0; i < n; i++) {
            printf ("%f ", b[i]);
        }
        printf("\n");
        */
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                temp [i][j] = A[i][j];
            }
        }

        decomposicaoLU(temp, p, n);
        solucao (temp, b, c, p, n);

        for (int i = 0; i < n; i++) {
            x[i] = x[i] + c[i];
        }

        printf ("Solucao do sistema %d:\n", iter);
        for (int i = 0; i < n; i++) printf("          %.8lf \n", x[i]);
        printf ("\n");
    }
    printf ("Numero de iteracoes: %d\n", iter);
    printf ("\n\n\n\n\n\n");
    /*
    printf ("Matriz A\n");
    ver(A, n);

    printf ("Vetor b inicial: ");
    for (int i = 0; i < n; i++) printf("%f ", b[i]);
    printf ("\n\n");

    /********DECOMPOR A MATRIZ A EM LU********/
    //decomposicaoLU(temp, p, n);

    /********SOLUCAO DO SISTEMA*********/
    //solucao (temp, b, c, p, n);*/    

     free(temp); free(c); free(p); free(A);
}

//CALCULANDO VETOR PCALC DE TAMNAHO nBarras
void calculaPcalc(int nBarras, double* tensao, double* vetorTeta, double** G, double** B, double* Pcalc) { //nPQ + nPV
    double soma = 0;
    for(int j = 0; j < nBarras; j++) {// consideramos todas as barras
        for(int k = 0; k < nBarras; k++) {// consideramos as barras Swing (nS)
            if(j != k) {
                soma = soma + tensao[k]*(G[j][k]*cos(vetorTeta[k]-vetorTeta[j]) - B[j][k]*sin(vetorTeta[k]-vetorTeta[j]));
            }
            else
                Pcalc[j] = (tensao[j]*tensao[j]*G[j][j]);
        }
        Pcalc[j] = Pcalc[j] + soma*tensao[j];
        soma = 0;
    }
}
//CALCULANDO VETOR QCALC DE TAMNAHO nBarras
void calculaQcalc(int nBarras, double* tensao, double* vetorTeta, double** G, double** B, double* Qcalc) { //nPQ + nPV
    double soma = 0;
    for(int j = 0; j < nBarras; j++) {// consideramos todas as barras
        for(int k = 0; k < nBarras; k++) {// consideramos as barras Swing (nS)
            if(j != k) {
                soma = soma + tensao[k]*(G[j][k]*sin(vetorTeta[k] - vetorTeta[j]) + B[j][k]*cos(vetorTeta[k]-vetorTeta[j]));
            }
            else
                Qcalc[j] = -(tensao[j]*tensao[j]*B[j][j]);
        }
        Qcalc[j] = Qcalc[j] - soma*tensao[j];
        soma = 0;
    }
}


void constroiDesvioDePotencia (int nPQ, int nPV, double* desvioP, double* PcalcOrdenado, 
                              double* QcalcOrdenado, double* campo4Ordenado) {
    for (int i = 0; i < nPQ + nPV; i++) {
        if (i < nPQ) {
            desvioP[i] = -PcalcOrdenado[i];
            desvioP[nPQ + nPV + i] = -QcalcOrdenado[i];
        }
        else desvioP[i] = -(PcalcOrdenado[i] - campo4Ordenado[i]);
    }
}

void constroiJacobiana (int nBarras, int nPQ, int nPV, double* tensao, double* vetorTeta, 
                        int* tipoDaBarra, double* Pcalc, double* Qcalc, double** G, 
                        double** B, double** J) {

double* fptheta = (double *)calloc(nBarras, sizeof(double*));
double* fpthetaOrdenado = (double *)calloc(nBarras, sizeof(double*));   
double* fpV = (double *)calloc(nBarras, sizeof(double*));
double* fpVOrdenado = (double *)calloc(nBarras, sizeof(double*));
double* fqtheta = (double *)calloc(nBarras, sizeof(double*)); 
double* fqthetaOrdenado = (double *)calloc(nBarras, sizeof(double*));
double* fqV = (double *)calloc(nBarras, sizeof(double*)); 
double* fqVOrdenado = (double *)calloc(nBarras, sizeof(double*));
    for(int j = 0; j < nBarras; j++) {// consideramos todas as barras
        for(int i = 0; i < nBarras; i++) {// consideramos as barras Swing (nS)
            if(j == i) {
                fptheta[i] = -Qcalc[j] - tensao[j]*tensao[j]*B[j][j];
                fpV[i] = Pcalc[j]/tensao[j] - tensao[j]*G[j][j];
                fqtheta[i] = Pcalc[j] - tensao[j]*tensao[j]*G[j][j];
                fqV[i] = Qcalc[j]/tensao[j] - tensao[j]*B[j][j];
                //soma = soma + tensao[i]*(G[j][i]*sin(vetorTeta[i] - vetorTeta[j]) + 
                //B[j][i]*cos(vetorTeta[i] - vetorTeta[j]));
            }
            else
                fptheta[i] = -(tensao[j]*tensao[i]*(G[j][i]*sin(vetorTeta[i] - vetorTeta[j]) + 
                B[j][i]*cos(vetorTeta[i] - vetorTeta[j])));
                
                fpV[i] = tensao[j]*((G[j][i]*cos(vetorTeta[i] - vetorTeta[j]) - 
                B[j][i]*sin(vetorTeta[i] - vetorTeta[j])));

                fqtheta[i] = -tensao[j]*tensao[i]*(G[j][i]);
        }

        ordenaVetor(nPQ, nPV, nBarras, tipoDaBarra, fpthetaOrdenado, fptheta);
        ordenaVetor(nPQ, nPV, nBarras, tipoDaBarra, fpVOrdenado, fptheta);
        ordenaVetor(nPQ, nPV, nBarras, tipoDaBarra, fqVOrdenado, fptheta);
        ordenaVetor(nPQ, nPV, nBarras, tipoDaBarra, fqVOrdenado, fptheta);

        for (int i = 0; i < nPQ + nPV; i++) {
            if (j < nPQ + nPV) J[j][i] = fptheta[i];
        }
    }

free(fptheta); free(fpthetaOrdenado); free(fpV); free(fpVOrdenado); free(fqtheta); 
free(fqthetaOrdenado); free(fqV); free(fqVOrdenado);
}

/*Esse algoritmo ordena um vetor na ordem PQ, PV*/
void ordenaVetor (int nPQ, int nPV, int nBarras, int* tipoDaBarra, double* ordenado, 
                double *original) {
    int jPQ = 0;
    int jPV = nPQ;
    int jS = nPV + nPQ;
    for(int i = 0; i < nBarras; i++) {
        if(tipoDaBarra[i] == 0) ordenado[jPQ] = original[i];
        else jPQ--;
        if (tipoDaBarra[i] == 1) ordenado[jPV] = original[i];
        else jPV--;
        if (tipoDaBarra[i] == 2) ordenado[jS] = original[i];
        else jS--;
        jPQ++; jPV++; jS++;
    }
}


void constroiTroca (int nPQ, int nPV, int nBarras, int* tipoDaBarra, int* troca) {
    int jPQ = 0;
    int jPV = nPQ;
    int jS = nPV + nPQ;
    for(int i = 0; i < nBarras; i++) {
        if(tipoDaBarra[i] == 0) troca[jPQ] = i;
        else jPQ--;
        if (tipoDaBarra[i] == 1) troca[jPV] = i;
        else jPV--;
        if (tipoDaBarra[i] == 2) troca[jS] = i;
        else jS--;
        jPQ++; jPV++; jS++;
    }
}