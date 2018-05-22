#include "funcoes.h"


void verMatriz (double **ver, int n) {
    printf ("\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf ("%.5lf \t", ver[i][j]);
        }
        printf ("\n");
    }
    printf("\n");
}

void verVetordouble (double *ver, int n) {
    for (int i = 0; i < n; i++) {
        printf ("%lf  ", ver[i]);
    }
    printf ("\n");
}

void verVetorint (int *ver, int n) {
    for (int i = 0; i < n; i++) {
        printf ("%d  ", ver[i]);
    }
    printf ("\n");
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

void metodoDeNewton(int nBarras, int nPQ, int nPV, int nS, double* tensao, double* vetorTeta, int* tipoDaBarra,
                    double* campo4, double** G, double** B, int* troca) {

double* c = (double *)calloc(2*nPQ + nPV, sizeof(double));
int* p = (int *)calloc(2*nPQ + nPV, sizeof(int));

double** fptheta = (double **)calloc(nBarras, sizeof(double*)); 
double** fpV = (double **)calloc(nBarras, sizeof(double*));
double** fqtheta = (double **)calloc(nBarras, sizeof(double*));
double** fqV = (double **)calloc(nBarras, sizeof(double*));   
for (int i = 0; i < nBarras; i++){//nLinhas =nColunas
fptheta[i] = (double *)calloc(nBarras, sizeof(double));
fpV[i] = (double *)calloc(nBarras, sizeof(double));
fqtheta[i] = (double *)calloc(nBarras, sizeof(double));
fqV[i] = (double *)calloc(nBarras, sizeof(double));
}
double* fpthetaOrdenado = (double *)calloc(nBarras, sizeof(double*));  
double* fpVOrdenado = (double *)calloc(nBarras, sizeof(double*));
double* fqthetaOrdenado = (double *)calloc(nBarras, sizeof(double*));
double* fqVOrdenado = (double *)calloc(nBarras, sizeof(double*));

double* Pcalc = (double *)calloc(nBarras, sizeof(double));    
double* PcalcOrdenado = (double *)calloc(nBarras, sizeof(double));
double* Qcalc = (double *)calloc(nBarras, sizeof(double));
double* QcalcOrdenado = (double *)calloc(nBarras, sizeof(double));
double* desvioP = (double *)calloc(2*nPQ + nPV, sizeof(double));
double* campo4Ordenado = (double *)calloc(nBarras, sizeof(double));
double** temp = (double **)calloc(2*nPQ+ nPV, sizeof(double*));
double** J = (double **)calloc(2*nPQ + nPV, sizeof(double *));
for (int i = 0; i < 2*nPQ + nPV; i++){//nLinhas =nColunas
    J[i] = (double *)calloc(2*nPQ + nPV, sizeof(double));
    temp[i] = (double *)calloc(2*nPQ + nPV, sizeof(double));
}
double* coluna = (double *)calloc(nBarras, sizeof(double*));
double* colunaOrdenado = (double *)calloc(nBarras, sizeof(double*));
double maxTheta, maxV;

    int iter = 0;
    do {
        iter++;
        printf("\n ********************** ITERAÇÃO %d ************************", iter);

        calculaPcalc(nBarras, tensao, vetorTeta, G, B, Pcalc);
        
        ordenaVetor(nPQ, nPV, nBarras, tipoDaBarra, PcalcOrdenado , Pcalc);  
        
        calculaQcalc(nBarras, tensao, vetorTeta, G, B, Qcalc);
        
        ordenaVetor(nPQ, nPV, nBarras, tipoDaBarra, QcalcOrdenado , Qcalc);       
        
        constroiJacobiana (nBarras, nPQ, nPV, tensao, vetorTeta, tipoDaBarra, Pcalc, Qcalc,
                           G, B, J, troca, fptheta, fpV, fqtheta, fqV, fpthetaOrdenado, fpVOrdenado,
                           fqthetaOrdenado, fqVOrdenado, coluna, colunaOrdenado);

        ordenaVetor(nPQ, nPV, nBarras, tipoDaBarra, campo4Ordenado, campo4);
        //Forma vetor de desvios já negativado!
        constroiDesvioDePotencia (nPQ, nPV, desvioP, PcalcOrdenado, QcalcOrdenado, campo4Ordenado);  

         for (int i = 0; i < 2*nPQ+ nPV; i++) {
            for (int j = 0; j < 2*nPQ+ nPV; j++) {
                temp [i][j] = J[i][j];
            }
         }
        printf("\nComeçou decomposicao LU\n");
        decomposicaoLU(temp, p, 2*nPQ + nPV);
        printf("\nTerminou decomposicao LU\n");
        printf("\nComeçou solucao\n");
        solucao (temp, desvioP, c, p, 2*nPQ + nPV);
        printf("\nTerminou solucao\n");

        for (int i = 0; i < nPQ + nPV; i++) {
            if (i < nPQ) {
                vetorTeta[troca[i]]= vetorTeta[troca[i]] +  c[i];
                tensao[troca[i]] = tensao[troca[i]] + c[i + nPQ + nPV];
            }
            else if (i >= nPQ) vetorTeta[troca[i]] = vetorTeta[troca[i]] + c[i];
        }
        maxTheta = fabs(c[0]);
        maxV = fabs(c[nPQ + nPV]);

        for (int i = 1; i < nPQ + nPV; i++) {
            if (fabs(c[i]) > fabs(maxTheta)) maxTheta = fabs(c[i]);
        }
        for (int i = nPQ + nPV + 1; i < 2*nPQ + nPV; i++) {
            if (fabs(c[i]) > fabs(maxV)) maxV = fabs(c[i]);
        }
    } while (fabs(maxV) > 1e-6 && fabs(maxTheta) > 1e-5);

    printf ("\n*************************** TÉRMINO *************************");
    printf ("\n************************* %d ITERAÇÕES ***********************", iter);
    printf ("\n\n\n");
    
    for (int i = 0; i < nBarras; i++) {
        free(fptheta[i]); 
        free(fpV[i]);
        free(fqtheta[i]);
        free(fqV[i]);
    }
    
    free(Pcalc); free(PcalcOrdenado); 
    free(Qcalc); free(QcalcOrdenado);
    for (int i = 0; i < 2*nPQ + nPV; i++) {free(J[i]);}
    free(campo4Ordenado);
    free(desvioP);
    free(J);
    free(c); free(p); 
    free(fpthetaOrdenado); free(fpV); free(fpVOrdenado); free(fqtheta); 
    free(fqthetaOrdenado); free(fqV); free(fqVOrdenado);
    free(coluna); free(colunaOrdenado);
}

//CALCULANDO VETOR PCALC DE TAMNAHO nBarras
void calculaPcalc(int nBarras, double* tensao, double* vetorTeta, double** G, double** B, double* Pcalc) { //nPQ + nPV
    double soma = 0;
    for(int j = 0; j < nBarras; j++) {// consideramos todas as barras
        for(int k = 0; k < nBarras; k++) {// consideramos as barras Swing (nS)
            if(j != k) {
                soma = soma + tensao[k]*(G[j][k]*cos(vetorTeta[k]-vetorTeta[j]) 
                       - B[j][k]*sin(vetorTeta[k]-vetorTeta[j]));
            }
            else
                Pcalc[j] = (tensao[j]*tensao[j]*G[j][j]);
        }
        Pcalc[j] = Pcalc[j] + soma*tensao[j];
        soma = 0;
    }

    printf("\nPcalc feito\n");
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
        printf("\nQcalc feito\n");
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
                        double** B, double** J, int* troca, double** fptheta, double** fpV,
                        double** fqtheta, double** fqV, double* fpthetaOrdenado, double* fpVOrdenado,
                        double* fqthetaOrdenado, double* fqVOrdenado, double* coluna, double* colunaOrdenado) {

    printf("\nCOmeçou a Jacobiana\n");
    
    for(int j = 0; j < nBarras; j++) {// consideramos todas as barras
        for(int i = 0; i < nBarras; i++) {// consideramos as barras Swing (nS)
            if(j == i) {
                fptheta[j][i] = -Qcalc[j] - tensao[j]*tensao[j]*B[j][j];
                fpV[j][i] = Pcalc[j]/tensao[j] + tensao[j]*G[j][j];
                fqtheta[j][i] = Pcalc[j] - tensao[j]*tensao[j]*G[j][j];
                fqV[j][i] = Qcalc[j]/tensao[j] - tensao[j]*B[j][j];
            }
            else {
                fptheta[j][i] = -(tensao[j]*tensao[i]*(G[j][i]*sin(vetorTeta[i] - vetorTeta[j]) 
                             + B[j][i]*cos(vetorTeta[i] - vetorTeta[j])));
                
                fpV[j][i] = tensao[j]*((G[j][i]*cos(vetorTeta[i] - vetorTeta[j]) 
                         - B[j][i]*sin(vetorTeta[i] - vetorTeta[j])));

                fqtheta[j][i] = -tensao[j]*tensao[i]*((G[j][i])*cos(vetorTeta[i] - vetorTeta[j])
                             - B[j][i]*sin(vetorTeta[i] - vetorTeta[j]));

                fqV[j][i] = -tensao[j]*(G[j][i]*sin(vetorTeta[i] - vetorTeta[j]) 
                         + B[j][i]*cos(vetorTeta[i] - vetorTeta[j]));
            }
        }

        ordenaVetor(nPQ, nPV, nBarras, tipoDaBarra, fpthetaOrdenado, fptheta[j]);
        for (int i = 0; i < nBarras; i++) fptheta[j][i] = fpthetaOrdenado[i];

        ordenaVetor(nPQ, nPV, nBarras, tipoDaBarra, fpVOrdenado, fpV[j]);
        for (int i = 0; i < nBarras; i++) fpV[j][i] = fpVOrdenado[i];

        ordenaVetor(nPQ, nPV, nBarras, tipoDaBarra, fqthetaOrdenado, fqtheta[j]);
        for (int i = 0; i < nBarras; i++) fqtheta[j][i] = fqthetaOrdenado[i];

        ordenaVetor(nPQ, nPV, nBarras, tipoDaBarra, fqVOrdenado, fqV[j]);
        for (int i = 0; i < nBarras; i++) fqV[j][i] = fqVOrdenado[i];
        
    }

    for (int i = 0; i < nBarras; i++) {
        for (int j = 0; j < nBarras; j++) coluna[j] = fptheta[j][i];
        ordenaVetor(nPQ, nPV, nBarras, tipoDaBarra, colunaOrdenado, coluna);
        for (int j = 0; j < nBarras; j++) fptheta[j][i] = colunaOrdenado[j];

        for (int j = 0; j < nBarras; j++) coluna[j] = fpV[j][i];
        ordenaVetor(nPQ, nPV, nBarras, tipoDaBarra, colunaOrdenado, coluna);
        for (int j = 0; j < nBarras; j++) fpV[j][i] = colunaOrdenado[j];

        for (int j = 0; j < nBarras; j++) coluna[j] = fqtheta[j][i];
        ordenaVetor(nPQ, nPV, nBarras, tipoDaBarra, colunaOrdenado, coluna);
        for (int j = 0; j < nBarras; j++) fqtheta[j][i] = colunaOrdenado[j];

        for (int j = 0; j < nBarras; j++) coluna[j] = fqV[j][i];
        ordenaVetor(nPQ, nPV, nBarras, tipoDaBarra, colunaOrdenado, coluna);
        for (int j = 0; j < nBarras; j++) fqV[j][i] = colunaOrdenado[j];
    }
    for (int i = 0; i < 2*nPQ + nPV; i++) {
        for (int j = 0; j < 2*nPQ + nPV; j++) {
            if ((i < nPQ + nPV) && (j < nPQ + nPV)) J[i][j] = fptheta[i][j];
            else if ((i >= nPQ + nPV) && (j < nPQ + nPV)) J[i][j] = fqtheta[i - (nPQ + nPV)][j];
            else if ((i < nPQ + nPV) && (j >= nPQ + nPV)) J[i][j] = fpV[i][j - (nPQ + nPV)];
            else if ((i >= nPQ + nPV) && (j >= nPQ + nPV)) J[i][j] = fqV[i - (nPQ + nPV)][j - (nPQ + nPV)];
        }
    }

    printf("\nTerminou a Jacobiana\n");    
}

void calculoDa_PotenciaAtiva_e_PerdaAtiva (int nBarras, int nPQ, int nPV, double* tensao, double* vetorTeta, 
                               int* tipoDaBarra, double** G, double** B, double** Sreal,
                               double** Simag, double** perdaAtiva) {
    double* Vreal = (double *)calloc(nBarras, sizeof(double*));
    double* Vimag = (double *)calloc(nBarras, sizeof(double*));
    double** Ireal = (double **)calloc(nBarras, sizeof(double*));
    double** Iimag = (double **)calloc(nBarras, sizeof(double*));   
    for (int i = 0; i < nBarras; i++) {
        Ireal[i] = (double *)calloc(nBarras, sizeof(double));
        Iimag[i] = (double *)calloc(nBarras, sizeof(double));
    }

    for (int i = 0; i < nBarras; i++) {
        Vreal[i] = tensao[i]*cos(vetorTeta[i]);
        Vimag[i] = tensao[i]*sin(vetorTeta[i]);
    }

    for (int i = 0; i < nBarras; i++) {
        for (int j = 0; j < nBarras; j++) {
            Ireal[i][j] = -(Vreal[i] - Vreal[j])*G[i][j] + (Vimag[i] - Vimag[j])*B[i][j]; 
            Iimag[i][j] = -(Vreal[i] - Vreal[j])*B[i][j] - (Vimag[i] - Vimag[j])*G[i][j];

            Sreal[i][j] = (0.001*(Ireal[i][j]*Vreal[i] + Vimag[i]*Iimag[i][j]))*3;
            Simag[i][j] = (0.001*(-Vreal[i]*Iimag[i][j] + Vimag[i]*Ireal[i][j]))*3;
        }
    }
    for (int i = 0; i < nBarras; i++) {
        for (int j = 0; j < nBarras; j++) {
            perdaAtiva[i][j] = ((Vreal[i] - Vreal[j])*(Vreal[i] - Vreal[j]) + 
                               (Vimag[i] - Vimag[j])*(Vimag[i] - Vimag[j]))*(-G[i][j]);
            perdaAtiva[i][j] = (0.001*perdaAtiva[i][j])*3;        
        }
    }


free(Vreal); free(Vimag);
for (int i = 0; i < nBarras; i++) {free(Ireal[i]); free(Iimag[i]);}
}

void calculoDaPotenciaAtivaDeCarga (int nBarras, double* tensao, double** G, double** B, 
                                    double* consumoDaCarga) {
    double somaReal = 0;
    double somaImag = 0;
    for (int i = 0; i < nBarras; i++) {
        somaReal = 0; somaImag = 0;
        for (int j = 0; j < nBarras; j++) {
            somaReal = somaReal + G[i][j];
        }
        consumoDaCarga[i] = 0.001*3*tensao[i]*tensao[i]*somaReal;
    }
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