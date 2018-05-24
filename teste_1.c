/******* Algoritmo de decomposição LU para o cálculo 
******** de método de Newton com várias variáveis
******** date: 12/05/2018
*/

#include <stdio.h>
#include <stdlib.h> 
#include <math.h>

#define e 2.71828182846

/*FUNCTIONS*/
/*IMPRIME UMA MATRIZ*/
void ver (double **ver, int n) {
    printf ("\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf ("%f  ", ver[i][j]);
        }
        printf ("\n");
    }
    printf("\n");
}

/*FAZ A DECOMPOSICAO LU DE UMA MATRIZ*/
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

void constroiJ (double** A, double* x, int n) {
    A[0][0] = 2; A[0][1] = 0;
    A[1][0] = 0; A[1][1] = 2;
}

void constroib (double* b, double *x, int n) {
   b[0] = -(2*x[0] - 4);
   b[1] = -(2*x[1] - 6); 
}

int main () {
int n = 2;
//n = n - 1;
double** A = (double **)calloc(n, sizeof(double*));
double** temp = (double **)calloc(n, sizeof(double*));
for (int i = 0; i < n; i++) {
    A[i] = (double *)calloc(n, sizeof(double*));
    temp[i] = (double *)calloc(n, sizeof(double*));
}

double* b = (double *)calloc(n, sizeof(double));
double* x = (double *)calloc(n, sizeof(double));
//double* xtemp = (double *)calloc(n, sizeof(double));
double* c = (double *)calloc(n, sizeof(double));
int* p = (int *)calloc(n, sizeof(int));
   

    printf ("\n\n\nValores de iniciais utilizados: \n");
    x[0] = rand()*7/3; 
    x[1] = rand()*2/7;
    printf ("x0: %5lf \n", x[0]);
    printf ("y0: %5lf \n", x[1]);
    printf ("\n");
    
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


    }
    printf ("Solucao do sistema \n");
    printf("x: %5.3lf \n", x[0]);
    printf("y: %5.3lf ", x[1]);
    printf ("\n\n");
    printf ("Numero de iteracoes: %d\n", iter);
    printf ("\n\n\n\n\n\n");

    free(A); free(temp); free(b); free(x); free(p); free(c); 
return 0;
}