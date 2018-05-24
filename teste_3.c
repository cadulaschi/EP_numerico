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
    A[0][0] = 2 - (pow(e, x[0]))/((n+1)*(n+1)); A[0][1] = -1;
    for (int i = 1; i < n - 1; i++) {
        for (int j = 0; j < n; j++) {
            if (j == i - 1) A[i][j] = -1;
            else if (j == i) A[i][j] = 2 - ((pow(e,x[i]))/(((n+1)*(n+1))));
            else if (j == i + 1) A[i][j] = -1;
        }
    }
    A[n-1][n-2] = -1; A[n-1][n-1] = 2 - (pow(e, x[n-1]))/((n+1)*(n+1));
}

void constroib (double* b, double *x, int n) {
    b[0] = -(2*x[0] - x[1] - (pow(e, x[0]))/((n+1)*(n+1))); 
    for (int i = 1; i < n - 1; i++) {
      b[i] = -(-x[i-1] + 2*x[i] - x[i+1] -((pow(e, x[i]))/((n+1)*(n+1))));
    }
    b[n-1] = -(-x[n-2] + 2*x[n-1] - (pow(e, x[n-1]))/((n+1)*(n+1)));   
}

int main () {
int n = 40;
n = n - 1;
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
double maxc;

    printf ("\nValores de iniciais utilizados: aproximação inicial nula");
    int iter = 0;
    do {
        iter++;
        constroiJ (A, x, n);
        constroib (b, x, n);
        
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

        maxc = fabs(c[0]);

        for (int i = 1; i < n; i++) {
            if (fabs(c[i]) > fabs(maxc)) maxc = fabs(c[i]);
        }
    } while (maxc > 1e-6);


    printf ("\n\nNumero de iteracoes: %d", iter);

     printf ("\n\nSolucao do sistema para n = %d\n", n + 1);
    for (int i = 0; i < n; i++) {
        printf("%e ", x[i]);
        printf ("\n");
        }
    printf ("\n");

    printf ("\n\n\n\n\n\n");   

    free(A); free(temp); free(b); free(x); free(p); free(c); 
return 0;
}