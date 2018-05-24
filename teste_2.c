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
    A[0][0] = 4 - x[3]; A[0][1] = -1;       A[0][2] = 1;        A[0][3] = -x[0];
    A[1][0] = -1;       A[1][1] = 3 - x[3]; A[1][2] = -2;       A[1][3] = -x[1];
    A[2][0] = 1;        A[2][1] = -2;       A[2][2] = 3 - x[3]; A[2][3] = -x[2];
    A[3][0] = 2*x[0];   A[3][1] = 2*x[1];   A[3][2] = 2*x[2];   A[3][3] = 0; 
}

void constroib (double* b, double *x, int n) {
    b[0] = -(4*x[0] - x[1] + x[2] - x[0]*x[3]);
    b[1] = -(-x[0] + 3*x[1] - 2*x[2] - x[1]*x[3]); 
    b[2] = -(x[0] - 2*x[1] + 3*x[2] - x[2]*x[3]);
    b[3] = -(x[0]*x[0] + x[1]*x[1] + x[2]*x[2] - 1);
}

int main () {
int n = 4;
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
double maxc;    

    x[0] = 1;
    x[1] = 1;
    x[2] = 1;
    x[3] = 1;  

    printf ("\n\n\n\nValores de iniciais utilizados: ");
    for (int i = 0; i < n; i++) printf ("| %.1lf ", x[i]);
        printf ("|\n\n");
    
    int iter = 0;
    do {
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

        maxc = fabs(c[0]);

        for (int i = 1; i < n; i++) {
            if (fabs(c[i]) > fabs(maxc)) maxc = fabs(c[i]);
        }

        printf ("Solucao do sistema %d:\n", iter);
        for (int i = 0; i < n; i++) printf("          %.8lf \n", x[i]);
        printf ("\n");

    } while (maxc > 1e-5);
    printf ("Numero de iteracoes: %d\n", iter);
    printf ("\n\n\n\n\n\n");  

    free(A); free(temp); free(b); free(x); free(p); free(c); 
return 0;
}