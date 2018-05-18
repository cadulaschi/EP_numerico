void metodoDeNewton (double** A, double* x, double* b, int n) {
double** temp = (double **)calloc(n, sizeof(double*));
for (int i = 0; i < n; i++) {
    temp[i] = (double *)calloc(n, sizeof(double*));
}

double* c = (double *)calloc(n, sizeof(double));
int* p = (int *)calloc(n, sizeof(int));

    printf ("\n\n\nValores de iniciais utilizados: ");
    for (int i = 0; i < n; i++) printf ("%.1lf | ", x[i]);
        printf ("\n\n");
    
    int iter = 0;
    while (iter < 10) {
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

        for (int i = 0; i < n; i++) x[i] = x[i] + c[i];

        printf ("Solucao do sistema %d:\n", iter);
        for (int i = 0; i < n; i++) printf("          %.8lf \n", x[i]);
        printf ("\n");
    }
    printf ("Numero de iteracoes: %d\n", iter);
    
    free(temp); free(c); free(p); free(A);
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