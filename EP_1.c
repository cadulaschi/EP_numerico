/******* Algoritmo de decomposição LU para o cálculo 
******** de método de Newton com várias variáveis
******** date: 12/05/2018
*/

#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#include "funcoes.h"

int main () {  
    int pergunta;

    //PEGANDO O TIPO DE REDE
    printf("Escolha a rede para analise: \n");
    printf("1 = Stevenson, 2 = Reticulad, 3 = Distribuicao Primaria ou 4 = Distribuicao primaria e secundaria\n");
	scanf ("%d", &pergunta);


    //DadosBarras
    int nBarras;
    int int1;
    int nPQ = 0;
    int nPV = 0;
    int nS = 0;

    //Vetores criados
 
	FILE *arqDadosBarra;

    //ABRINDO ARQUIVO DA REDE ESCOLIDA
	if (pergunta == 1) arqDadosBarra = fopen("1_Stevenson_DadosBarras.txt", "r");
	else if (pergunta == 2) arqDadosBarra = fopen("2_Reticulada_DadosBarras.txt", "r");
	else if (pergunta == 3) arqDadosBarra = fopen("3_Distribuicao_Primaria_DadosBarras.txt", "r");
    else if (pergunta == 4) arqDadosBarra = fopen("4_Distribuicao_Primaria_Secundaria_DadosBarras.txt", "r");
    printf("\n");

	if(arqDadosBarra == NULL) {
			printf("Erro, nao foi possivel abrir o arquivo DadosBarras do tipo %d\n", pergunta);
            return 0;
    }
        fscanf(arqDadosBarra,"%d\n",&nBarras);

    int* tipoDaBarra = (int *)calloc(nBarras, sizeof(int));
    double* tensao = (double *)calloc(nBarras, sizeof(double));
    double* tensaoNominal = (double *)calloc(nBarras, sizeof(double));
    double* campo4 = (double *)calloc(nBarras, sizeof(double));
    double* campo5 = (double *)calloc(nBarras, sizeof(double));



    //GERANDO AS MATRIZES COM ALOCAÇÃO DINAMICA
    double** B = (double **)calloc(nBarras, sizeof(double *));
    double** G = (double **)calloc(nBarras, sizeof(double *));
    for (int i = 0; i < nBarras; i++){//nLinhas =nColunas
        B[i] = (double *)calloc(nBarras, sizeof(double));
        G[i] = (double *)calloc(nBarras, sizeof(double));
    }

        
    for (int i=0; i < nBarras; i++){
        fscanf(arqDadosBarra,"%d %d %lf %lf %lf\n", &int1, &tipoDaBarra[i], &tensao[i],
                                                    &campo4[i], &campo5[i]);
        if(tipoDaBarra[i]==0) nPQ ++;
        if(tipoDaBarra[i]==1) nPV ++;
        if(tipoDaBarra[i]==2) nS ++;
        //TESTES OK
        //printf("%d %d %.3lf %.3lf %.3lf\n", int1, tipoDaBarra[i], tensao[i],
        //                                    campo4[i], campo5[i]);
    }
    //printf ("\nPQ: %d  ; PV: %d ; Swing: %d\n", nPQ, nPV, nS);

    for (int i = 0; i < nBarras; i++) tensaoNominal[i] = tensao[i];


    //CRIAMOS O VETOR TETA DE TAMANHO NBARRAS 
    //COMO USAMOS CALLOC, ELE SEMPRE SERÁ INICIADO EM ZERO EM TODOS OS VALORES
    double* vetorTeta = (double *)calloc(nBarras, sizeof(double)); 
    
    for(int i = 0; i < nBarras; i ++ /*jS++*/) {
        if (tipoDaBarra[i] == 2) vetorTeta[i] = (PI*campo5[i])/180;
    }
    fclose(arqDadosBarra);

    FILE *arqYNodal;

    if (pergunta == 1) arqYNodal = fopen("1_Stevenson_Ynodal.txt", "r");
	else if (pergunta == 2) arqYNodal = fopen("2_Reticulada_Ynodal.txt", "r");
	else if (pergunta == 3) arqYNodal = fopen("3_Distribuicao_Primaria_Ynodal.txt", "r");
    else if (pergunta == 4) arqYNodal = fopen("4_Distribuicao_Primaria_Secundaria_Ynodal.txt", "r");

    int nLinhasY;
    int nDaLinhaYNodal;
    int nDaColunaYNodal;


	if(arqYNodal == NULL) {
			printf("Erro, nao foi possivel abrir o arquivo do tipo %d\n",pergunta);
            return 0;
    }
        fscanf(arqYNodal,"%d\n",&nLinhasY);
        //ERRO
        if (B == NULL || G == NULL ){
            printf("malloc B ou G devolveu NULL!");
            exit (EXIT_FAILURE);
        }

        //GRAVANDO NAS MATRIZES
        for (int i=0; i<nLinhasY; i++){
            fscanf(arqYNodal, "%d", &nDaLinhaYNodal);
            fscanf(arqYNodal, "%d", &nDaColunaYNodal);
            fscanf(arqYNodal, "%lf", &G[nDaLinhaYNodal][nDaColunaYNodal]);
            fscanf(arqYNodal, "%lf\n", &B[nDaLinhaYNodal][nDaColunaYNodal]);
        }
    
    fclose(arqYNodal);

    int* troca = (int *)calloc(nBarras, sizeof(int));
    constroiTroca (nPQ, nPV, nBarras, tipoDaBarra, troca);  

    printf ("\n START \n");

    metodoDeNewton(nBarras, nPQ, nPV, nS, tensao, vetorTeta, tipoDaBarra, campo4, G, B, troca); 

    double** Sreal = (double **)calloc(nBarras, sizeof(double *));
    double** Simag = (double **)calloc(nBarras, sizeof(double *));
    double** perdaAtiva = (double **)calloc(nBarras, sizeof(double *));
    for (int i = 0; i < nBarras; i++){//nLinhas =nColunas
        Sreal[i] = (double *)calloc(nBarras, sizeof(double));
        Simag[i] = (double *)calloc(nBarras, sizeof(double));
        perdaAtiva[i] = (double *)calloc(nBarras, sizeof(double));
    }


    calculoDa_PotenciaAtiva_e_PerdaAtiva (nBarras, nPQ, nPV, tensao, vetorTeta, tipoDaBarra, G, B, Sreal, 
                                          Simag, perdaAtiva);

    
    double* consumoDaCarga = (double *)calloc(nBarras, sizeof(double));    
    double potenciaAtivaTotalDaCarga = 0;
    
    calculoDaPotenciaAtivaDeCarga (nBarras, tensao, G, B, consumoDaCarga);
    for(int i = 0; i < nBarras; i++) potenciaAtivaTotalDaCarga = potenciaAtivaTotalDaCarga + consumoDaCarga[i];
    //printf("Potência ativa total de carga (absorvida):  %.6lf \n\n\n", potenciaAtivaTotalDaCarga);
    
    double perdaAtivaTotal = 0;

    printf ("Barra          Tensao (pu)          Angulo          Módulo da Tensao\n");
    if (pergunta == 1) {
        for (int i = 0; i < nBarras; i++) {
            printf("%3d", i);
            printf ("%21.6lf", tensao[i]/tensaoNominal[i]);
            printf("%18.4lf", (vetorTeta[i]*180)/PI);
            printf("%23.3lf\n", tensao[i]);
        }

        printf ("\n         Trecho                   Potencia Ativa            Perda Ativa\n");
        printf ("Barra inicial  Barra final\n");
        for (int i = 0; i < nBarras; i++) {
            for (int j = 0; j < nBarras; j++) {
                if((i == 0 && j == 1) || (i == 0 && j == 4) || (i == 1 && j == 2) || (i == 2 && j == 3)
                   || (i == 2 && j == 4) || (i == 3 && j == 4) ) {
                printf ("%7d%14d", i, j);
                printf ("%25.3lf", Sreal[i][j]);
                printf ("%24.3lf", perdaAtiva[i][j]);
                printf("\n");
                }
            if (i < j) perdaAtivaTotal = perdaAtivaTotal + perdaAtiva[i][j];
            }
        }
        printf ("\n");
        printf("Potencia ativa total gerada: %30.3lf\n", perdaAtivaTotal + potenciaAtivaTotalDaCarga);
        printf("Potencia ativa total de carga (absorvida): %16.3lf\n", potenciaAtivaTotalDaCarga);
        printf("Perda ativa total: %40.3lf\n", perdaAtivaTotal);
        printf ("\n\n");
    }

    else if (pergunta == 2) {
        for (int i = 0; i < nBarras; i++) {
            if (i == 2 || i == 11 || i == 25 || i == 28 || i == 30 || i == 42 || i == 43 || i == 47
                || i == 48 || i == 49) {
             printf("%3d", i);
            printf ("%21.6lf", tensao[i]/tensaoNominal[i]);
            printf("%18.4lf", (vetorTeta[i]*180)/PI);
            printf("%23.3lf\n", tensao[i]);
            }
        }
        printf ("\n         Trecho                   Potencia Ativa            Perda Ativa\n");
        printf ("Barra inicial  Barra final\n");
        for (int i = 0; i < nBarras; i++) {
            for (int j = 0; j < nBarras; j++) {
                if((i == 3 && j == 4) || (i == 6 && j == 5) || (i == 12 && j == 28) || (i == 13 && j == 9)
                   || (i == 17 && j == 1) || (i == 18 && j == 2) || (i == 19 && j == 20) 
                   || (i == 24 && j == 52) || (i == 60 && j == 62) || (i == 75 && j == 2))  {
                printf ("%7d%14d", i, j);
                printf ("%25.3lf", Sreal[i][j]);
                printf ("%24.3lf", perdaAtiva[i][j]);
                printf("\n");
                }
                if (i < j) perdaAtivaTotal = perdaAtivaTotal + perdaAtiva[i][j];     
            }
        }
        printf ("\n");
        printf("Potencia ativa total gerada: %30.3lf\n", perdaAtivaTotal + potenciaAtivaTotalDaCarga);
        printf("Potencia ativa total de carga (absorvida): %16.3lf\n", potenciaAtivaTotalDaCarga);
        printf("Perda ativa total: %40.3lf\n", perdaAtivaTotal);
        printf ("\n\n");
    }

    else if (pergunta == 3) {
        for (int i = 0; i < nBarras; i++) {
            if (i == 0 || i == 1 || i == 47 || i == 633 || i == 1414 || i == 1429 || i == 1528 || i == 1607
                || i == 1609 || i == 1636) {
             printf("%6d", i);
            printf ("%18.6lf", tensao[i]/tensaoNominal[i]);
            printf("%18.4lf", (vetorTeta[i]*180)/PI);
            printf("%23.3lf\n", tensao[i]);
            }
        }
        
        printf ("\n         Trecho                   Potencia Ativa            Perda Ativa\n");
        printf ("Barra inicial  Barra final\n");
        for (int i = 0; i < nBarras; i++) {
            for (int j = 0; j < nBarras; j++) {
                if((i == 0 && j == 1185) || (i == 1 && j == 2) || (i == 1 && j == 92) || (i == 47 && j == 6)
                   || (i == 47 && j == 31) || (i == 633 && j == 632) || (i == 633 && j == 634) 
                   || (i == 1414 && j == 1415) || (i == 1607 && j == 286) || (i == 1621 && j == 1622))  {
                printf ("%7d%14d", i, j);
                printf ("%25.3lf", Sreal[i][j]);
                printf ("%24.3lf", perdaAtiva[i][j]);
                printf("\n");
                }
                if (i < j) perdaAtivaTotal = perdaAtivaTotal + perdaAtiva[i][j];     
            }
        }
        printf ("\n");
        printf("Potencia ativa total gerada: %30.3lf\n", perdaAtivaTotal + potenciaAtivaTotalDaCarga);
        printf("Potencia ativa total de carga (absorvida): %16.3lf\n", potenciaAtivaTotalDaCarga);
        printf("Perda ativa total: %40.3lf\n", perdaAtivaTotal);
        printf ("\n\n");
    }

    else if (pergunta == 4) {
        for (int i = 0; i < nBarras; i++) {
            if (i == 3 || i == 990 || i == 1310 || i == 1466 || i == 3947 || i == 4015 || i == 4188 
                || i == 5820 || i == 5830 || i == 5840) {
            printf("%6d", i);
            printf ("%18.6lf", tensao[i]/tensaoNominal[i]);
            printf("%18.4lf", (vetorTeta[i]*180)/PI);
            printf("%23.3lf\n", tensao[i]);
            }
        }

        printf ("\n         Trecho                   Potencia Ativa            Perda Ativa\n");
        printf ("Barra inicial  Barra final\n");
        for (int i = 0; i < nBarras; i++) {
            for (int j = 0; j < nBarras; j++) {
                if((i == 0 && j == 1185) || (i == 710 && j == 543) || (i == 776 && j == 1748) 
                   ||(i == 1543 && j == 1542) || (i == 1600 && j == 1387) || (i == 1631 && j == 1630) 
                   ||(i == 1748 && j == 776)  || (i == 2867 && j == 2868) || (i == 2878 && j == 2877) 
                   || (i == 3640 && j == 3947))  {
                printf ("%7d%14d", i, j);
                printf ("%25.3lf", Sreal[i][j]);
                printf ("%24.3lf", perdaAtiva[i][j]);
                printf("\n");
                }
                if (i < j) perdaAtivaTotal = perdaAtivaTotal + perdaAtiva[i][j];     
            }
        }
        printf ("\n");
        printf("Potencia ativa total gerada: %30.3lf\n", perdaAtivaTotal + potenciaAtivaTotalDaCarga);
        printf("Potencia ativa total de carga (absorvida): %16.3lf\n", potenciaAtivaTotalDaCarga);
        printf("Perda ativa total: %40.3lf\n", perdaAtivaTotal);
        printf ("\n\n");
    }

free(tipoDaBarra); free(tensao); free(campo4); free(campo5); 
for (int i = 0; i < nBarras; i++) {
    free(B[i]); free(G[i]);
    free(Sreal[i]);
    free(Simag[i]);
    free(perdaAtiva[i]);
}
free(Sreal);
free(Simag);
free(perdaAtiva);
free(troca);
free(B); free(G);
free(tensaoNominal);
free(vetorTeta);
return 0;
}