/******* Algoritmo de decomposição LU para o cálculo 
******** de método de Newton com várias variáveis
******** date: 12/05/2018
*/

#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#include "funcoes.h"


#define tensaoNom1 = 132790.56;
#define tensaoNom21 = 7967.434;
#define tensaoNom22 = 127.017;
#define tensaoNom3 =  7967.434;
#define tensaoNom41 = 7967.434;
#define tensaoNom42 = 254.034;
#define tensaoNom43 = 219.393;


int main () {
/*int n = 40;
n = n - 1;
double** A = (double **)calloc(n, sizeof(double*));
double* b = (double *)calloc(n, sizeof(double));
double* x = (double *)calloc(n, sizeof(double));

metodoDeNewton (A, x, b, n);
*/  
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
        printf("%d %d %.3lf %.3lf %.3lf\n", int1, tipoDaBarra[i], tensao[i],
                                            campo4[i], campo5[i]);
    }
    printf ("\nPQ: %d  ; PV: %d ; Swing: %d\n", nPQ, nPV, nS);

    //Consegui um vetor de tensões iniciais como as tensões PQ e PV
    //nas posições corretas (VPQ1 VPQ2......VPQn VPV1 VPV2......VPVn VS1 VS2 .... VSn)
    //double* vetorV = (double *)calloc(nBarras, sizeof(double));


    //CRIAMOS O VETOR TETA DE TAMANHO NBARRAS 
    //COMO USAMOS CALLOC, ELE SEMPRE SERÁ INICIADO EM ZERO EM TODOS OS VALORES
    double* vetorTeta = (double *)calloc(nBarras, sizeof(double)); 
                        //testes ok
    //jS = nPV + nPQ;
    for(int i = 0; i < nBarras; i ++ /*jS++*/) {
        if (tipoDaBarra[i] == 2) vetorTeta[i] = campo5[i];
        //else vetorTeta[i] = rand();
        //else jS--;
    }
    //TESTE
    printf ("\nPQ: %d  ; PV: %d ; Swing: %d\n", nPQ, nPV, nS);

    //TESTE
    printf("Tensoes\n");
    for (int i = 0; i < nBarras; i++) {
        printf ("%.5lf ", tensao[i]);
    }
    /*
    printf("\n");
    for (int i = 0; i < nBarras; i++) {
        printf ("%.5lf ", vetorV[i]);
    }*/
    printf("\n");
    printf("Tetas\n");
    for (int i = 0; i < nBarras; i++) {
        printf ("%.5lf ", vetorTeta[i]);
    }
    printf("\n");

    fclose(arqDadosBarra);
	
    //printf("\n");
  //YNodal

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
        printf("\n\n Matriz B");
        verMatriz(B, nBarras);
        printf("\n\n Matriz G");
        verMatriz(G, nBarras);

    fclose(arqYNodal);

    double* Pcalc = (double *)calloc(nBarras, sizeof(double));
    calculaPcalc(nBarras, tensao, vetorTeta, G, B, Pcalc);
    printf ("Pcalc antes:     \n");
    verVetordouble(Pcalc, nBarras);
    
    double* PcalcOrdenado = (double *)calloc(nBarras, sizeof(double));
    ordenaVetor(nPQ, nPV, nBarras, tipoDaBarra, PcalcOrdenado , Pcalc);
    printf ("Pcalc depois:     \n");
    verVetordouble(PcalcOrdenado, nBarras);

    double* Qcalc = (double *)calloc(nBarras, sizeof(double));
    calculaQcalc(nBarras, tensao, vetorTeta, G, B, Qcalc);
    printf ("Qcalc antes:     \n");
    verVetordouble(Qcalc, nBarras);

    double* QcalcOrdenado = (double *)calloc(nBarras, sizeof(double));
    ordenaVetor(nPQ, nPV, nBarras, tipoDaBarra, QcalcOrdenado , Qcalc);
    printf ("Qcalc depois:     \n");
    verVetordouble(QcalcOrdenado, nBarras);

    int* troca = (int *)calloc(nBarras, sizeof(int));
    constroiTroca (nPQ, nPV, nBarras, tipoDaBarra, troca);
    printf("Troca: \n");
    verVetorint(troca, nBarras);

    double* desvioP = (double *)calloc(2*nPQ + nPV, sizeof(double));
    double* campo4Ordenado = (double *)calloc(nBarras, sizeof(double));
    ordenaVetor(nPQ, nPV, nBarras, tipoDaBarra, campo4Ordenado, campo4);
    constroiDesvioDePotencia (nPQ, nPV, desvioP, PcalcOrdenado, QcalcOrdenado, campo4Ordenado);
    printf("Vetor desvios: \n");
    verVetordouble(desvioP, 2*nPQ + nPV);

    
    double** J = (double **)calloc(2*nPQ + nPV, sizeof(double *));
    for (int i = 0; i < 2*nPQ + nPV; i++){//nLinhas =nColunas
        J[i] = (double *)calloc(2*nPQ + nPV, sizeof(double));
    }

    constroiJacobiana (nBarras, nPQ, nPV, tensao, vetorTeta, tipoDaBarra, Pcalc, Qcalc, G, B, J);
    
    printf("Jacobiana\n");
    verMatriz(J, 2*nPQ+nPV);
    /*double* desvioPotencia = (double *)calloc(2*nPQ + nPV, sizeof(double));
    double*  = (double *)calloc(2*nPQ + nPV, sizeof(double));
    calculaDesvioDePotencia ();
    */

//free(A); free(x); free(b);
//free(tipoDaBarra); free(tensao);
free(tipoDaBarra); free(tensao); free(campo4); free(campo5); 
for (int i = 0; i < nBarras; i++) {free(B[i]); free(G[i]);}
free(B); free(G);
//free(vetorV); 
free(vetorTeta);
free(Pcalc); free(PcalcOrdenado); 
free(Qcalc);


//free(vetorV); free(vetorVPQ); free(vetorVPQ);

return 0;
}