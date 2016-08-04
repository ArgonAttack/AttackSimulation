
//  Change Log:
// Modified construction of set S to save factor of 2 
//  (key notion: advancing in segment)
// Minor: don't need pebble at begining of layer
// Call balloon phase at curlayer*rootd steps before the end of a light phase (we don't need d steps until the end when curlayer = rootd)
// moved core helper functions (generating graphs, generating sets S etc...) to BuildGraph.c
// added core function int * selectPbSetS( int n, struct bnode* G, int parallelism)
// modified function templates to take pair (sizeOfLayer, gapInLayer) instead of rootd --- previously when we had no bound on parallelism both values could be derived from rootd

#include "Attack.h"
#include "BuildGraph.h"
#include "statistics.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define DEBUG
#ifdef DEBUG
# define DEBUG_PRINT(fmt, args)    printf(fmt, ## args)
#else
# define DEBUG_PRINT(x) do {} while (0)
#endif



// Forward Declarations
void testParStats(int exponent, int passes, double R, int n, int window, int parUbound);
void testSearchbForSegSize(int exponent, int passes, double R, int n, int window, int parUBound, int g);
int testSearchbForg(int exponent, int passes, double R, int n, int window, int parUBoundd);
void testSearchbFord(int exponent, int passes, double R, int n, int window);
void testSelectSetS(int exponent, int passes, double R, int n, int window, int parUBound);
void testRecursiveAttack(int exponent, int passes, double R, int n, int window, int rootd1, int rootd2, int g1, int g2);

void generatePlot(int *nPows, int nLength, int *passes, int passesLength, int numSamples, FILE * fp, FILE * fp2, double R);
void generatePlotOriginal(int *nPows, int nLength, int *passes, int passesLength, int numSamples, FILE * fp, FILE * fp2, double R);
void generatePlotBalloon(int *nPows, int nLength, int *passes, int passesLength, int numSamples, FILE * fp, FILE * fp2, double R);
void generateParallelismPlot(int nPow, int *parBounds, int parBoundsLength, int *passes, int passesLength, int numSamples, FILE *fp, FILE *fp2,  double R);
void generateParallelismPlotFixedMemory(int mPow, int *parBounds, int parBoundsLength, int *passes, int passesLength, int numSamples, FILE *fp, FILE *fp2,  double R, int type);
void generatePlotArgon4i(int *nPows,int nLength, int *passesArray, int passesLength, int samples, FILE* fp, FILE* fp2,double R);
void generatePlotArgon4iWithParallelism(int * nPows,int nLength, int * passesArray, int passesLength, int numSamples, FILE *fp, FILE *fp2, double R,int parallelism);
void generatePlotBalloonHighMemory(int *nPows, int nLength, int *passes, int passesLength, int numSamples, FILE * fp, FILE * fp2, double R);

// Main
int main(int argc, char* argv[])
{

    if (argc > 1) printf("Ignoring arguments %s,....\n", argv[1]);
	// Instance Parameters

	double R = 3000.0;   // core-memory energy ratio

    FILE * fp;
    FILE * fp2;


    fp = fopen ("results.txt", "w+");
    fp2 = fopen("params.txt", "w+");
    int nOverP = pow(2,20);
    int window = pow(2,21);
    int parallelism = 8;
    int numLayers = 128;
    int gap = 32;
    int g = pow(nOverP,0.85);

    printf("Honest Cost: (nodes = %d, window = %d, parallelism = %d) \n  ---- cost = %f\n", nOverP*parallelism, window,parallelism, honestECompParallelism(nOverP, window, R, parallelism));
    printf("Honest Cost: (nodes = %d, window = %d, parallelism = %d) \n  ---- cost = %f\n", nOverP*parallelism, window/2, parallelism, honestECompParallelism(nOverP, window/2, R, parallelism));
    printf("Honest Cost: (nodes = %d, window = %d, parallelism = %d) \n  ---- cost = %f\n", nOverP*parallelism, window/4, parallelism, honestECompParallelism(nOverP, window/4, R, parallelism));


    struct bnode* G = generateArgon2iBGraphParallelism(nOverP,window,123456,parallelism);
    int * S =  selectbSetGapSParallelism(gap, numLayers, nOverP, parallelism, G);
    double attackCost = attackParallelLanes(G,S,nOverP/256,gap,nOverP, g,window,R,parallelism);



    printf("Parallelism: %d, ", parallelism); printf("nodes %d, memory %d\n", nOverP*parallelism, window);
    printf("Honest Cost: %f \n" ,honestECompParallelism(nOverP, window, R, parallelism));
    printf("Attack Cost: %f \n,", attackCost);
    printf("Attack Quality: %f \n",honestECompParallelism(nOverP, window, R, parallelism)/attackCost);
    int nPows[8] = {17,18,19,20,21,22,23,24};
    int nLength = 8;
    int passesArray[6] =  {1,3,4,6,10,20}; //  {2,3};
    int passesLength = 6;
    generatePlotArgon4iWithParallelism(nPows,nLength, passesArray, passesLength, 10, fp, fp2, R,parallelism);
   // generatePlotOriginal(nPows, nLength, passesArray, passesLength, 10, fp,fp2,R );
   // generatePlotBalloon(nPows, nLength, passesArray, passesLength, 10, fp,fp2,R );
   // generatePlotBalloonHighMemory(nPows, nLength, passesArray, passesLength, 10, fp,fp2,R );

    int parBoundsPasses[2/*4*/]={/*1,4,*/7,10};
    int parArrayPassesLength = 2;//4;
    int parBounds[7] = {25,50,100,200,500,750, 1000};
    int parBoundsLength = 7;
    int tempPar[3] = {500,750,1000};
    int temp4Pass[1] = {4};
    generateParallelismPlotFixedMemory(20,tempPar, 3,temp4Pass,1,10,fp,fp2,R, 1);
    generateParallelismPlotFixedMemory(20,parBounds, parBoundsLength,parBoundsPasses,parArrayPassesLength,10,fp,fp2,R, 1);

   // generateParallelismPlot(20,parBounds, parBoundsLength,parBoundsPasses,parArrayPassesLength,10,fp,fp2,R);
    fclose(fp); fclose(fp2);






	getchar();

	return 0;
}


void testParStats(int exponent, int passes, double R, int n, int window, int parUBound)
{
	// Sample Argon3i graph and choose g and parallelism upper-bound
	struct bnode* G = generateBHGraph(n, window, 123456);
	parStats(G, n, window, R, parUBound);
}

void testSearchbForSegSize(int exponent, int passes, double R, int n, int window, int parUBound, int g)
{
	// Sample Argon3i graph and choose g and parallelism upper-bound
	struct bnode* G = generateBHGraph(n, window, 123456);
	//int g = 33847;
	int granularity = 5;   // number of iterations to search for the optimal g.
	//selectPbSetS(n, G, 340, 19066/56);
	//selectbSetS(55, n, G);
	struct sSet *S = searchbForSegSize(G, n, window, R, g, parUBound);
	searchbForg(G, S->S, n, window, R, S->sizeOfLayer, S->segSize, granularity);
	free(S);
	freeGraph(G, n);
}

int testSearchbForg(int exponent, int passes, double R, int n, int window, int parUBound)
{
    // Parameters
	int granularity = 5;   // number of iterations to search for the optimal g.

	struct bnode* G = generateBHGraph(n, window, 123456);   // sample Argon3i graph
	int segSize = maximum(1,(int)(n/(pow(parUBound,1.5)*1.1)));        // optimal #layers is a little bit bigger than sqrt(p) whenever p^2 << n
	struct sSet* S = selectPbSetS(n, G, parUBound, segSize);

	printf("Running searchbForg with segSize = %d\t sizeOfLayer = %d\n", S->segSize, S->sizeOfLayer);
	// Search for optimal g
	int g= searchbForg(G, S->S, n, window, R, S->sizeOfLayer, S->segSize, granularity);
	free(S);
	freeGraph(G,n);
	return g;
}

void testSearchbFord(int exponent, int passes, double R, int n, int window)
{
	int rootd = (int)pow(2, exponent / 4.0);

	// Sample Argon3i graph
	struct bnode* G = generateBHGraph(n, window, 123456);

	//Find a node set S cutting depth down to d = rootd^2
	int* S = selectbSetS(rootd, n, G);

	int g = searchbForg(G, S, n, window, R, (int)ceil(n*1.0 / (rootd*1.0)), rootd, 5); // 33847; 
	free(S);
	searchbFord(G, &S, n, window, R, g);

	freeGraph(G,n);
}



void testSelectSetS(int exponent, int passes, double R, int n, int window, int parUBound)
{
	int rootd = (int)pow(2, exponent / 4.0);
	int segSize = rootd;

	// Sample Argon3i graph
	struct bnode* G = generateBHGraph(n, window, 123456);

	//Find a node set S cutting depth down to d = rootd^2
	int* S = selectbSetS(rootd, n, G);
	struct sSet *sSet = selectPbSetS(n, G, (int)1024, segSize);

	free(S);
	freeGraph(G, n);
}

// int type = 0 -- Argon2iA with uniform edge distribution (original version from Password Hashing Competition)
// int type > 1 --- Argon2iB with non-uniform edge distribution and XOR blocks before writing into memory (latest version)
void generateParallelismPlotFixedMemory(int mPow, int *parBounds, int parBoundsLength, int *passes, int passesLength, int numSamples, FILE *fp, FILE *fp2,  double R, int type)
{
   fprintf(fp, "\\%%M=2^%d\n", mPow);
   fprintf(fp2, "\\%%M=2^%d\n", mPow);
   int window = (int)pow(2,mPow);
   int i = 0;
   for(i=0; i < passesLength; i++)
   {
	   int numPasses= passes[i];
	   int n = window*numPasses;


	   fprintf(fp,"\\%% passes =  %d \n", numPasses);
	   fprintf(fp2,"\\%% passes =  %d \n", numPasses);
       int j = 0;

       for( j=0; j < parBoundsLength; j++)
       {
          int parBound = parBounds[j];

          double averageCost = 0.0;
          int * S = NULL;
          int g, gap, layers, sizeOfLayer;
          int k;
          for( k =0; k < numSamples; k++)
          {
            	int randSeed = i*numSamples+j*passesLength*numSamples+k;
            	struct bnode* G = (type==0)?generateArgon2iAGraph(n, window, randSeed):generateArgon2iBGraph(n, window, randSeed);
				if (k==0){
					struct sSet* data = NULL;
					g = (int)pow(n,0.75);

					data = searchbForSegSize(G, n, window, R, g, parBound);

					S = data->S;

					gap = data->segSize;

					layers = data->numLayers;
						//printf("gap %d\n, layers %d\n", gap, layers);

					sizeOfLayer = (int)ceil(n*1.0/(1.0*(layers)));


					g= searchbForg(G, S, n, window, R,  sizeOfLayer, (gap), 5);

					fprintf(fp,"\\addplot coordinates {\n");
					fprintf(fp,"%%n=%d, g = %d, gap=%d, NumLayers=%d, SizeOfLayer =%d \n", n, g, gap, layers, sizeOfLayer);
					fflush(fp);
					fprintf(fp2,"g=%d\n gap=%d\n NumLayers=%d\n SizeOfLayer=%d\n parallelism=%d\n NumLightPhaseUnits=%d\n",g,gap,layers,sizeOfLayer, (int)ceil(sizeOfLayer/(1.0*(gap+1))), (int)ceil(g/(1.0*(gap*layers))));
					fprintf(fp2, "HonestCost: %f\n", honestEComp(n,window,R));
					fflush(fp2);
					free(data->S);
					free(data);
				}
				S = selectbSetGapS(gap,layers,n,G);
				int sSize=0;
				int ii = 0;
                for(ii = 0; ii < n; ii++) sSize+=S[ii];


				double cost = attackBH(G, S, sizeOfLayer, gap, n, g, window,  R);
                averageCost+= cost/(1.0*numSamples);
                fprintf(fp2,"Instance %d: |S|=%d, cost=%f \n",k, sSize, cost);
                fflush(fp2);
					//printf("Cost: %f \n",cost);
					//printf("Honest Cost: %f \n", honestCost);
					//printf("Quality: %f \n" , honestCost/cost);

				    freeGraph(G,n);
				    free(S);
          }
		  double honestCost = honestEComp(n,window,R);
		  double attackQuality = honestCost/averageCost;
          fprintf(fp, " (%d,%f) ", parBound, attackQuality);
       }
       fflush(fp);
       fflush(fp2);
   }
}

void generateParallelismPlot(int nPow, int *parBounds, int parBoundsLength, int *passes, int passesLength, int numSamples, FILE *fp, FILE *fp2,  double R)
{
   fprintf(fp, "\\%%n=2^%d\n", nPow);
   fprintf(fp2, "\\%%n=2^%d\n", nPow);
	int n = (int)pow(2,nPow);

   int i = 0;
   for(i=0; i < passesLength; i++)
   {
	   int numPasses= passes[i];
	   int window = n/numPasses;
	   fprintf(fp,"\\%% passes =  %d \n", numPasses);
	   fprintf(fp2,"\\%% passes =  %d \n", numPasses);
       int j = 0;

       for( j=0; j < parBoundsLength; j++)
       {
          int parBound = parBounds[j];

          double averageCost = 0.0;
          int * S = NULL;
          int g, gap, layers, sizeOfLayer;
          int k;
          for( k =0; k < numSamples; k++)
          {
            	int randSeed = i*numSamples+j*passesLength*numSamples+k;
            	struct bnode* G = generateBHGraph(n, window, randSeed);
				if (k==0){
					struct sSet* data = NULL;
					g = (int)pow(n,0.75);

					data = searchbForSegSize(G, n, window, R, g, parBound);

					S = data->S;

					gap = data->segSize;

					layers = data->numLayers;
						//printf("gap %d\n, layers %d\n", gap, layers);

					sizeOfLayer = (int)ceil(n*1.0/(1.0*(layers)));


					g= searchbForg(G, S, n, window, R,  sizeOfLayer, (gap), 5);

					fprintf(fp,"\\addplot coordinates {\n");
					fprintf(fp,"%% g = %d, gap=%d, NumLayers=%d, SizeOfLayer =%d \n", g, gap, layers, sizeOfLayer);
					fprintf(fp2,"g=%d\n gap=%d\n NumLayers=%d\n SizeOfLayer=%d\n parallelism=%d\n NumLightPhaseUnits=%d\n",g,gap,layers,sizeOfLayer, (int)ceil(sizeOfLayer/(1.0*(gap+1))), (int)ceil(g/(1.0*(gap*layers))));
					fprintf(fp2, "HonestCost: %f\n", honestEComp(n,window,R));
					free(data->S);
					free(data);
				}
				S = selectbSetGapS(gap,layers,n,G);
				int sSize=0;
				int ii = 0;
                for(ii = 0; ii < n; ii++) sSize+=S[ii];


				double cost = attackBH(G, S, sizeOfLayer, gap, n, g, window,  R);
                averageCost+= cost/(1.0*numSamples);
                fprintf(fp2,"Instance %d: |S|=%d, cost=%f \n",k, sSize, cost);
                fflush(fp2);
					//printf("Cost: %f \n",cost);
					//printf("Honest Cost: %f \n", honestCost);
					//printf("Quality: %f \n" , honestCost/cost);

				    freeGraph(G,n);
				    free(S);
          }
		  double honestCost = honestEComp(n,window,R);
		  double attackQuality = honestCost/averageCost;
          fprintf(fp, " (%d,%f) ", parBound, attackQuality);
       }
       fflush(fp);
       fflush(fp2);
   }
}

void generatePlotBalloonHighMemory(int *nPows, int nLength, int *passes, int passesLength, int numSamples, FILE * fp, FILE * fp2, double R)
{
    int i = 0;
    for(i=0; i < passesLength; i++)
    {
    	int numPasses = passes[i];
    	fprintf(fp, "%%  passes = %d, samples =  %d\n", numPasses, numSamples);
    	fprintf(fp,"\\addplot coordinates {");
    	fflush(fp);
        int j = 0;
        for( j=0; j < nLength; j++){
        	int n = (int)pow(2,nPows[j]);
        	int window = n/numPasses;

        	int k = 0;
        	double averageCost = 0.0;
        	int * S = NULL;
        	int g, gap, layers, sizeOfLayer;
        	fprintf(fp2, "\n___________________n=%d, passes = %d\n___________________\n",n,numPasses);
        	double honestCost = honestEComp(n,window,R);
            printf("n=%d\nHonest Cost: %f \n",n,honestCost);

        	for( k=0; k < numSamples; k++)
			{

	        	int randSeed = i*numSamples+j*passesLength*numSamples+k;
	        	printf("Building Graph\n");
	        	struct bnodeCompressed* G = generateBHLinGraphNoMalloc(n, window, randSeed);
                printf("Built Graph\n");
				if (k==0){
					struct searchData* data = NULL;
					g = (int)pow(n,0.75);
					printf("Searching for Best Depth Parameters\n");
					data = searchbForBestDepthParametersCompressed(G, n, window, R, g);

					S = data->S;

					gap = data->gapInLayer;

					layers = data->numLayers;

					printf("Found Them: gap=%d, layers = %d\n", gap, layers);
						//printf("gap %d\n, layers %d\n", gap, layers);

					sizeOfLayer = (int)ceil(n*1.0/(1.0*(layers)));

                    printf("Searching for g\n");
					g= searchbForgCompressed(G, S, n, window, R,  sizeOfLayer, (gap), 5);
					//fprintf(fp,"\\%% g = %d, gap=%d, NumLayers=%d, SizeOfLayer =%d \n", g, gap, layers, sizeOfLayer);

					fprintf(fp2,"g=%d\n gap=%d\n NumLayers=%d\n SizeOfLayer=%d\n parallelism=%d\n NumLightPhaseUnits=%d\n",g,gap,layers,sizeOfLayer, (int)ceil(sizeOfLayer/(1.0*(gap+1))), (int)ceil(g/(1.0*(gap*layers))));
					fprintf(fp2, "HonestCost: %f\n", honestEComp(n,window,R));
					fflush(fp2);
					free(data->S);
					free(data);
				}
				printf("Building S\n");
				S = selectbSetGapSCompressed(gap,layers,n,window, G);
				printf("Built S\n");
				int sSize=0;
				int ii;
                for(ii = 0; ii < n; ii++) sSize+=S[ii];

                printf("Attacking\n");
				double cost = attackiXOR(G, S, sizeOfLayer, gap, n, g, window,  R);
				printf("Attacked\n");
                averageCost+= cost/(1.0*numSamples);
                fprintf(fp2,"Instance %d:  |S|=%d, cost=%f \n",k, sSize, cost);
                fflush(fp2);
					//printf("Cost: %f \n",cost);
					//printf("Honest Cost: %f \n", honestCost);
					//printf("Quality: %f \n" , honestCost/cost);

				    //freeGraph(G,n);
				    free(S);
				    free(G);
			}


			double attackQuality = honestCost/averageCost;
            fprintf(fp, " (%d,%f) ", nPows[j], attackQuality);
            fflush(fp);
        }
        fprintf(fp,"};\n ");
        fflush(fp);
        fflush(fp2);
    }
}

void generatePlotBalloon(int *nPows, int nLength, int *passes, int passesLength, int numSamples, FILE * fp, FILE * fp2, double R)
{
    int i = 0;
    for(i=0; i < passesLength; i++)
    {
    	int numPasses = passes[i];
    	fprintf(fp, "%%  passes = %d, samples =  %d\n", numPasses, numSamples);
    	fprintf(fp,"\\addplot coordinates {");
    	fflush(fp);
        int j = 0;
        for( j=0; j < nLength; j++){
        	int n = (int)pow(2,nPows[j]);
        	int window = n/numPasses;

        	int k = 0;
        	double averageCost = 0.0;
        	int * S = NULL;
        	int g, gap, layers, sizeOfLayer;
        	fprintf(fp2, "\n___________________n=%d, passes = %d\n___________________\n",n,numPasses);

        	for( k=0; k < numSamples; k++)
			{

	        	int randSeed = i*numSamples+j*passesLength*numSamples+k;
	        	struct bnode* G = generateiXORGraph(n, window, randSeed);

				if (k==0){
					struct searchData* data = NULL;
					g = (int)pow(n,0.75);
					data = searchbForBestDepthParameters(G, n, window, R, g);

					S = data->S;

					gap = data->gapInLayer;

					layers = data->numLayers;
						//printf("gap %d\n, layers %d\n", gap, layers);

					sizeOfLayer = (int)ceil(n*1.0/(1.0*(layers)));


					g= searchbForg(G, S, n, window, R,  sizeOfLayer, (gap), 5);
					//fprintf(fp,"\\%% g = %d, gap=%d, NumLayers=%d, SizeOfLayer =%d \n", g, gap, layers, sizeOfLayer);

					fprintf(fp2,"g=%d\n gap=%d\n NumLayers=%d\n SizeOfLayer=%d\n parallelism=%d\n NumLightPhaseUnits=%d\n",g,gap,layers,sizeOfLayer, (int)ceil(sizeOfLayer/(1.0*(gap+1))), (int)ceil(g/(1.0*(gap*layers))));
					fprintf(fp2, "HonestCost: %f\n", honestEComp(n,window,R));
					fflush(fp2);
					free(data->S);
					free(data);
				}
				S = selectbSetGapS(gap,layers,n,G);
				int sSize=0;
				int ii;
                for(ii = 0; ii < n; ii++) sSize+=S[ii];


				double cost = attackBH(G, S, sizeOfLayer, gap, n, g, window,  R);
                averageCost+= cost/(1.0*numSamples);
                fprintf(fp2,"Instance %d:  |S|=%d, cost=%f \n",k, sSize, cost);

					//printf("Cost: %f \n",cost);
					//printf("Honest Cost: %f \n", honestCost);
					//printf("Quality: %f \n" , honestCost/cost);

				    freeGraph(G,n);
				    free(S);
			}

			double honestCost = honestEComp(n,window,R);
			double attackQuality = honestCost/averageCost;
            fprintf(fp, " (%d,%f) ", nPows[j], attackQuality);
            fflush(fp);
        }
        fprintf(fp,"};\n ");
        fflush(fp);
        fflush(fp2);
    }
}

void generatePlot(int *nPows, int nLength, int *passes, int passesLength, int numSamples, FILE * fp, FILE * fp2, double R)
{
    int i = 0;
    for(i=0; i < passesLength; i++)
    {
    	int numPasses = passes[i];
    	fprintf(fp, "\\%%  passes = %d, samples =  %d\n", numPasses, numSamples);
    	fprintf(fp,"\\addplot coordinates {");
        int j = 0;
        for( j=0; j < nLength; j++){
        	int n = (int)pow(2,nPows[j]);
        	int window = n/numPasses;

        	int k = 0;
        	double averageCost = 0.0;
        	int * S = NULL;
        	int g, gap, layers, sizeOfLayer;
        	fprintf(fp2, "\n___________________n=%d, passes = %d\n___________________\n",n,numPasses);

        	for( k=0; k < numSamples; k++)
			{

	        	int randSeed = i*numSamples+j*passesLength*numSamples+k;
	        	struct bnode* G = generateBHGraph(n, window, randSeed);

				if (k==0){
					struct searchData* data = NULL;
					g = (int)pow(n,0.75);
					data = searchbForBestDepthParameters(G, n, window, R, g);

					S = data->S;

					gap = data->gapInLayer;

					layers = data->numLayers;
						//printf("gap %d\n, layers %d\n", gap, layers);

					sizeOfLayer = (int)ceil(n*1.0/(1.0*(layers)));


					g= searchbForg(G, S, n, window, R,  sizeOfLayer, (gap), 5);
					//fprintf(fp,"\\%% g = %d, gap=%d, NumLayers=%d, SizeOfLayer =%d \n", g, gap, layers, sizeOfLayer);

					fprintf(fp2,"g=%d\n gap=%d\n NumLayers=%d\n SizeOfLayer=%d\n parallelism=%d\n NumLightPhaseUnits=%d\n",g,gap,layers,sizeOfLayer, (int)ceil(sizeOfLayer/(1.0*(gap+1))), (int)ceil(g/(1.0*(gap*layers))));
					fprintf(fp2, "HonestCost: %f\n", honestEComp(n,window,R));
					fflush(fp2);
					free(data->S);
					free(data);
				}
				S = selectbSetGapS(gap,layers,n,G);
				int sSize=0;
				int ii;
                for(ii = 0; ii < n; ii++) sSize+=S[ii];


				double cost = attackBH(G, S, sizeOfLayer, gap, n, g, window,  R);
                averageCost+= cost/(1.0*numSamples);
                fprintf(fp2,"Instance %d: |S|=%d, cost=%f \n",k, sSize, cost);

					//printf("Cost: %f \n",cost);
					//printf("Honest Cost: %f \n", honestCost);
					//printf("Quality: %f \n" , honestCost/cost);

				    freeGraph(G,n);
				    free(S);
			}

			double honestCost = honestEComp(n,window,R);
			double attackQuality = honestCost/averageCost;
            fprintf(fp, " (%d,%f) ", nPows[j], attackQuality);
        }
        fprintf(fp,"};\n ");
        fflush(fp);
        fflush(fp2);
    }
}
void generatePlotArgon4iWithParallelism(int * nPows,int nLength, int * passesArray, int passesLength, int numSamples, FILE *fp, FILE *fp2, double R,int parallelism)
{
    int i = 0;
    for(i=0; i < passesLength; i++)
    {
    	int numPasses = passesArray[i];
    	fprintf(fp, "\\%%  passes = %d, samples =  %d\n", numPasses, numSamples);
    	fprintf(fp,"\\addplot coordinates {");
        int j = 0;
        for( j=0; j < nLength; j++){
        	int n = (int)pow(2,nPows[j]);
        	int window = n*parallelism/numPasses;

        	int k = 0;
        	double averageCost = 0.0;
        	int * S = NULL;
        	int g, gap, layers, sizeOfLayer;
        	fprintf(fp2, "\n___________________n=%d, passes = %d\n___________________\n",n,numPasses);

        	for( k=0; k < numSamples; k++)
			{

	        	int randSeed = i*numSamples+j*passesLength*numSamples+k;
	        	struct bnode* G = generateArgon2iBGraphParallelism(n, window, randSeed,parallelism);

				if (k==0){
					struct searchData* data = NULL;
					g = (int)pow(n,0.75);
					data = searchbForBestDepthParameters(G, n, window, R, g);

					S = data->S;

					gap = data->gapInLayer;

					layers = data->numLayers;
						//printf("gap %d\n, layers %d\n", gap, layers);

					sizeOfLayer = (int)ceil(n*1.0/(1.0*(layers)));


					g= searchbForg(G, S, n, window, R,  sizeOfLayer, (gap), 5);
					//fprintf(fp,"\\%% g = %d, gap=%d, NumLayers=%d, SizeOfLayer =%d \n", g, gap, layers, sizeOfLayer);

					fprintf(fp2,"g=%d\n gap=%d\n NumLayers=%d\n SizeOfLayer=%d\n parallelism=%d\n NumLightPhaseUnits=%d\n",g,gap,layers,sizeOfLayer, (int)ceil(sizeOfLayer/(1.0*(gap+1))), (int)ceil(g/(1.0*(gap*layers))));
					fprintf(fp2, "HonestCost: %f\n", honestEComp(n,window,R));
					fflush(fp2);
					free(data->S);
					free(data);
				}
				S = selectbSetGapS(gap,layers,n,G);
				int sSize=0;
				int ii;
                for(ii = 0; ii < n; ii++) sSize+=S[ii];


				double cost = attackBH(G, S, sizeOfLayer, gap, n, g, window,  R);
                averageCost+= cost/(1.0*numSamples);
                fprintf(fp2,"Instance %d: |S|=%d, cost=%f \n",k, sSize, cost);

					//printf("Cost: %f \n",cost);
					//printf("Honest Cost: %f \n", honestCost);
					//printf("Quality: %f \n" , honestCost/cost);

				    freeGraph(G,n);
				    free(S);
			}

			double honestCost = honestEComp(n,window,R);
			double attackQuality = honestCost/averageCost;
            fprintf(fp, " (%d,%f) ", nPows[j], attackQuality);
            fflush(fp);
        }
        fprintf(fp,"};\n ");
        fflush(fp);
        fflush(fp2);
    }
}
void generatePlotArgon4i(int *nPows, int nLength, int *passes, int passesLength, int numSamples, FILE * fp, FILE * fp2, double R)
{
    int i = 0;
    for(i=0; i < passesLength; i++)
    {
    	int numPasses = passes[i];
    	fprintf(fp, "%%  passes = %d, samples =  %d\n", numPasses, numSamples);
    	fprintf(fp,"\\addplot coordinates {");
        int j = 0;
        for( j=0; j < nLength; j++){
        	int n = (int)pow(2,nPows[j]);
        	int window = n/numPasses;

        	int k = 0;
        	double averageCost = 0.0;
        	int * S = NULL;
        	int g, gap, layers, sizeOfLayer;
        	fprintf(fp2, "\n___________________n=%d, passes = %d\n___________________\n",n,numPasses);

        	for( k=0; k < numSamples; k++)
			{

	        	int randSeed = i*numSamples+j*passesLength*numSamples+k;
	        	struct bnode* G = generateArgon2iBGraph(n, window, randSeed);

				if (k==0){
					struct searchData* data = NULL;
					g = (int)pow(n,0.75);
					data = searchbForBestDepthParameters(G, n, window, R, g);

					S = data->S;

					gap = data->gapInLayer;

					layers = data->numLayers;
						//printf("gap %d\n, layers %d\n", gap, layers);

					sizeOfLayer = (int)ceil(n*1.0/(1.0*(layers)));


					g= searchbForg(G, S, n, window, R,  sizeOfLayer, (gap), 5);
					//fprintf(fp,"\\%% g = %d, gap=%d, NumLayers=%d, SizeOfLayer =%d \n", g, gap, layers, sizeOfLayer);

					fprintf(fp2,"g=%d\n gap=%d\n NumLayers=%d\n SizeOfLayer=%d\n parallelism=%d\n NumLightPhaseUnits=%d\n",g,gap,layers,sizeOfLayer, (int)ceil(sizeOfLayer/(1.0*(gap+1))), (int)ceil(g/(1.0*(gap*layers))));
					fprintf(fp2, "HonestCost: %f\n", honestEComp(n,window,R));
					fflush(fp2);
					free(data->S);
					free(data);
				}
				S = selectbSetGapS(gap,layers,n,G);
				int sSize=0;
				int ii;
                for(ii = 0; ii < n; ii++) sSize+=S[ii];


				double cost = attackBH(G, S, sizeOfLayer, gap, n, g, window,  R);
                averageCost+= cost/(1.0*numSamples);
                fprintf(fp2,"Instance %d: |S|=%d, cost=%f \n",k, sSize, cost);

					//printf("Cost: %f \n",cost);
					//printf("Honest Cost: %f \n", honestCost);
					//printf("Quality: %f \n" , honestCost/cost);

				    freeGraph(G,n);
				    free(S);
			}

			double honestCost = honestEComp(n,window,R);
			double attackQuality = honestCost/averageCost;
            fprintf(fp, " (%d,%f) ", nPows[j], attackQuality);
            fflush(fp);
        }
        fprintf(fp,"};\n ");
        fflush(fp);
        fflush(fp2);
    }
}


void generatePlotOriginal(int *nPows, int nLength, int *passes, int passesLength, int numSamples, FILE * fp, FILE * fp2, double R)
{
    int i = 0;
    for(i=0; i < passesLength; i++)
    {
    	int numPasses = passes[i];
    	fprintf(fp, "\\%%  passes = %d, samples =  %d\n", numPasses, numSamples);
    	fprintf(fp,"\\addplot coordinates {");
        int j = 0;
        for( j=0; j < nLength; j++){
        	int n = (int)pow(2,nPows[j]);
        	int window = n/numPasses;

        	int k = 0;
        	double averageCost = 0.0;
        	int * S = NULL;
        	int g, gap, layers, sizeOfLayer;
        	fprintf(fp2, "\n___________________n=%d, passes = %d\n___________________\n",n,numPasses);

        	for( k=0; k < numSamples; k++)
			{

	        	int randSeed = i*numSamples+j*passesLength*numSamples+k;
	        	struct bnode* G = generateArgon2iAGraph(n, window, randSeed);

				if (k==0){
					struct searchData* data = NULL;
					g = (int)pow(n,0.75);
					data = searchbForBestDepthParameters(G, n, window, R, g);

					S = data->S;

					gap = data->gapInLayer;

					layers = data->numLayers;
						//printf("gap %d\n, layers %d\n", gap, layers);

					sizeOfLayer = (int)ceil(n*1.0/(1.0*(layers)));


					g= searchbForg(G, S, n, window, R,  sizeOfLayer, (gap), 5);
					//fprintf(fp,"\\%% g = %d, gap=%d, NumLayers=%d, SizeOfLayer =%d \n", g, gap, layers, sizeOfLayer);

					fprintf(fp2,"g=%d\n gap=%d\n NumLayers=%d\n SizeOfLayer=%d\n parallelism=%d\n NumLightPhaseUnits=%d\n",g,gap,layers,sizeOfLayer, (int)ceil(sizeOfLayer/(1.0*(gap+1))), (int)ceil(g/(1.0*(gap*layers))));
					fprintf(fp2, "HonestCost: %f\n", honestEComp(n,window,R));
					fflush(fp2);
					free(data->S);
					free(data);
				}
				S = selectbSetGapS(gap,layers,n,G);
				int sSize=0;
				int ii;
                for(ii = 0; ii < n; ii++) sSize+=S[ii];


				double cost = attackBH(G, S, sizeOfLayer, gap, n, g, window,  R);
                averageCost+= cost/(1.0*numSamples);
                fprintf(fp2,"Instance %d: |S|=%d, cost=%f \n",k, sSize, cost);

					//printf("Cost: %f \n",cost);
					//printf("Honest Cost: %f \n", honestCost);
					//printf("Quality: %f \n" , honestCost/cost);

				    freeGraph(G,n);
				    free(S);
			}

			double honestCost = honestEComp(n,window,R);
			double attackQuality = honestCost/averageCost;
            fprintf(fp, " (%d,%f) ", nPows[j], attackQuality);
        }
        fprintf(fp,"};\n ");
        fflush(fp);
        fflush(fp2);
    }
}

