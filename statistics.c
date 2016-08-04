#include "statistics.h"
#include "BuildGraph.h"
#include "Attack.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


/* Outputs vector of complexities of attack for increasingly loose upper bounds on parallelism
    For each choice of parallelism do:
     1) Choose initial numLayers (and set S).
     2) Repeat 3 times:
       A) Find optimal g for current numLayers.
       B) Find optimal numLayers (and set S) for current g.
     3) Return final complexity */
void parStats(struct bnode *G, int n, int window, double R, int parUBound)
{
    // Parameters
	const int minParUBound = 50;    // lowest considered parallelism upper bound = 25
	const int numDataPoints = 1;    /* number of upper bounds we will test. each subsequent upper bound is double
                                       the previous => maxParUBound = minParUBound * 2^numDataPoints = 25 * 2^8 = 6400 */
	struct dataPoint {int parVal; double optCost;};
	struct dataPoint *points = (struct dataPoint *) malloc(numDataPoints * sizeof(struct dataPoint));

	// Populate upper bounds (aka DataPoints)
	int i;
	points[0].parVal = minParUBound;
	printf("---------------------------------------------------\n");
	for (i = 1; i < numDataPoints; i++)
	{
		points[i].parVal = minimum((int) (sqrt(n)*1.1), points[i-1].parVal*2);  // double the previous upper bound
	}

	// Find best cost for each data point.
	for (i = 0; i < numDataPoints; i++)
	{
		// 1) choose initial segSize (and matching set S)
		int segSize = (int)(1.1*sqrt(points[i].parVal * 1.0));
		struct sSet *curS = selectPbSetS(n, G, points[i].parVal, segSize);
		int j, optg;
        printf("In parStats beginning parUbound = %d with segSize = %d\n", points[i].parVal, segSize);
		// 2) find optimal g and segSize 3
		for (j=0; j<2; j++) {
			optg = searchbForg(G, curS->S, n, window, R, curS->sizeOfLayer, curS->segSize, 5);
			if (curS) free(curS);
			curS = searchbForSegSize(G, n, window, R, optg, points[i].parVal);
		}
		// 3) Store resulting cost after searching for optimal g and segSize
		points[i].optCost = attackBH(G, curS->S, curS->sizeOfLayer, curS->segSize, n, optg, window, R);
		if (curS) free(curS);
		printf("---------------------------------------------------\n");
	}

	double honestEcomp = n*R + (window + 1.0)*(1.0*window) / 2.0 + (n - window)*1.0*window;
	printf("Parallelism : Quality\n");
	for (i = 0; i < numDataPoints; i++)
	{
		printf("%d : %f\n", points[i].parVal, honestEcomp / points[i].optCost);
	}
}



int searchbForg(struct bnode* G, int* S, int n, int window, double R, int sizeOfLayer, int segSize, int granularity)
{
	// Search Parameters
	const int numDataPoints = 10;

	struct gData { int g; double cost; };
	struct gData* array = malloc(sizeof(struct gData) * numDataPoints);
	int theoryOpt = (int)pow(n, 0.75);
	int d =  (int)ceil(1.0*n / (1.0*sizeOfLayer)) *segSize;
	int gMin = maximum(d, theoryOpt/4);
	int gMax = n;//minimum(theoryOpt * 50,n);
	int gStep = (gMax - gMin) / (numDataPoints-1);
	int g = gMin;
	int i, j;

	int optG = gMin;
	double optCost = INFINITY;

	//int besti = 0;                                                      // DEBUG CODE

	for (i = 0; g <= gMax && i < numDataPoints; i++)
	{
//		printf("g %d, window %d, segSize %d", g, window, segSize);
		array[i].cost = attackBH(G, S, sizeOfLayer, segSize, n, g, window, R);

		array[i].g = g;

		if (optCost > array[i].cost) {
			optCost = array[i].cost;
			optG = g;
			//besti = i;                                                  // DEBUG CODE
			//printf("Current new best i = %d\n", i);                     // DEBUG CODE
		}
		g += gStep;
	}

	//printf("Best i = %d.\n", besti);                                // DEBUG CODE

	for (j = 0; j < granularity; j++)
	{
		g = gMin = maximum(d, optG - gStep);
		gMax = minimum(n, optG + gStep);
		gStep = (gMax - gMin) / (numDataPoints-1);
		//besti = 0;                                                      // DEBUG CODE
		for (i = 0; g <= gMax && i < numDataPoints; i++)
		{
			array[i].cost = attackBH(G, S, sizeOfLayer, segSize, n, g, window, R);

			if (optCost > array[i].cost) {
				optCost = array[i].cost;
				optG = g;
				//besti = i;                                              // DEBUG CODE
				//				printf("Current new best i = %d\n", i);                 // DEBUG CODE
			}
			g += gStep;
		}
		//printf("Best i = %d.\n", besti);                                // DEBUG CODE

	}

	// Noisy output
	printf("Finished SearchbForg with cost %f\n", optCost);           // DEBUG CODE
	double honestEcomp = n*R + (window + 1.0)*(1.0*window) / 2.0 + (n - window)*1.0*window;
	//printf("Honest Ecomp   %f\n", honestEcomp);
	printf("Quality   %f \n", honestEcomp / optCost);
	printf("Opt g %d\n\n", optG);
	free(array);
	return optG;
}


int searchbForgCompressed(struct bnodeCompressed* G, int* S, int n, int window, double R, int sizeOfLayer, int segSize, int granularity)
{
	// Search Parameters
	const int numDataPoints = 10;

	struct gData { int g; double cost; };
	struct gData* array = malloc(sizeof(struct gData) * numDataPoints);
	int theoryOpt = (int)pow(n, 0.75);
	int d =  (int)ceil(1.0*n / (1.0*sizeOfLayer)) *segSize;
	int gMin = maximum(d, theoryOpt/4);
	int gMax = n;//minimum(theoryOpt * 50,n);
	int gStep = (gMax - gMin) / (numDataPoints-1);
	int g = gMin;
	int i, j;

	int optG = gMin;
	double optCost = INFINITY;

	//int besti = 0;                                                      // DEBUG CODE

	for (i = 0; g <= gMax && i < numDataPoints; i++)
	{
//		printf("g %d, window %d, segSize %d", g, window, segSize);
		array[i].cost = attackiXOR(G, S, sizeOfLayer, segSize, n, g, window, R);

		array[i].g = g;

		if (optCost > array[i].cost) {
			optCost = array[i].cost;
			optG = g;
			//besti = i;                                                  // DEBUG CODE
			//printf("Current new best i = %d\n", i);                     // DEBUG CODE
		}
		g += gStep;
	}

	//printf("Best i = %d.\n", besti);                                // DEBUG CODE

	for (j = 0; j < granularity; j++)
	{
		g = gMin = maximum(d, optG - gStep);
		gMax = minimum(n, optG + gStep);
		gStep = (gMax - gMin) / (numDataPoints-1);
		//besti = 0;                                                      // DEBUG CODE
		for (i = 0; g <= gMax && i < numDataPoints; i++)
		{
			array[i].cost = attackiXOR(G, S, sizeOfLayer, segSize, n, g, window, R);

			if (optCost > array[i].cost) {
				optCost = array[i].cost;
				optG = g;
				//besti = i;                                              // DEBUG CODE
				//				printf("Current new best i = %d\n", i);                 // DEBUG CODE
			}
			g += gStep;
		}
		//printf("Best i = %d.\n", besti);                                // DEBUG CODE

	}

	// Noisy output
	printf("Finished SearchbForg with cost %f\n", optCost);           // DEBUG CODE
	double honestEcomp = n*R + (window + 1.0)*(1.0*window) / 2.0 + (n - window)*1.0*window;
	//printf("Honest Ecomp   %f\n", honestEcomp);
	printf("Quality   %f \n", honestEcomp / optCost);
	printf("Opt g %d\n\n", optG);
	free(array);
	return optG;
}


// Find optimal segSize for remaining given parameters.
struct sSet *searchbForSegSize(struct bnode* G, int n, int window, double R, int g, int parUBound)
{
	// Search Parameters
	const int numDataPoints = 9;   // keep this odd.
	const int granularity = 5;

	struct dData {int segSize; int segNum; double cost; };
	struct dData* array = malloc(sizeof(struct dData) * numDataPoints);

    // Search range for segSize
    int segSizeTheoryOpt = maximum(1, (int) ((1.0)*n / (1.1*pow(parUBound,1.5))));    // Optimal #layers is a little bit bigger than sqrt(p) whenever p^2 << n
	int segSizeMin = maximum(1, (int) segSizeTheoryOpt / 10);
	int segSizeMax = minimum(n, segSizeTheoryOpt * 10);
	int segSizeStep = (segSizeMax - segSizeMin) / (numDataPoints-1);
	int segSize = segSizeMin;
	int optSegSize = segSizeMin;

	// Search range for segNum
//	int segNumTheoryOpt = maximum(1, );

	double optCost = INFINITY;

	int i, j;
	struct sSet * curS = 0;

	int besti = 0;                                                      // DEBUG CODE

	for (i = 0; segSize <= segSizeMax && i < numDataPoints; i++)
	{
		if (curS) {
			if (curS->S) { free(curS->S); } free(curS);
		}
		curS = selectPbSetS(n, G, parUBound, segSize);

		array[i].cost = attackBH(G, curS->S, curS->sizeOfLayer, curS->segSize, n, g, window, R);
		array[i].segSize = segSize;
		printf("segSize = %d\t cost = %f\n",array[i].segSize, array[i].cost);  // DEBUG CODE

		if (optCost > array[i].cost) {
			optCost = array[i].cost;
			optSegSize = segSize;
			besti = i;                                                  // DEBUG CODE
			//printf("Current new best i = %d\n", i);                     // DEBUG CODE
		}
		segSize += segSizeStep;
	}

	printf("Best i = %d\n", besti);                                     // DEBUG CODE

	for (j = 0; j < granularity-1; j++)
	{
		segSize = segSizeMin = maximum(1, optSegSize - segSizeStep);
		segSizeMax = minimum(n, optSegSize + segSizeStep);
		segSizeStep = (segSizeMax - segSizeMin) / (numDataPoints-1);
		besti = 0;                                                      // DEBUG CODE
		for (i = 0; segSize <= segSizeMax && i < numDataPoints; i++)
		{
			if (curS) {
				if (curS->S) { free(curS->S); } free(curS);
			}
			curS = selectPbSetS(n, G, parUBound, segSize);

			array[i].cost = attackBH(G, curS->S, curS->sizeOfLayer, curS->segSize, n, g, window, R);
			array[i].segSize = segSize;
			printf("segSize = %d\t cost = %f\n",array[i].segSize, array[i].cost);  // DEBUG CODE

			if (optCost > array[i].cost) {
				optCost = array[i].cost;
				optSegSize = segSize;
				besti = i;                                              // DEBUG CODE
				//printf("Current new best i = %d\n", i);                 // DEBUG CODE
			}
			segSize += segSizeStep;
		}
		printf("Best i = %d.\n", besti);                                // DEBUG CODE

	}

	// Output results and return optimal set S.
	if (curS) {
		if (curS->S) { free(curS->S); } free(curS);
	}
	printf("\nFinished SearchForSegSize with cost %f\n", optCost);
	double honestEcomp = n*R + (window + 1.0)*(1.0*window) / 2.0 + (n - window)*1.0*window;
	//printf("Honest Ecomp   %f\n", honestEcomp);
	printf("Quality   %f \n", honestEcomp / optCost);
	int numLayers =  (int)ceil(n*1.0 / (parUBound*(optSegSize+1)*1.0));
	int layerSize = n/numLayers;
	printf("Opt segSize %d\tnumLayers = %d\tlayerSize = %d\tnumSeg = %d\n", optSegSize, numLayers, layerSize, layerSize/optSegSize);

	return selectPbSetS(n, G, parUBound, optSegSize);
}
struct searchData *searchbForBestDepthParametersCompressed(struct bnodeCompressed* G,   int n, int window, double R, int g)
{
	int numDataPointsX = 5;
	int numDataPointsY = 5;
	int granularity =3;

	int *S=0;
	int i=0;
	int j = 0;
	int theoryOptGap = (int)pow(n, 0.25);
	int theoryOptLayers = (int)pow(n,0.25);
	int GapMin = theoryOptGap/4;
	int GapMax = (int)floor(sqrt(g*1.0));

	int LayerMin = theoryOptLayers/4;
	int LayerMax= (int)floor(sqrt(g*1.0));

	int GapStep =  (GapMax - GapMin) / (numDataPointsX-1);
	int LayerStep = (LayerMax-LayerMin)/(numDataPointsY-1);
	int gap = GapMin, layers = LayerMin;
    printf("GapStep %d\n GapMax %d\n GapMin %d\n",GapStep,GapMax,GapMin);
    printf("LayerStep %d\n LayerMax %d\n LayerMin %d\n",LayerStep,LayerMax,LayerMin);
	double optCost = INFINITY;
	double bestLayers, bestGap;
	for (i = 0; gap <= GapMax && i < numDataPointsX; i++)
	{
        layers = LayerMin;
		for(j=0; layers<=LayerMax && j < numDataPointsY; j++)
		{
			if (S) free(S);
			printf("Building S\n");
			S = selectbSetGapSCompressed(gap,layers, n,window, G);
			printf("Built S\n");
			int sizeOfLayer = (int)ceil(n*1.0 / (layers*1.0));
			printf("Attacking with S\n");
			double curCost = attackiXOR(G, S, sizeOfLayer, gap, n, g, window, R);
            printf("Attacked\n");


			if (optCost > curCost) {
				bestLayers = layers;
				bestGap = gap;
				optCost = curCost;
				//besti = i;                                                  // DEBUG CODE
				//			printf("Current new best i = %d\n", i);                     // DEBUG CODE
			}
			layers += LayerStep;
		}
		gap += GapStep;
	}

    int k;
	for (k = 0; k < granularity; k++)
	{
		GapMin= maximum(1,bestGap-GapStep);
		gap = GapMin;
		LayerMin = maximum(LayerMin,bestLayers-LayerStep);
		//layers = bestLayers - LayerStep;
		GapMax = minimum(GapMax, gap + GapStep);
		LayerMax = minimum(LayerMax, bestLayers+LayerStep);
		GapStep = (GapMax - gap) / (numDataPointsX-1);
		LayerStep = (LayerMax-LayerMin)/(numDataPointsY-1);
	    printf("GapStep %d\n GapMax %d\n GapMin %d\n",GapStep,GapMax,gap);
	    printf("LayerStep %d\n LayerMax %d\n LayerMin %d\n",LayerStep,LayerMax,LayerMin);
		//besti = 0;                                                      // DEBUG CODE
		for (i = 0; gap <= GapMax && i < numDataPointsX; i++)
		{
			layers = LayerMin;
			for(j=0; layers<=LayerMax && j < numDataPointsY; j++)
			{
				if (S) free(S);
				S = selectbSetGapSCompressed(gap,layers, n, window, G);

				int sizeOfLayer = (int)ceil(n*1.0 / (layers*1.0));
				double curCost = attackiXOR(G, S, sizeOfLayer, gap, n, g, window, R);



				if (optCost > curCost) {
					bestLayers = layers;
					bestGap = gap;
					optCost = curCost;
					//besti = i;                                                  // DEBUG CODE
					//			printf("Current new best i = %d\n", i);                     // DEBUG CODE
				}
				layers += LayerStep;
			}
			gap += GapStep;
		}		//printf("Best i = %d.\n", besti);                                // DEBUG CODE


	}
	//printf("debug here 3\n");
	struct searchData* data = malloc(sizeof(struct searchData));

	if (S) free(S);
	data->g = g;
	data->gapInLayer = bestGap;
	data->numLayers = bestLayers;
	data->S = selectbSetGapSCompressed(bestGap,bestLayers, n, window, G);
	//*Sout =
	//printf("debug here 3\n");
	//*gapOut = bestGap;
	//printf("debug here 4\n");

	//*layersOut = bestLayers;

	return data;

}
struct searchData *searchbForBestDepthParameters(struct bnode* G,  int n, int window, double R, int g)
{
	int numDataPointsX = 5;
	int numDataPointsY = 5;
	int granularity =3;

	int *S=0;
	int i=0;
	int j = 0;
	int theoryOptGap = (int)pow(n, 0.25);
	int theoryOptLayers = (int)pow(n,0.25);
	int GapMin = theoryOptGap/4;
	int GapMax = (int)floor(sqrt(g*1.0));

	int LayerMin = theoryOptLayers/4;
	int LayerMax= (int)floor(sqrt(g*1.0));

	int GapStep =  (GapMax - GapMin) / (numDataPointsX-1);
	int LayerStep = (LayerMax-LayerMin)/(numDataPointsY-1);
	int gap = GapMin, layers = LayerMin;
    printf("GapStep %d\n GapMax %d\n GapMin %d\n",GapStep,GapMax,GapMin);
    printf("LayerStep %d\n LayerMax %d\n LayerMin %d\n",LayerStep,LayerMax,LayerMin);
	double optCost = INFINITY;
	double bestLayers, bestGap;
	for (i = 0; gap <= GapMax && i < numDataPointsX; i++)
	{
        layers = LayerMin;
		for(j=0; layers<=LayerMax && j < numDataPointsY; j++)
		{
			if (S) free(S);
			S = selectbSetGapS(gap,layers, n, G);
			int sizeOfLayer = (int)ceil(n*1.0 / (layers*1.0));
			double curCost = attackBH(G, S, sizeOfLayer, gap, n, g, window, R);



			if (optCost > curCost) {
				bestLayers = layers;
				bestGap = gap;
				optCost = curCost;
				//besti = i;                                                  // DEBUG CODE
				//			printf("Current new best i = %d\n", i);                     // DEBUG CODE
			}
			layers += LayerStep;
		}
		gap += GapStep;
	}

    int k;
	for (k = 0; k < granularity; k++)
	{
		GapMin= maximum(1,bestGap-GapStep);
		gap = GapMin;
		LayerMin = maximum(LayerMin,bestLayers-LayerStep);
		//layers = bestLayers - LayerStep;
		GapMax = minimum(GapMax, gap + GapStep);
		LayerMax = minimum(LayerMax, bestLayers+LayerStep);
		GapStep = (GapMax - gap) / (numDataPointsX-1);
		LayerStep = (LayerMax-LayerMin)/(numDataPointsY-1);
	    printf("GapStep %d\n GapMax %d\n GapMin %d\n",GapStep,GapMax,gap);
	    printf("LayerStep %d\n LayerMax %d\n LayerMin %d\n",LayerStep,LayerMax,LayerMin);
		//besti = 0;                                                      // DEBUG CODE
		for (i = 0; gap <= GapMax && i < numDataPointsX; i++)
		{
			layers = LayerMin;
			for(j=0; layers<=LayerMax && j < numDataPointsY; j++)
			{
				if (S) free(S);
				S = selectbSetGapS(gap,layers, n, G);

				int sizeOfLayer = (int)ceil(n*1.0 / (layers*1.0));
				double curCost = attackBH(G, S, sizeOfLayer, gap, n, g, window, R);



				if (optCost > curCost) {
					bestLayers = layers;
					bestGap = gap;
					optCost = curCost;
					//besti = i;                                                  // DEBUG CODE
					//			printf("Current new best i = %d\n", i);                     // DEBUG CODE
				}
				layers += LayerStep;
			}
			gap += GapStep;
		}		//printf("Best i = %d.\n", besti);                                // DEBUG CODE


	}
	//printf("debug here 3\n");
	struct searchData* data = malloc(sizeof(struct searchData));

	if (S) free(S);
	data->g = g;
	data->gapInLayer = bestGap;
	data->numLayers = bestLayers;
	data->S = selectbSetGapS(bestGap,bestLayers, n, G);
	//*Sout =
	//printf("debug here 3\n");
	//*gapOut = bestGap;
	//printf("debug here 4\n");

	//*layersOut = bestLayers;

	return data;

}
struct searchData *searchbForBestDepthParametersBddParallelism(struct bnode* G,  int n, int window, double R, int g,int parBound)
{
	int numDataPointsX = 5;
	int numDataPointsY = 5;
	int granularity =3;

	int *S=0;
	int i=0;
	int j = 0;
	int theoryOptGap = maximum(1, (int) ((1.0)*n / (1.1*pow(parBound,1.5)))); // Optimal #layers is a little bit bigger than sqrt(p) whenever p^2 << n
	int theoryOptLayers = (int)pow(n,0.25);
	int GapMin = theoryOptGap/4;
	int GapMax = (int)floor(2*sqrt(g*1.0));



	int LayerMin = theoryOptLayers/4;
	int LayerMax= (int)ceil(g*1.0/(1.0*GapMin));

	int GapStep =  (GapMax - GapMin) / (numDataPointsX-1);
	int LayerStep = (LayerMax-LayerMin)/(numDataPointsY-1);
	int gap = GapMin, layers = LayerMin;
    printf("GapStep %d\n GapMax %d\n GapMin %d\n",GapStep,GapMax,GapMin);
   double optCost = INFINITY;
	double bestLayers, bestGap;
	for (i = 0; gap <= GapMax && i < numDataPointsX; i++)
	{
		int maxLayerSizeAtGap = parBound*gap;
		int layerLowerBound = (int)ceil(n*1.0/(1.0*maxLayerSizeAtGap));

		LayerMin = maximum(LayerMin,layerLowerBound);
		LayerStep = (LayerMax-LayerMin)/(numDataPointsY-1);
		layers = LayerMin;
		printf("LayerStep %d\n LayerMax %d\n LayerMin %d\n",LayerStep,LayerMax,LayerMin);
		for(j=0; layers<=LayerMax && j < numDataPointsY; j++)
		{
			if (S) free(S);
			S = selectbSetGapS(gap,layers, n, G);
			int sizeOfLayer = (int)ceil(n*1.0 / (layers*1.0));
			double curCost = attackBH(G, S, sizeOfLayer, gap, n, g, window, R);



			if (optCost > curCost) {
				bestLayers = layers;
				bestGap = gap;
				optCost = curCost;
				//besti = i;                                                  // DEBUG CODE
				//			printf("Current new best i = %d\n", i);                     // DEBUG CODE
			}
			layers += LayerStep;
		}
		gap += GapStep;
	}

    int k;
	for (k = 0; k < granularity; k++)
	{
		gap = bestGap - GapStep;
		LayerMin = maximum(LayerMin,bestLayers-LayerStep);
		//layers = bestLayers - LayerStep;
		GapMax = minimum(GapMax, gap + GapStep);
		LayerMax = minimum(LayerMax, bestLayers+LayerStep);
		GapStep = (GapMax - gap) / (numDataPointsX-1);
		LayerStep = (LayerMax-LayerMin)/(numDataPointsY-1);
	    printf("GapStep %d\n GapMax %d\n GapMin %d\n",GapStep,GapMax,gap);
	    printf("LayerStep %d\n LayerMax %d\n LayerMin %d\n",LayerStep,LayerMax,LayerMin);
		//besti = 0;                                                      // DEBUG CODE
		for (i = 0; gap <= GapMax && i < numDataPointsX; i++)
		{
			int maxLayerSizeAtGap = parBound*gap;
			int layerLowerBound = (int)ceil(n*1.0/(1.0*maxLayerSizeAtGap));

			LayerMin = maximum(LayerMin,layerLowerBound);
			LayerStep = (LayerMax-LayerMin)/(numDataPointsY-1);

			layers = LayerMin;
			for(j=0; layers<=LayerMax && j < numDataPointsY; j++)
			{
				if (S) free(S);
				S = selectbSetGapS(gap,layers, n, G);

				int sizeOfLayer = (int)ceil(n*1.0 / (layers*1.0));
				double curCost = attackBH(G, S, sizeOfLayer, gap, n, g, window, R);



				if (optCost > curCost) {
					bestLayers = layers;
					bestGap = gap;
					optCost = curCost;
					//besti = i;                                                  // DEBUG CODE
					//			printf("Current new best i = %d\n", i);                     // DEBUG CODE
				}
				layers += LayerStep;
			}
			gap += GapStep;
		}		//printf("Best i = %d.\n", besti);                                // DEBUG CODE


	}
	//printf("debug here 3\n");
	struct searchData* data = malloc(sizeof(struct searchData));

	if (S) free(S);
	data->g = g;
	data->gapInLayer = bestGap;
	data->numLayers = bestLayers;
	data->S = selectbSetGapS(bestGap,bestLayers, n, G);
	//*Sout =
	//printf("debug here 3\n");
	//*gapOut = bestGap;
	//printf("debug here 4\n");

	//*layersOut = bestLayers;

	return data;

}
// Return rootd, set *Sout = S
int searchbFord(struct bnode* G, int ** Sout, int n, int window, double R, int g)
{
	// Search Parameters
	int numDataPoints = 9;
	int granularity = 5;

	struct dData { int rootd; double cost; int * S; };
	struct dData* array = malloc(sizeof(struct dData) * numDataPoints);
	int* S = 0;
	int theoryOpt = (int)pow(n, 0.25);
	int rootdMin = theoryOpt / 4;
	int rootdMax = (int)floor(sqrt(g*1.0));//minimum(theoryOpt * 50,n);
	int rootdStep = (rootdMax - rootdMin) / (numDataPoints-1);
	int rootd = rootdMin;

	int i, j;

	int optrootD = rootdMin;
	double optCost = INFINITY;

	//int besti = 0;                                                      // DEBUG CODE

	for (i = 0; rootd <= rootdMax && i < numDataPoints; i++)
	{
		if (S) free(S);
		S = selectbSetS(rootd, n, G);

		int sizeOfLayer = (int)ceil(n*1.0 / (rootd*1.0));
		int gapInLayer = rootd;
		array[i].cost = attackBH(G, S, sizeOfLayer, gapInLayer, n, g, window, R);
		array[i].rootd = rootd;
		array[i].S = 0;

		if (optCost > array[i].cost) {
			optCost = array[i].cost;
			optrootD = rootd;
			//besti = i;                                                  // DEBUG CODE
			//			printf("Current new best i = %d\n", i);                     // DEBUG CODE
		}
		rootd += rootdStep;
	}


	for (j = 0; j < granularity; j++)
	{
		rootd = optrootD - rootdStep;
		rootdMax = minimum(rootdMax, optrootD + rootdStep);
		rootdStep = (rootdMax - rootd) / (numDataPoints-1);
		//besti = 0;                                                      // DEBUG CODE
		for (i = 0; rootd <= rootdMax && i < numDataPoints; i++)
		{
			int sizeOfLayer = (int)ceil(n*1.0 / (rootd*1.0));
			int segSize = rootd;
			if (S) free(S);
			S = selectbSetS(rootd, n, G);
			array[i].cost = attackBH(G, S, sizeOfLayer, segSize, n, g, window, R);

			if (optCost > array[i].cost) {
				optCost = array[i].cost;
				optrootD = rootd;
				//besti = i;                                              // DEBUG CODE
				//				printf("Current new best i = %d\n", i);                 // DEBUG CODE
			}
			rootd += rootdStep;
		}
		//printf("Best i = %d.\n", besti);                                // DEBUG CODE

	}
	if (S) free(S);
	*Sout = selectbSetS(optrootD, n, G);
	printf("\nSearchbFord done with cost %f\n", optCost);
	double honestEcomp = n*R + (window + 1.0)*(1.0*window) / 2.0 + (n - window)*1.0*window;
	//printf("Honest Ecomp   %f\n", honestEcomp);
	printf("Quality   %f \n", honestEcomp / optCost);
	int numLayers =  optrootD;
	int layerSize = (int)ceil( n*1.0/(1.0*optrootD));
	printf("Opt rootd %d\tnumLayers = %d\tlayerSize = %d\tsegNum =%d\n", optrootD, numLayers, layerSize, layerSize/(optrootD+1));
	//printf("n =   %d\n", n);
	return optrootD;

}


// ================
// Helper Functions
//=================
//
// Returns maximum of two integer arguments
int maximum(int a, int b) { return (a > b) ? a : b; }
