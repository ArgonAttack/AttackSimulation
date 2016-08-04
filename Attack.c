// ConsoleApplication7.cpp : Defines the entry point for the console application.
//
#include "Attack.h"
#include "BuildGraph.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// potentially can be optimized further
double balloonPhaseCostParallel(struct bnode* G, int* S,  int sizeOfLayer, int gapInLayer, int nOverP, int w, double R, int startTime, int parallelism)
{
	//int sizeOfLayer = (int)ceil(n*1.0 / (rootd*1.0));
	//int gapInLayer = rootd;
	int currentLayer = (int)floor(startTime*1.0 / (1.0*sizeOfLayer));
	int balloonPhaseTime = (currentLayer + 1)*gapInLayer;              //total time balloon phase runs under heuristic
	int endTime = startTime + balloonPhaseTime;                                //last step of balloon phase  under heuristic

	double cost = 0.0;
	int j;
	int iL, lane;
	//calculates cost up until round i
	for (iL = 0; iL < startTime; iL++)
	{
		for(lane = 0; lane < parallelism; lane++)
		{
			j = iL*parallelism+lane;
			if (S[j] == 0){           // already accounted for cost of pebbles in S or in parents
				//Get the layer of node j, the index within layer and the place of node j within its segment on that layer
				int layerOfNode = (int)floor(iL*1.0 / (1.0*sizeOfLayer));
				int layerOfNodePlusW = (int)floor((iL + w/parallelism)*1.0 / (1.0*sizeOfLayer));
				int endOfWindow = (layerOfNodePlusW + 1)*gapInLayer;    //(upper bound on) time until node j is outside of memory window
				//int placeInLayer = iL % sizeOfLayer;
				int placeInSegment = iL % (gapInLayer + 1);
				int timeJFirstPebbled = layerOfNode*gapInLayer + placeInSegment;

				// either keep pebble till end of balloon phase or end of window
				int timeOnGraph = minimum(balloonPhaseTime - timeJFirstPebbled + 1, endOfWindow - timeJFirstPebbled + 1);

				//cost to place pebble plus cost to store pebble
				cost += 1.0*R + timeOnGraph;
			}
		}
	}

	int lastBalloonPhaseRound = minimum(endTime, nOverP - 1);

	// calculates cost on nodes [i,lastBalloonPhaseNode]
	for (iL = startTime; iL <= lastBalloonPhaseRound; iL++)
	{
		for(lane=0; lane < parallelism; lane++) {
			j = iL*parallelism+lane;
			if (S[j] == 0){ // already accounted for cost of pebbles in S or in parents
				cost += 1.0*R + endTime-iL+1;
			}
		}
	}

	return cost;
}


double balloonPhaseCostiXOR(struct bnodeCompressed* G, int* S, int sizeOfLayer, int gapInLayer, int n, int w, double R, int i)
{
	int currentLayer = nodeLayer(i, sizeOfLayer); // (int)floor(i*1.0 / (1.0*sizeOfLayer));
	int balloonPhaseTime = (currentLayer + 1)*gapInLayer;              //total time balloon phase runs under heuristic
	int endTime = i + balloonPhaseTime;                         //last step of balloon phase  under heuristic

	double cost = 0.0;
	int j;
	//calculates cost up until node i
	for (j = 0; j < i; j++)
	{
		if (S[j] == 0 ){           // already accounted for cost of pebbles in S or in parents
			//Get the layer of node j, the index within layer and the place of node j within its segment on that layer
			int layerOfNodeJ = nodeLayer(j, sizeOfLayer); // (int)floor(j*1.0 / (1.0*sizeOfLayer));
			int layerOfNodeJPlusW = nodeLayer(j + w, sizeOfLayer);//  (int)floor((j + w)*1.0 / (1.0*sizeOfLayer));
			int endOfWindow = (layerOfNodeJPlusW + 1)*gapInLayer;    //(upper bound on) time until node j is outside of memory window
			int placeInLayer = nodePlaceInLayer(j, sizeOfLayer); //				j % sizeOfLayer;
			int placeInSegment = nodePlaceInGap(j, sizeOfLayer, gapInLayer); // placeInLayer % (gapInLayer + 1);
			int timeJFirstPebbled = layerOfNodeJ*gapInLayer + placeInSegment;

			// either keep pebble till end of balloon phase or end of window
			int timeOnGraph = minimum(balloonPhaseTime - timeJFirstPebbled + 1, endOfWindow - timeJFirstPebbled + 1);

			//cost to place pebble plus cost to store pebble
			cost += 1.0*R + timeOnGraph;
		}
	}

	int lastBalloonPhaseNode = minimum(endTime, n - 1);

	// calculates cost on nodes [i,lastBalloonPhaseNode]
	for (j = i; j <= lastBalloonPhaseNode; j++)
	{
		if (S[j] == 0  ){ // already accounted for cost of pebbles in S or in parents
			cost += 1.0*R + endTime-j+1;
		}
	}

	return cost;
}

// potentially can be optimized further
double balloonbPhaseCost(struct bnode* G, int* S, int sizeOfLayer, int gapInLayer, int n, int w, double R, int i)
{
	//int sizeOfLayer = (int) ceil(n*1.0 / (rootd*1.0));
	//int gapInLayer = rootd;
	int currentLayer = nodeLayer(i, sizeOfLayer); // (int)floor(i*1.0 / (1.0*sizeOfLayer));
	int balloonPhaseTime = (currentLayer + 1)*gapInLayer;              //total time balloon phase runs under heuristic
	int endTime = i + balloonPhaseTime;                         //last step of balloon phase  under heuristic

	double cost = 0.0;
	int j;
	//calculates cost up until node i
	for (j = 0; j < i; j++)
	{
		if (S[j] == 0 && IsbParent(G,j,i,endTime)==0){           // already accounted for cost of pebbles in S or in parents
			//Get the layer of node j, the index within layer and the place of node j within its segment on that layer
			int layerOfNodeJ = nodeLayer(j, sizeOfLayer); // (int)floor(j*1.0 / (1.0*sizeOfLayer));
			int layerOfNodeJPlusW = nodeLayer(j + w, sizeOfLayer);//  (int)floor((j + w)*1.0 / (1.0*sizeOfLayer));
			int endOfWindow = (layerOfNodeJPlusW + 1)*gapInLayer;    //(upper bound on) time until node j is outside of memory window
			int placeInLayer = nodePlaceInLayer(j, sizeOfLayer); //				j % sizeOfLayer;
			int placeInSegment = nodePlaceInGap(j, sizeOfLayer, gapInLayer); // placeInLayer % (gapInLayer + 1);
			int timeJFirstPebbled = layerOfNodeJ*gapInLayer + placeInSegment;

			// either keep pebble till end of balloon phase or end of window
			int timeOnGraph = minimum(balloonPhaseTime - timeJFirstPebbled + 1, endOfWindow - timeJFirstPebbled + 1);

			//cost to place pebble plus cost to store pebble
			cost += 1.0*R + timeOnGraph;
		}
	}

	int lastBalloonPhaseNode = minimum(endTime, n - 1);

	// calculates cost on nodes [i,lastBalloonPhaseNode]
	for (j = i; j <= lastBalloonPhaseNode; j++)
	{
		if (S[j] == 0 && IsbParent(G, j, i, endTime) == 0 ){ // already accounted for cost of pebbles in S or in parents
			cost += 1.0*R + endTime-j+1;
		}
	}

	return cost;
}


//Output: the energy cost of attacking G
// Parameter R: core memory ratio
// Parameter rootd: rootd^2 = d >= depth(G-S). We have rootd layers, and each segment in a layer has size rootd.
double attackBH(struct bnode* G, int* S, int sizeOfLayer, int segSize, int n, int g, int w, double R)
{
	if (n <= 0) return -1000; /* Should only call with n > 0*/
	int i;
    resetGraph(G,n);
	//int sizeOfLayer = (int)ceil(n*1.0 / (rootd*1.0));

	/*total memory used keeping pebbles on S in pebbling*/
	double costOnS = 0;
	double costOnP = 0;
	/* number RO queries */

	int lightPhaseEnd = g;
	/*current number of pebbles on S at step i*/
	int pebblesOnSRightNow = 0;
	int numberPebblesonParentsRightNow = 0;

	double costOfBalloonPhases = 0.0;
	int totalBalloonPhaseRound = 0;

	for (i = 1; i < n; i++)
	{
		pebblesOnSRightNow = pebblesOnSRightNow + S[i];
		costOnS = costOnS + pebblesOnSRightNow;

		costOnP = costOnP + numberPebblesonParentsRightNow;

		if (i % g == 0 && i > 0 )
		{
			numberPebblesonParentsRightNow = bparents(G, S, i, i + g, n);
			lightPhaseEnd += g;
		}
		else // if (i >= g)
		{
			numberPebblesonParentsRightNow = updatebParents(G, S, numberPebblesonParentsRightNow, i, lightPhaseEnd);

		}
	}

	for (i =g; i < n; i+=g)
	{
		int currentLayer = nodeLayer(i, sizeOfLayer);  // (int)floor(i * 1.0 / (1.0 * sizeOfLayer));
		int timeBalloonPhaseWouldTake = currentLayer*segSize + nodePlaceInGap(i,sizeOfLayer,segSize);
		costOfBalloonPhases +=  balloonbPhaseCost(G, S, sizeOfLayer, segSize, n, w, R, i - timeBalloonPhaseWouldTake);
		totalBalloonPhaseRound += timeBalloonPhaseWouldTake;
	}

	double cost = costOnS + costOnP + (n - totalBalloonPhaseRound)*R + costOfBalloonPhases;

	// DEBUG CODE
	printf("g = %d\t cost = %0.0f\t costOnP = %0.1f\t costOfBalloonPhases = %0.0f\t totalBalloonPhaseRound = %d\n",  \
	         g, cost, costOnP, costOfBalloonPhases, totalBalloonPhaseRound);

	resetGraph(G, n);
	return cost;
}


//Output: the energy cost of attacking G
// Parameter R: core memory ratio
// Parameter rootd: rootd^2 = d >= depth(G-S). We have n/sizeOfLayer layers, and each segment in a layer has size segSize.
double attackiXOR(struct bnodeCompressed* G, int* S, int sizeOfLayer, int segSize, int n, int g, int w, double R)
{
	if (n <= 0) return -1000; /* Should only call with n > 0*/
	int i;

	//int sizeOfLayer = (int)ceil(n*1.0 / (rootd*1.0));

	/*total memory used keeping pebbles on S in pebbling*/
	double costOnS = 0;
	double costOnP = 0;
	/* number RO queries */

	int lightPhaseEnd = g;
	/*current number of pebbles on S at step i*/
	int pebblesOnSRightNow = 0;
	int numberPebblesonParentsRightNow = 1;

	double costOfBalloonPhases = 0.0;
	int totalBalloonPhaseRound = 0;

	for (i = 1; i < n; i++)
	{
		pebblesOnSRightNow = pebblesOnSRightNow + S[i];
		costOnS = costOnS + pebblesOnSRightNow;
		costOnP = costOnP + numberPebblesonParentsRightNow;
        if (i < g)
        {
        	numberPebblesonParentsRightNow = maximum(i,g-i);
        }
        else if (i % g == 0 && i > 0 )
		{
			numberPebblesonParentsRightNow = g;
			lightPhaseEnd += g;
		}
		else // if (i >= g)
		{
			numberPebblesonParentsRightNow = lightPhaseEnd-i;

		}

	}

	for (i =g; i < n; i+=g)
	{
		int currentLayer = nodeLayer(i, sizeOfLayer);  // (int)floor(i * 1.0 / (1.0 * sizeOfLayer));
		int timeBalloonPhaseWouldTake = currentLayer*segSize + nodePlaceInGap(i,sizeOfLayer,segSize);
		costOfBalloonPhases +=  balloonPhaseCostiXOR(G, S, sizeOfLayer, segSize, n, w, R, i - timeBalloonPhaseWouldTake);
		totalBalloonPhaseRound += timeBalloonPhaseWouldTake;
	}

	double cost = costOnS + costOnP + (n - totalBalloonPhaseRound)*R + costOfBalloonPhases;

	// DEBUG CODE
	printf("g = %d\t cost = %0.0f\t costOnP = %0.1f\t costOfBalloonPhases = %0.0f\t totalBalloonPhaseRound = %d\n",  \
	         g, cost, costOnP, costOfBalloonPhases, totalBalloonPhaseRound);


	return cost;
}

double attackParallelLanes(struct bnode* G, int* S, int sizeOfLayer, int gapInLayer, int nOverP, int g, int w, double R, int parallelism)
{
	if (nOverP <= 0 || parallelism <= 0) return -1000; /* Should only call with n > 0*/
	int i,iL, lane;
    resetGraph(G,nOverP*parallelism);
	//int sizeOfLayer = (int)ceil(n*1.0 / (rootd*1.0));

	/*total memory used keeping pebbles on S in pebbling*/
	double costOnS = 0;
	double costOnP = 0;
	/* number RO queries */

	int lightPhaseEnd = g;
	/*current number of pebbles on S at step i*/
	int pebblesOnSRightNow = 0;
	int numberPebblesonParentsRightNow = 0;

	double costOfBalloonPhases = 0.0;
	int totalBalloonPhaseRound = 0;

	for (iL = 1; iL < nOverP; iL++)
	{
		if (iL % g == 0 && iL > 0 )
		{
			numberPebblesonParentsRightNow = bparents(G, S, iL*parallelism, iL*parallelism + g*parallelism, nOverP*parallelism);
			lightPhaseEnd += g;
		}
		else // if (i >= g)
		{
			for(lane = 0; lane < parallelism; lane++)
			{
			numberPebblesonParentsRightNow = updatebParents(G, S, numberPebblesonParentsRightNow, iL*parallelism+lane, lightPhaseEnd*parallelism);
			}

		}
        for(lane = 0; lane < parallelism; lane++)
        {
			i = iL*parallelism+lane;
		    pebblesOnSRightNow = pebblesOnSRightNow + S[i];
		}
	    costOnS = costOnS + pebblesOnSRightNow;

	    costOnP = costOnP + numberPebblesonParentsRightNow;
	}

	for (iL =g; iL < nOverP; iL+=g)
	{
		int currentLayer = nodeLayer(iL, sizeOfLayer);  // (int)floor(i * 1.0 / (1.0 * sizeOfLayer));
		int timeBalloonPhaseWouldTake = currentLayer*gapInLayer + nodePlaceInGap(iL,sizeOfLayer,gapInLayer);

		costOfBalloonPhases +=  balloonPhaseCostParallel(G,  S, sizeOfLayer, gapInLayer, nOverP, w, R, iL - timeBalloonPhaseWouldTake,parallelism);
		totalBalloonPhaseRound += timeBalloonPhaseWouldTake;
	}

	double cost = costOnS + costOnP + (nOverP - totalBalloonPhaseRound)*parallelism*R + costOfBalloonPhases;

	// DEBUG CODE
	printf("g = %d\t cost = %0.0f\t costOnP = %0.1f\t costOfBalloonPhases = %0.0f\t totalBalloonPhaseRound = %d\n",  \
	         g, cost, costOnP, costOfBalloonPhases, totalBalloonPhaseRound);

	resetGraph(G, nOverP*parallelism);
	return cost;
}


