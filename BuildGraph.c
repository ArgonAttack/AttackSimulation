#include "BuildGraph.h"
#include "Attack.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


// ================
// Helper Functions
//=================
//
// Returns minimum of two integer arguments
int minimum(int a, int b) { return (a < b) ? a : b; }

// returns new integer array consisting of src with newItem appended and and frees src
int * intArrayAppend(int * src, int len, int newItem)
{
	int * p = (int *)malloc((len + 1) * sizeof(int));
	memcpy(p, src, len * sizeof(int));
	p[len] = newItem;
	if (src) free(src);
	return p;
}


// ========================
// Graph Sampling Functions
// ========================
//
// Sample an Argon2i Graph as a bnode
struct bnode* generateArgon2iAGraph(int n, int w, int seed)
{
	srand(seed);
	struct bnode* graph = malloc(n*sizeof(struct bnode));
	int i, modulus;
	graph[0].childCount = 0;
	graph[0].children = 0;
	graph[0].parentCount = 0;
	graph[0].parents = 0;

	for (i = 1; i < n; i++)
	{
		graph[i].childCount = 0;
		graph[i].peb = 0;
		graph[i].children = 0;
		graph[0].parentCount = 0;
		graph[0].parents = 0;
		modulus = (i < w) ? i : w;
		int p = i-(randomNumber() % modulus);


		// Add p to parent list of node i.
		graph[i].parents = intArrayAppend(graph[i].parents, graph[i].parentCount, p);
		graph[i].parentCount += 1;
		// Add node i to children list of p.
		graph[p].children = intArrayAppend(graph[p].children, graph[p].childCount, i);
		graph[p].childCount += 1;
		/*graph[p].childCount++;
		int count = graph[p].childCount;
		int* temp = graph[p].children;
		graph[p].children = malloc(count*sizeof(int));
		int j;
		for (j = 0; j < count - 1; j++)
		{
			graph[p].children[j] = temp[j];
		}
		graph[p].children[count - 1] = i;
		if (temp) free(temp); */

	}
	return graph;
}

int randomNumber()
{
    int r1 = rand();
    int r2 = rand();
    return (r1 << 15) + r2;
}
// nOverP must be divisible by parallelism and w
// w must be a divisible 4*parallelism
struct bnode* generateArgon2iBGraphParallelism(int nOverP, int w, int seed, int parallelism)
{
	srand(seed);
	struct bnode* graph = malloc(nOverP*parallelism*sizeof(struct bnode));
	int i, j;
	int lane;
	int indeg = 3;
	graph[0].childCount = 0;
	graph[0].children = 0;
	graph[0].parentCount = 0;
	graph[0].parents = 0;
    int sliceSizeInLane = w/(4*parallelism);
	int iL;
	for (iL = 1; iL < nOverP; iL++)
	{
		int node = iL*parallelism;
		int slice = (int)floor(4.0*w/(1.0*node));

		for (lane = 0 ; lane < parallelism; lane++)
		{

			i = node+lane;
			graph[i].childCount = 0;
			graph[i].peb = 0;
			graph[i].children = 0;
			graph[i].parentCount = 0;
			graph[i].parents = 0;
			int numParentsExcludingXOR = minimum(i, indeg-2);	//Number of parents for node i excluding XOR node (i-w) and node (i-1).
			// Populate parent list and add i to children list of its parents.
			for (j = 0; j < numParentsExcludingXOR; j++) {
				int p;
				// Sample parent
				int rL = randomNumber() % parallelism;
				if(iL < sliceSizeInLane || rL == lane)
				{
					int modulus = (iL*parallelism > w? w/parallelism: iL);//
					int r = randomNumber();
					int x = floor(r*1.0*r*1.0/pow(2.0,30.0));
					int y = floor((modulus*1.0*x*1.0)/pow(2.0,30.0));
					int z = modulus-1-y;

					p = (iL-modulus + z)*parallelism+rL;
				}
				else
				{

					int iLinSlice = iL % sliceSizeInLane;
					int sliceStartInLane = iL - iLinSlice;
					int modulus = (sliceStartInLane*parallelism > w? w/parallelism: sliceStartInLane);
					int r = randomNumber();
					int x = floor(r*1.0*r*1.0/pow(2.0,30.0));
					int y = floor((modulus*1.0*x*1.0)/pow(2.0,30.0));
					int z = modulus-1-y;

					p = (sliceStartInLane-modulus + z)*parallelism+rL;
				}
				//int p = i-(randomNumber() % modulus);
				// Add p to parent list of node i.
				graph[i].parents = intArrayAppend(graph[i].parents, graph[i].parentCount, p);
				graph[i].parentCount += 1;
				// Add node i to children list of p.
				graph[p].children = intArrayAppend(graph[p].children, graph[p].childCount, i);
				graph[p].childCount += 1;
			}
			// Add edge (i-w,i) if the current node i>=w
			if (i >= w)
			{
				int p = i - w;
				// Add p to parent list of node i.
				graph[i].parents = intArrayAppend(graph[i].parents, graph[i].parentCount, p);
				graph[i].parentCount += 1;
				// Add node i to children list of p.
				graph[p].children = intArrayAppend(graph[p].children, graph[p].childCount, i);
				graph[p].childCount += 1;
			}


		}
    }
	return graph;

}
struct bnode* generateArgon2iBGraph(int n, int w, int seed)
{
	srand(seed);
	struct bnode* graph = malloc(n*sizeof(struct bnode));
	int i, j;
	int indeg = 3;
	graph[0].childCount = 0;
	graph[0].children = 0;
	graph[0].parentCount = 0;
	graph[0].parents = 0;

	for (i = 1; i < n; i++)
	{
		graph[i].childCount = 0;
		graph[i].peb = 0;
		graph[i].children = 0;
		graph[i].parentCount = 0;
		graph[i].parents = 0;
		int numParentsExcludingXOR = minimum(i, indeg-2);	//Number of parents for node i excluding XOR node (i-w) and node (i-1).
		// Populate parent list and add i to children list of its parents.
		for (j = 0; j < numParentsExcludingXOR; j++) {
			// Sample parent
			int modulus = (i < w) ? i : w;
			int r = randomNumber();
			int x = floor(r*1.0*r*1.0/pow(2.0,30.0));
			int y = floor((modulus*1.0*x*1.0)/pow(2.0,30.0));
			int z = modulus-1-y;
			int p = i-modulus + z;
			//int p = i-(randomNumber() % modulus);
			// Add p to parent list of node i.
			graph[i].parents = intArrayAppend(graph[i].parents, graph[i].parentCount, p);
			graph[i].parentCount += 1;
			// Add node i to children list of p.
			graph[p].children = intArrayAppend(graph[p].children, graph[p].childCount, i);
			graph[p].childCount += 1;
		}
		// Add edge (i-w,i) if the current node i>=w
		if (i >= w)
		{
			int p = i - w;
			// Add p to parent list of node i.
			graph[i].parents = intArrayAppend(graph[i].parents, graph[i].parentCount, p);
			graph[i].parentCount += 1;
			// Add node i to children list of p.
			graph[p].children = intArrayAppend(graph[p].children, graph[p].childCount, i);
			graph[p].childCount += 1;
		}


	}
	return graph;

}

// Sample Argon3i (e.g., latest `tweak') Hashing Linear Graph
struct bnode* generateBHGraph(int n, int w, int seed)
{
	srand(seed);
	struct bnode* graph = malloc(n*sizeof(struct bnode));
	int i, j;
	int indeg = 3;
	graph[0].childCount = 0;
	graph[0].children = 0;
	graph[0].parentCount = 0;
	graph[0].parents = 0;

	for (i = 1; i < n; i++)
	{
		graph[i].childCount = 0;
		graph[i].peb = 0;
		graph[i].children = 0;
		graph[i].parentCount = 0;
		graph[i].parents = 0;
		int numParentsExcludingXOR = minimum(i, indeg-2);	//Number of parents for node i excluding XOR node (i-w) and node (i-1).
		// Populate parent list and add i to children list of its parents.
		for (j = 0; j < numParentsExcludingXOR; j++) {
			// Sample parent
			int modulus = (i < w) ? i : w;
			int p = i-(randomNumber() % modulus);
			// Add p to parent list of node i.
			graph[i].parents = intArrayAppend(graph[i].parents, graph[i].parentCount, p);
			graph[i].parentCount += 1;
			// Add node i to children list of p.
			graph[p].children = intArrayAppend(graph[p].children, graph[p].childCount, i);
			graph[p].childCount += 1;
		}
		// Add edge (i-w,i) if the current node i>=w
		if (i >= w)
		{
			int p = i - w;
			// Add p to parent list of node i.
			graph[i].parents = intArrayAppend(graph[i].parents, graph[i].parentCount, p);
			graph[i].parentCount += 1;
			// Add node i to children list of p.
			graph[p].children = intArrayAppend(graph[p].children, graph[p].childCount, i);
			graph[p].childCount += 1;
		}


	}
	return graph;
}

// Sample iXOR Graph
struct bnode* generateiXORGraph(int n, int w, int seed)
{
	srand(seed);
	struct bnode* graph = malloc(n*sizeof(struct bnode));
	int i, j;
	int indeg = 20;
	graph[0].childCount = 0;
	graph[0].children = 0;
	graph[0].parentCount = 0;
	graph[0].parents = 0;

	for (i = 1; i < n; i++)
	{
		graph[i].childCount = 0;
		graph[i].peb = 0;
		graph[i].children = 0;
		graph[i].parentCount = 0;
		graph[i].parents = 0;
		int numParents = minimum(i, indeg);	//Number of parents for node i.
		// Populate parent list and add i to children list of its parents.
		for (j = 0; j < numParents; j++) {
			// Sample parent
			int modulus = (i < w) ? i : w;
			int p = i-(randomNumber() % modulus);
		    // Add p to parent list of node i.
			graph[i].parents = intArrayAppend(graph[i].parents, graph[i].parentCount, p);
			graph[i].parentCount += 1;
			// Add node i to children list of p.
			graph[p].children = intArrayAppend(graph[p].children, graph[p].childCount, i);
			graph[p].childCount += 1;
		}
	}
	return graph;
}

// Sample Balloon Hashing Linear Graph
struct bnodeCompressed* generateBHLinGraphNoMalloc(int n, int w, int seed)
{
	srand(seed);
	struct bnodeCompressed* graph =malloc(n*sizeof(struct bnodeCompressed));
	int i, j;
	int indeg = 20;

	graph[0].parentCount = 0;

    printf("starting to build the graph, n = %d\n",n);
	for (i = 1; i < n; i++)
	{
		graph[i].peb = 0;
		graph[i].parentCount = 0;
		int numParents = minimum(i, indeg);	//Number of parents for node i.
		// Populate parent list and add i to children list of its parents.
		graph[i].parentCount = numParents;
        graph[i].parentsSeedForRand = randomNumber();


	}

	return graph;
}
void resetGraph(struct bnode* G,int n){
	int i;
	for (i = 0; i < n; i++)
	{
		G[i].peb = 0;
	}
}
// Build node set S which reduces depth of graph.
struct sSet *selectPbSetS(int n, struct bnode* G, int parUBound, int segSize)
{
	int sizeOfLayer = maximum(1, minimum(n, parUBound * (segSize+1)));
	int numLayers = (int)ceil(n*1.0 / (1.0*sizeOfLayer));
	int * S = malloc(n*sizeof(int));
	memset(S, 0, n*sizeof(int));

	int i, j, k;
	int count = 0;

	for (i = 0; i < numLayers; i++)
	{
		int layerStart = i*sizeOfLayer;
		// Handle Last Layer (Special Size due to ceiling)
		int currentLayerSize = (i<numLayers - 1) ? sizeOfLayer : n - (numLayers - 1)*sizeOfLayer;
		for (j = 0; j < currentLayerSize; j++)
		{
			int added = 0;
			int curNode = layerStart + j;
			int segPos = nodePlaceInGap(curNode, sizeOfLayer, segSize); // j % (segSize + 1);
			if (segPos == 0 && j > 0) { S[curNode] = 1; added = 1; }  // "j>0" because no need to pebble first node in a layer

			for (k = 0; k < G[curNode].parentCount; k++)
			{
				int parent = G[curNode].parents[k];
				int parentSegPos = nodePlaceInGap(parent, sizeOfLayer, segSize); // ((parent) % sizeOfLayer) % (segSize + 1); // parents position in segment
				if (parent >= layerStart && segPos <= parentSegPos)        // ensure that edge either leaves current layer, or makes progress in segment position
				{
					// If edge doesn't leave layer and/or make progress then add to set S
					S[curNode] = 1;
					added = 1;
				}

			}
			count += added;
		}
	}

	printf("Size of S: %d\n", count);				// DEBUG CODE

	struct sSet *retVal = (struct sSet *)malloc(sizeof(struct sSet));
    retVal->numLayers   = numLayers;
	retVal->sizeOfLayer = sizeOfLayer;
	retVal->segSize     = segSize;
	retVal->S = S;

	return retVal;
}

double honestEComp(int n, int window,double R)
{
	double eComp= n*R + (window + 1.0)*(1.0*window) / 2.0 + (n - window)*1.0*window;
	return eComp;
}


int * selectbSetGapSCompressed(int gap, int numLayers, int n,int w, struct bnodeCompressed* G)
{
	//int numLayers = rootd;
	//int gap = rootd;
	int sizeOfLayer =  (int)ceil(n*1.0 / (numLayers*1.0));

	int * S = malloc(n*sizeof(int));
	memset(S, 0, n*sizeof(int));

	int i, j, k;
	int count = 0;

	for (i = 0; i < numLayers; i++)
	{
		int layerStart = i*sizeOfLayer;
		// Handle Last Layer (Special Size due to ceiling)
		int currentLayerSize = (i<numLayers - 1) ? sizeOfLayer : n - (numLayers - 1)*sizeOfLayer;
		for (j = 0; j < currentLayerSize; j++)
 	{
			int added = 0;
			int curNode = layerStart + j;
			int SegPos = nodePlaceInGap(curNode, sizeOfLayer, gap); // j % (rootd + 1);
			if (SegPos == 0 && j > 0) { S[curNode] = 1; added = 1; }  // "j>0" because no need to pebble first node in a layer
			srand(G[curNode].parentsSeedForRand);
			for (k = 0; k < G[curNode].parentCount; k++)
			{

				int modulus = (curNode < w) ? curNode : w;

				int parent =curNode-(randomNumber() % modulus);
				int parentSegPos = nodePlaceInGap(parent, sizeOfLayer, gap); // (parent % sizeOfLayer) % (rootd + 1);    // parents position in segment
				if (parent >= layerStart && SegPos <= parentSegPos)      // Ensure that edge either leaves current layer, or makes progress in segment position
				{
					// If pebble doesn't leave layer and/or make progress then add to set S
					S[curNode] = 1;
					added = 1;
				}

			}
			count += added;
		}
	}

	printf("Size of S: %d     Layers: %d     Gap: %d\n", count,numLayers,gap);                  //DEBUG CODE
	return S;
}

double honestECompParallelism(int noverp, int window,double R, int parallelism)
{
	printf("noverp-window/parallelism %d\n", noverp-window/parallelism);
   return noverp*parallelism*R+ noverp*1.0*window-(window*1.0/(1.0*parallelism))*window*1.0/2.0;
}

// Build node set S which reduces depth of graph.
int * selectbSetGapSParallelism(int gap, int numLayers, int nOverP, int parallelism, struct bnode* G)
{
	//int numLayers = rootd;
	//int gap = rootd;
	int sizeOfLayer =  (int)ceil(nOverP*1.0 / (numLayers*1.0));

	int * S = malloc(nOverP*parallelism*sizeof(int));
	memset(S, 0, nOverP*parallelism*sizeof(int));

	int i, j, k;
	int lane;
	int count = 0;

	for (i = 0; i < numLayers; i++)
	{
		int layerStart = i*sizeOfLayer;
		// Handle Last Layer (Special Size due to ceiling)
		int currentLayerSize = (i<numLayers - 1) ? sizeOfLayer : nOverP - (numLayers - 1)*sizeOfLayer;
		for (j = 0; j < currentLayerSize; j++)
		{
			for(lane=0; lane<parallelism; lane++)
			{
				int added = 0;
				int curNodePlaceInLane = layerStart + j;
				int curNode = curNodePlaceInLane*parallelism + lane;
				int SegPos = nodePlaceInGap(curNodePlaceInLane, sizeOfLayer, gap); // j % (rootd + 1);
				if (SegPos == 0 && j > 0) { S[curNode] = 1; added = 1; }  // "j>0" because no need to pebble first node in a layer

				for (k = 0; k < G[curNode].parentCount; k++)
				{
					int parent = G[curNode].parents[k];
					int parentLane = parent % parallelism;
					int parentPlaceInLane = (parent-parentLane)/parallelism;
					int parentSegPos = nodePlaceInGap(parentPlaceInLane, sizeOfLayer, gap); // (parent % sizeOfLayer) % (rootd + 1);    // parents position in segment
					if (parentPlaceInLane >= layerStart && SegPos <= parentSegPos)      // Ensure that edge either leaves current layer, or makes progress in segment position
					{
						// If pebble doesn't leave layer and/or make progress then add to set S
						S[curNode] = 1;
						added = 1;
					}

				}
				count += added;
			}
		}
	}

	printf("Size of S: %d     Layers: %d     Gap: %d\n", count,numLayers,gap);                  //DEBUG CODE
	return S;
}

// Build node set S which reduces depth of graph.
int * selectbSetGapS(int gap, int numLayers, int n, struct bnode* G)
{
	//int numLayers = rootd;
	//int gap = rootd;
	int sizeOfLayer =  (int)ceil(n*1.0 / (numLayers*1.0));

	int * S = malloc(n*sizeof(int));
	memset(S, 0, n*sizeof(int));

	int i, j, k;
	int count = 0;

	for (i = 0; i < numLayers; i++)
	{
		int layerStart = i*sizeOfLayer;
		// Handle Last Layer (Special Size due to ceiling)
		int currentLayerSize = (i<numLayers - 1) ? sizeOfLayer : n - (numLayers - 1)*sizeOfLayer;
		for (j = 0; j < currentLayerSize; j++)
		{
			int added = 0;
			int curNode = layerStart + j;
			int SegPos = nodePlaceInGap(curNode, sizeOfLayer, gap); // j % (rootd + 1);
			if (SegPos == 0 && j > 0) { S[curNode] = 1; added = 1; }  // "j>0" because no need to pebble first node in a layer

			for (k = 0; k < G[curNode].parentCount; k++)
			{
				int parent = G[curNode].parents[k];
				int parentSegPos = nodePlaceInGap(parent, sizeOfLayer, gap); // (parent % sizeOfLayer) % (rootd + 1);    // parents position in segment
				if (parent >= layerStart && SegPos <= parentSegPos)      // Ensure that edge either leaves current layer, or makes progress in segment position
				{
					// If pebble doesn't leave layer and/or make progress then add to set S
					S[curNode] = 1;
					added = 1;
				}

			}
			count += added;
		}
	}

	printf("Size of S: %d     Layers: %d     Gap: %d\n", count,numLayers,gap);                  //DEBUG CODE
	return S;
}


// Build node set S which reduces depth of graph.
int * selectbSetS(int rootd, int n, struct bnode* G)
{
	int numLayers = rootd;
	int gap = rootd;
	int sizeOfLayer =  (int)ceil(n*1.0 / (rootd*1.0));

	int * S = malloc(n*sizeof(int));
	memset(S, 0, n*sizeof(int));

	int i, j, k;
	int count = 0;

	for (i = 0; i < numLayers; i++)
	{
		int layerStart = i*sizeOfLayer;
		// Handle Last Layer (Special Size due to ceiling)
		int currentLayerSize = (i<numLayers - 1) ? sizeOfLayer : n - (numLayers - 1)*sizeOfLayer;
		for (j = 0; j < currentLayerSize; j++)
		{
			int added = 0;
			int curNode = layerStart + j;
			int SegPos = nodePlaceInGap(curNode, sizeOfLayer, gap); // j % (rootd + 1);
			if (SegPos == 0 && j > 0) { S[curNode] = 1; added = 1; }  // "j>0" because no need to pebble first node in a layer

			for (k = 0; k < G[curNode].parentCount; k++)
			{
				int parent = G[curNode].parents[k];
				int parentSegPos = nodePlaceInGap(parent, sizeOfLayer, gap); // (parent % sizeOfLayer) % (rootd + 1);    // parents position in segment
				if (parent >= layerStart && SegPos <= parentSegPos)      // Ensure that edge either leaves current layer, or makes progress in segment position
				{
					// If pebble doesn't leave layer and/or make progress then add to set S
					S[curNode] = 1;
					added = 1;
				}

			}
			count += added;
		}
	}

	printf("Size of S: %d\n", count);                  //DEBUG CODE
	return S;
}





//returns 1 if u is parent of a node v in [start,end] in G; otherwise 0
int IsbParent(struct bnode* G, int u, int start, int end)
{
	int numChildren = G[u].childCount;
	int j;
	for (j = 0; j < numChildren; j++)
	{
		if (G[u].children[j] <= end && G[u].children[j] >= start) return 1;
	}
	return 0;
}



// Output: the number of pebbles on nodes in the set Want(i,j) := "Parents(start,end) intersect [1, start-1] setminus S"
int bparents(struct bnode* G, int *S, int start, int end, int n)
{
	resetGraph(G, n);
	int size = 0;
	int i = start;
	int j;
	int parent;
	int*par = malloc(n* sizeof(int));   // Keeps track of parents already added. (In case u and v in [start, end] have common pebbled parent.)
	memset(par, 0, n * sizeof(int));
	int min = minimum(end, n - 1);
	for (i = start; i <= min; i++)
	{
		for (j = 0; j < G[i].parentCount; j++)
		{
			parent = G[i].parents[j];
			if (S[parent] == 0 && par[parent] != 1 && parent < start)
			{
				size++;
				G[i].peb = 1;
				par[i] = 1;            // Make sure dont double count this guy.
			}
		}

	}
	free(par);
	return size;
}

__inline  int nodeLayer(int node, int layerSize)
{
	return (int)floor(node*1.0 / (1.0*layerSize));
}

__inline  int nodePlaceInLayer(int node, int layerSize)
{
	return node- nodeLayer(node,layerSize)*layerSize;
}
// returns r s.t.  (node-r) is in the set S
// layerSize, gapInLayer (parameters which specify how S is constructed)
__inline  int nodePlaceInGap(int node, int layerSize, int gapInLayer)
{
	return nodePlaceInLayer(node, layerSize) % (gapInLayer + 1);
}


// Output: the number of pebbles on nodes in the set Want(j,end) := "Parents(j,end) intersect [1, j-1] setminus S"
// Precondition: size = |Want(j-1,end)|
int updatebParents(struct bnode* G, int *S, int size, int j, int end)
{
	int update=0;
	// Check if we need to add j.
	if (S[j] == 0 && G[j].peb ==0) {
		update++;
		G[j].peb = 1;
		
	}

	int k, l;
	//int* alreadyRemoved = malloc(sizeof(int) * n);
	// Iterate over parents of node j and check if they are still needed
	for (l = 0; l < G[j].parentCount; l++)
	{
		// Can we remove pebble from parent of j?
		int parent = G[j].parents[l];
		int count = G[parent].childCount;
		int *theChildren = G[parent].children;
		int canBeRemoved = 1;   //Flag denoting if current parent can be removed.

		// If parent l of node j is in S ignore it (as we account for S separately).
		if (S[parent] == 1 || G[parent].peb == 0)
		{
			canBeRemoved = 0;
		}
		else
		{
			// Check if any other children of parent(j) are in [j+1,end]. If so we still need a pebble on parent(j).
			for (k = 0; k < count; k++)
			{
				if (theChildren[k] > j && theChildren[k] <= end) 
				{
					canBeRemoved = 0; 
					break;
				}
			}
		}
		if (G[parent].peb == 1 && canBeRemoved == 1) {
			update--;
			G[parent].peb = 0;
		}

	}

	if (j > 0 && S[j - 1] == 0)
	{
		int parent = j - 1;
		int count = G[parent].childCount;
		int *theChildren = G[parent].children;
		int canBeRemoved = 1;   //Flag denoting if current parent can be removed.

		// If parent l of node j is in S ignore it (as we account for S separately).
		if (S[parent] == 1 || G[parent].peb == 0)
		{
			canBeRemoved = 0;
		}
		else
		{
			// Check if any other children of parent(j) are in [j+1,end]. If so we still need a pebble on parent(j).
			for (k = 0; k < count; k++)
			{
				if (theChildren[k] > j && theChildren[k] <= end) canBeRemoved = 0;
			}
		}
		if (G[parent].peb == 1 && canBeRemoved == 1) {
			update--;
			G[parent].peb = 0;
		}

	}
	
	return size + update;
}




void freeGraph(struct bnode* G, int n)
{
	int i;
	for (i = 0; i < n; i++)
	{
		int* children = G[i].children;
		int* parents = G[i].parents;
		if (children && G[i].childCount > 0) free(children);
		if (parents && G[i].parentCount > 0) free(parents);
	}
	free(G);
}

int parallelism(int layerSize, int gapInLayer)
{
	return 1 + (int)floor(layerSize*1.0 / (1.0*gapInLayer + 1.0));
}
