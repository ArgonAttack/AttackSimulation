#ifndef STATISTICS_H_INCLUDED
#define STATISTICS_H_INCLUDED

struct bnode;
struct bnodeCompressed;
struct sSet;

struct searchData {
	int g;
	int numLayers;
	int gapInLayer;
    int* S;
};

// Outputs vector of complexities of attack for increasingly tight bounds on parallelism
//  For each choice of parallelism do:
//   1) Choose initial numLayers (and set S).
//   2) Repeat 3 times:
//     A) Find optimal g for current numLayers.
//     B) Find optimal numLayers (and set S) for current g.
//   3) Return final complexity
void parStats(struct bnode *G, int n, int window, double R, int parUBound);

// Searches for optimal value of g for attack on graph G with set S
// Returns optimal g
int searchbForg(struct bnode* G, int* S, int n, int window, double R, int sizeOfLayer, int gapInLayer, int granularity);
int searchbForgCompressed(struct bnodeCompressed* G, int* S, int n, int window, double R, int sizeOfLayer, int segSize, int granularity);


// Deprecated: Unlike searchForg, searchFord generates the set S and returns optimal d, assumes gap+1=numlayers=rootd
// Returns optimal value of rootd

int searchbFord(struct bnode* G, int ** Sout, int n, int window, double R, int g);

struct searchData *searchbForBestDepthParameters(struct bnode* G,   int n, int window, double R, int g);
struct searchData *searchbForBestDepthParametersCompressed(struct bnodeCompressed* G,   int n, int window, double R, int g);
struct searchData *searchbForBestDepthParametersBddParallelism(struct bnode* G,  int n, int window, double R, int g,int parBound);


// Unlike searchForg, searchForNumLayers generates the set S and returns optimal numLayers
// Returns optimal value of numLayers
int searchForNumLayers(struct node* G, int ** S, int n, int window, double R, int g, int parUBound);
struct sSet *searchbForSegSize(struct bnode* G, int n, int window, double R, int g, int parUBound);



// Returns the minimum (resp. maximum) of two integers a and b
int maximum(int a, int b);
int minimum(int a, int b);


#endif
