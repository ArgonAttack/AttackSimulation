#ifndef ATTACK_H_INCLUDED
#define ATTACK_H_INCLUDED


// BHash Node
struct bnode {
	short peb;
	int parentCount;
	int* parents;
	int childCount;
	int* children;
};

// iXOR Node (saves memory by stored graph in compresses format) instead of explicitly storing pointers to 20 different parents, store the random
// seed necessary to recompute these parent pointers on the fly
struct bnodeCompressed {
	short peb;
	int parentCount;
	int parentsSeedForRand;
};

// Function Declarations

// Computes the cost of executing a single balloon phase of AB16 attack (with additional optimizations)
// Attack on Argon2i (version A or B) with 1-lane
// i current pebbling round when balloon phase starts
// S, sizeOfLAyer, gapInLayer are the parameters of the depth reducing set
// n number of nodes in the graph
// w memory window (e.g., (n/w)-passes through memory)
// R core memory-energy ratio which gives cost of computer hash function once (typically R=3000)
double balloonbPhaseCost(struct bnode* G, int* S, int sizeOfLayer, int gapInLayer, int n, int w, double R, int i);

// Computes the cost of executing a single balloon phase of AB16 attack (with additional optimizations)
// Attack is tailored to the iXOR construction with the compression trick.
// i current pebbling round when balloon phase starts
// S, sizeOfLAyer, gapInLayer are the parameters of the depth reducing set
// n number of nodes in the graph
// w memory window (e.g., (n/w)-passes through memory)
// R core memory-energy ratio which gives cost of computer hash function once (typically R=3000)
double balloonPhaseCostiXOR(struct bnodeCompressed* G, int* S, int sizeOfLayer, int gapInLayer, int n, int w, double R, int i);

// Computes the cost of executing a single balloon phase of AB16 attack (with additional optimizations)
// Attack on Argon2i (version A or B) with parallelism-lane
// i current pebbling round when balloon phase starts
// S, sizeOfLAyer, gapInLayer are the parameters of the depth reducing set
// n number of nodes in the graph
// w memory window (e.g., (n/w)-passes through memory)
// R core memory-energy ratio which gives cost of computer hash function once (typically R=3000)
double balloonPhaseCostParallel(struct bnode* G, int* S,  int sizeOfLayer, int gapInLayer, int nOverP, int w, double R, int startTime, int parallelism);

// Attack an Argon2i graph with parallelism lanes and memory window g
//  S, sizeOfLayer, gapInLayer define the depth-reducing set
//  g is an additional attack parameter
//  R is the core memory-energy ratio (i.e., cost to call compression function once)
// Simulates attack and returns an upper bound on the attack cost (upper bound is almost exact cost)
double attackParallelLanes(struct bnode* G, int* S, int sizeOfLayer, int gapInLayer, int nOverP, int g, int w, double R, int parallelism);

// Attack an Argon2i-A or Argon2i-B graph with memory window g
//  S, sizeOfLayer, gapInLayer define the depth-reducing set
//  g is an additional attack parameter
//  R is the core memory-energy ratio (i.e., cost to call compression function once)
// Simulates attack and returns an upper bound on the attack cost (upper bound is almost exact cost)
double attackBH(struct bnode* G, int* S, int sizeOfLayer, int gapInLayer, int n, int g, int w, double R);

// Attack an iXOR graph with memory window g
//  S, sizeOfLayer, gapInLayer define the depth-reducing set
//  g is an additional attack parameter
//  R is the core memory-energy ratio (i.e., cost to call compression function once)
// Simulates attack and returns an upper bound on the attack cost (upper bound is almost exact cost)
double attackiXOR(struct bnodeCompressed* G, int* S, int sizeOfLayer, int segSize, int n, int g, int w, double R);


#endif
