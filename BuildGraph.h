#ifndef BUILDGRAPH_H_INCLUDED
#define BUILDGRAPH_H_INCLUDED

// Type Forward Declarations
//struct node;
struct bnode;

// Type Forward Declarations
struct bnodeCompressed;

// Type Definitions: S set
struct sSet {
	int *S;
	int numLayers;
    int	sizeOfLayer;
    int segSize;
};

// Function Declarations

//________________________________________________________________________________________________________________________
// Argon2iA is the version of Argon2i that won the password hashing competition. Uniform edge distribution.
// Argon2iB is an updated version of Argon2i in response to recent attacks.
//      Non-uniform edge distribution + the label for node v+w is XORed with label v before the label v is
//       overwritted in memory
// iXOR is a hypothe
//______________________________________________

// Samples a  random (n/w)-pass Argon2iA graph with 1-lane
// int seed is the seed for the PRG
struct bnode* generateArgon2iAGraph(int n, int w, int seed);


// Samples a  random (n/w)-pass Argon2iA graph with parallelism lanes
// int seed is the seed for the PRG
struct bnode* generateArgon2iAGraphParallelism(int nOverP, int w, int seed, int parallelism);

// Generates a (n/w)-pass balloon hashing graph. Edge distributions are uniform, but
// the node v-w is also a parent of the node v
struct bnode* generateBHGraph(int n, int w, int seed);

// Samples a random (n/w)-pass Argon2iB graph with 1-lane, seed is the seed for the PRG
struct bnode* generateArgon2iBGraph(int n, int w, int seed);

// Samples a random (n/w)-pass Argon2iB graph with parallelism-lanes, seed is the seed for the PRG
struct bnode* generateArgon2iBGraphParallelism(int nOverP, int w, int seed, int parallelism);
struct bnode* generateiXORGraph(int n, int w, int seed);
struct bnodeCompressed* generateBHLinGraphNoMalloc(int n, int w, int seed);

// selects a depth reducing set S for the DAG G (assumed to be an iXOR instance)
//   divides graph into layers of (n/numLayers) consecutive nodes each
//   divides layers into segments of size gap+1
// Output: S such that depth(G-S) <= gap*numLayers
int * selectbSetGapS(int gapInLayer, int numLayers, int n, struct bnode* G);

// Similar to selectbSetGapS except that gapInLayer+1=numLayers=rootd
// Output: S such that  depth of G-S is at most (gap)*numlayers < d
int * selectSetS(int rootd, int n, struct node* G);
int * selectbSetS(int rootd, int n, struct bnode* G);

// selects a depth reducing set S for the DAG G (assumed to be an iXOR instance)
//   divides graph into layers of (n/numLayers) consecutive nodes each
//   divides layers into segments of size gap+1
// Output: S such that  depth of G-S is at most (gap)*numlayers
int * selectbSetGapSCompressed(int gap, int numLayers, int n,int w, struct bnodeCompressed* G);

// selects a depth reducing set S for the DAG G (assumed to be an Argon2i instance with parallelism different lanes)
// Output: S such that  depth of G-S is at most (gap)*numlayers
int * selectbSetGapSParallelism(int gap, int numLayers, int nOverP, int parallelism, struct bnode* G);

// Honest energy cost to compute Argon2i-A or Argon2i-B with n nodes and memory window w ... i.e. tau=(n/w)-pass Argon2i
// Assumes: Argon2i uses 1-lane
// R is the core energy/memory ratio (typical value R=3000)
double honestEComp(int n, int window,double R);

// Honest energy cost to compute Argon2i-A or Argon2i-B with n nodes and memory window w ... i.e. tau=(n/w)-pass Argon2i
// with #lanes = parallelism
//
// R is the core energy/memory ratio (typical value R=3000)
double honestECompParallelism(int noverp, int window,double R, int parallelism);


// Construct S when we have bounded parallelism (and hence can only reduce "exploitable depth" to d=n/parallelism)
// Output: sSet which contains
///          set S such that  depth of G-S is at most (gap)*numlayers, where gap+1 is implicitly defined as n/(parallelism*numLayers)
//           int parallelism
//           int gap
//           int numLayers
struct sSet *selectPbSetS(int n, struct bnode* G, int parallelism, int segSize);


// Local helper functions

int IsbParent(struct bnode* G, int u, int start, int end);


// Output: the number of pebbles on nodes in the set Want(start,end) := "Parents(PN(start,end)) intersect [1, j-1] setminus (S1 cup S2)" 
// here PN(start,end) denotes the nodes pebbled during the (outer) balloon phase between steps start and end
// start and end are defined by the current layer, current place in layer, and the end layer/endplace in endlayer. 
int bparents(struct bnode* G, int *S, int start, int end, int n);
// Output: the number of pebbles on nodes in the set Want(j,end) := "Parents(j,end) intersect [1, j-1] setminus S"
// Precondition: size = |Want(j-1,end)|
int updatebParents(struct bnode* G, int *S, int size, int j, int end);




// ================
// Helper Functions
//=================
//
// Returns minimum of two integer arguments
int minimum(int a, int b);
int maximum(int a, int b);

// Returns a random number in the range 0 to 2^30
// Necessary because the standard c random number generator returns numbers in the range 0 to 2^{15}
int randomNumber();

// returns new integer array consisting of src with newItem appended and and frees src
int * intArrayAppend(int * src, int len, int newItem);

// If we divide n nodes into (hypothetical) layers of size layerSize, which layer does int node fall into? Returns the answer.
__inline  int nodeLayer(int node, int layerSize);
// If we divide n nodes into (hypothetical) layers of size layerSize what is the index of int node in its layer? Returns the answer.
__inline  int nodePlaceInLayer(int node, int layerSize);
// If we divide n nodes into (hypothetical) layers of size layerSize, and we further divide  each layer into segments  of size gap+1.
// What is the index of int node in its segment? Returns the answer.
__inline  int nodePlaceInGap(int node, int layerSize, int gapInLayer);

// Given that we construct S with parameters layerSize and gapInLayer, how much parallelism is needed during balloon phase?
// Returns the answer.
int parallelism(int layerSize, int gapInLayer);


// Free's all allocated memory associated with G
void freeGraph(struct bnode* G, int n);
// reset G[i].peb = 0 for all i < n
void resetGraph(struct bnode* G,int n);

#endif
