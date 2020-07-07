#include <stdio.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <stdlib.h>
#include <string.h>

//#define VERBOSE
#define PROGRESS
#define RECOVER_FLOW
#define STATS

#define EPSILON 1e-10
#define INFINITY 1.79769e+308

typedef long long int llint;

static char *
getNextWord (char *line, char *word)
{
	int wordlength = 0;

	while (((*line) == ' ') || ((*line) == '\t') || ((*line) == '\n') || ((*line) == '\0'))
	{
		++ line;
	}

	while (((*line) != ' ') && ((*line) != '\t') && ((*line) != '\n') && ((*line) != '\0'))
	{
		word[wordlength] = (*line);
		++ wordlength;
		++ line; 
	}
	
	word[wordlength] = '\0';
	return line;
}

double 
timer (void)
{
  struct rusage r;

  getrusage(0, &r);
  return (double) (r.ru_utime.tv_sec + r.ru_utime.tv_usec / (double)1000000);
}

struct node;

typedef struct arc 
{
	struct node *from;
	struct node *to;
	double flow;
	double capacity;
    double base_capacity;
	int direction;
	//int *capacities;
} Arc;

typedef struct edge 
{
	int from;
    int to;
	double capacity;
} Edge;

typedef struct node 
{
	int visited;
	int numAdjacent;
	int number;
	int label;
	double excess;
	struct node *parent;
	struct node *childList;
	struct node *nextScan;
	int numOutOfTree;
	Arc **outOfTree;
	int nextArc;
	Arc *arcToParent;
	struct node *next;
	struct node *prev;
	//int breakpoint;
} Node;


typedef struct root 
{
	Node *start;
	Node *end;
} Root;

//---------------  Global variables ------------------
static int numNodes = 0;
static int numArcs = 0;
static int source = 0;
static int sink = 0;
static Edge *initEdgeList = NULL;



static int numParams = 0;

static int highestStrongLabel = 1;

static Node *adjacencyList = NULL;
static Root *strongRoots = NULL;
static int *labelCount = NULL;
static Arc *arcList = NULL;
//-----------------------------------------------------

#ifdef STATS
static llint numPushes = 0;
static int numMergers = 0;
static int numRelabels = 0;
static int numGaps = 0;
static llint numArcScans = 0;
#endif

static void
initializeNode (Node *nd, const int n)
{
	nd->label = 0;
	nd->excess = 0;
	nd->parent = NULL;
	nd->childList = NULL;
	nd->nextScan = NULL;
	nd->nextArc = 0;
	nd->numOutOfTree = 0;
	nd->arcToParent = NULL;
	nd->next = NULL;
	nd->prev = NULL;
	nd->visited = 0;
	nd->numAdjacent = 0;
	nd->number = n;
	nd->outOfTree = NULL;
	//nd->breakpoint = (numParams+1);
}

static void
initializeRoot (Root *rt) 
{
	rt->start = (Node *) malloc (sizeof(Node));
	rt->end = (Node *) malloc (sizeof(Node));

	if ((rt->start == NULL) || (rt->end == NULL))
	{
		printf ("%s Line %d: Out of memory\n", __FILE__, __LINE__);
		exit (1);
	}

	initializeNode (rt->start, 0);
	initializeNode (rt->end, 0);

	rt->start->next = rt->end;
	rt->end->prev = rt->start;
}


static void
freeRoot (Root *rt) 
{
	free(rt->start);
	rt->start = NULL;

	free(rt->end);
	rt->end = NULL;
}

static void
//liftAll (Node *rootNode, const int theparam) 
liftAll (Node *rootNode) 
{
	Node *temp, *current=rootNode;

	current->nextScan = current->childList;

	-- labelCount[current->label];
	current->label = numNodes;	
	//current->breakpoint = (theparam+1);

	for ( ; (current); current = current->parent)
	{
		while (current->nextScan) 
		{
			temp = current->nextScan;
			current->nextScan = current->nextScan->next;
			current = temp;
			current->nextScan = current->childList;

			-- labelCount[current->label];
			current->label = numNodes;
		}
	}
}

static void
addToStrongBucket (Node *newRoot, Node *rootEnd) 
{
	newRoot->next = rootEnd;
	newRoot->prev = rootEnd->prev;
	rootEnd->prev = newRoot;
	newRoot->prev->next = newRoot;
}

static void
createOutOfTree (Node *nd)
{
	if (nd->numAdjacent)
	{
		if ((nd->outOfTree = (Arc **) malloc (nd->numAdjacent * sizeof (Arc *))) == NULL)
		{
			printf ("%s Line %d: Out of memory\n", __FILE__, __LINE__);
			exit (1);
		}
	}
}

static void
initializeArc (Arc *ac)
{
	//int i;

	ac->from = NULL;
	ac->to = NULL;
	ac->capacity = 0;
    ac->base_capacity = 0;
	ac->flow = 0;
	ac->direction = 1;
}

static void
addOutOfTreeNode (Node *n, Arc *out) 
{
	n->outOfTree[n->numOutOfTree] = out;
	++ n->numOutOfTree;
}

static void
reCreateGraph(const int num_edges, const double lambda, const double c) 
{
	int i, capacity, from, to, first=0;
	Arc *ac = NULL;


	if ((adjacencyList = (Node *) malloc (numNodes * sizeof (Node))) == NULL)
	{
		printf ("%s, %d: Could not allocate memory.\n", __FILE__, __LINE__);
		exit (1);
	}

	if ((strongRoots = (Root *) malloc (numNodes * sizeof (Root))) == NULL)
	{
		printf ("%s, %d: Could not allocate memory.\n", __FILE__, __LINE__);
		exit (1);
	}

	if ((labelCount = (int *) malloc (numNodes * sizeof (int))) == NULL)
	{
		printf ("%s, %d: Could not allocate memory.\n", __FILE__, __LINE__);
		exit (1);
	}

	if ((arcList = (Arc *) malloc (numArcs * sizeof (Arc))) == NULL)
	{
		printf ("%s, %d: Could not allocate memory.\n", __FILE__, __LINE__);
		exit (1);
	}
          
	for (i=0; i<numNodes; ++i)
	{
		initializeRoot (&strongRoots[i]);
		initializeNode (&adjacencyList[i], (i+1));
		labelCount[i] = 0;
	}

	for (i=0; i<numArcs; ++i)
	{
		initializeArc (&arcList[i]);
	}
    for (i=0; i<numArcs; ++i)
    {

        from = initEdgeList[i].from;		
        to = initEdgeList[i].to;

        ac = &arcList[first];

        ac->from = &adjacencyList[from-1];
        ac->to = &adjacencyList[to-1];
              
        
        ac->capacity = initEdgeList[i].capacity;
        ac->base_capacity = ac->capacity;
        
        if (from == source)
        {
            if (to <= num_edges)
            {
                ac->capacity = c;
            } else {
                ac->capacity = lambda;
            }
        }
        ++ first;

        ++ ac->from->numAdjacent;
        ++ ac->to->numAdjacent;
    }
		

	for (i=0; i<numNodes; ++i) 
	{
		createOutOfTree (&adjacencyList[i]);
	}

	for (i=0; i<numArcs; i++) 
	{
		to = arcList[i].to->number;
		from = arcList[i].from->number;
		capacity = arcList[i].capacity;

		if (!((source == to) || (sink == from) || (from == to))) 
		{
			if ((source == from) && (to == sink)) 
			{
				arcList[i].flow = capacity;
			}
			else if (from == source)
			{
				addOutOfTreeNode (&adjacencyList[from-1], &arcList[i]);
			}
			else if (to == sink)
			{
				addOutOfTreeNode (&adjacencyList[to-1], &arcList[i]);
			}
			else
			{
				addOutOfTreeNode (&adjacencyList[from-1], &arcList[i]);
			}
		}
	}
}


static void
readDimacsFileCreateList (const char *filename) 
{
	int lineLength=32768, i, capacity, numLines = 0, from, to, first=0;
	char *line, *word, ch, ch1, *tmpline;
	Arc *ac = NULL;
    printf ("Read filename %s\n", filename);
    FILE *file = fopen(filename, "r");

	if ((line = (char *) malloc ((lineLength+1) * sizeof (char))) == NULL)
	{
		printf ("%s, %d: Could not allocate memory.\n", __FILE__, __LINE__);
		exit (1);
	}

	if ((word = (char *) malloc ((lineLength+1) * sizeof (char))) == NULL)
	{
		printf ("%s, %d: Could not allocate memory.\n", __FILE__, __LINE__);
		exit (1);
	}

	while (fgets (line, lineLength, file))
	{
		++ numLines;

#ifdef VERBOSE
		printf ("Read line %s\n", line);
#endif

		switch (*line)
		{
		case 'p':

			//sscanf (line, "%c %s %d %d %d", &ch, word, &numNodes, &numArcs, &numParams);
            sscanf (line, "%c %s %d %d", &ch, word, &numNodes, &numArcs);

#ifdef VERBOSE
			printf ("numNodes: %d, numArcs: %d, numParams: %d\n", 
				numNodes,
				numArcs,
				numParams);
#endif

			if ((adjacencyList = (Node *) malloc (numNodes * sizeof (Node))) == NULL)
			{
				printf ("%s, %d: Could not allocate memory.\n", __FILE__, __LINE__);
				exit (1);
			}

			if ((strongRoots = (Root *) malloc (numNodes * sizeof (Root))) == NULL)
			{
				printf ("%s, %d: Could not allocate memory.\n", __FILE__, __LINE__);
				exit (1);
			}

			if ((labelCount = (int *) malloc (numNodes * sizeof (int))) == NULL)
			{
				printf ("%s, %d: Could not allocate memory.\n", __FILE__, __LINE__);
				exit (1);
			}

			if ((arcList = (Arc *) malloc (numArcs * sizeof (Arc))) == NULL)
			{
				printf ("%s, %d: Could not allocate memory.\n", __FILE__, __LINE__);
				exit (1);
			}
            
            if ((initEdgeList = (Edge *) malloc (numArcs * sizeof (Edge))) == NULL)
			{
				printf ("%s, %d: Could not allocate memory.\n", __FILE__, __LINE__);
				exit (1);
			}
           
			for (i=0; i<numNodes; ++i)
			{
				initializeRoot (&strongRoots[i]);
				initializeNode (&adjacencyList[i], (i+1));
				labelCount[i] = 0;
			}

			for (i=0; i<numArcs; ++i)
			{
				initializeArc (&arcList[i]);
			}

			break;

		case 'a':

			tmpline = line;
			++ tmpline;

			tmpline = getNextWord (tmpline, word);			
			from = (int) atoi (word);

			tmpline = getNextWord (tmpline, word);			
			to = (int) atoi (word);

			ac = &arcList[first];

			ac->from = &adjacencyList[from-1];
			ac->to = &adjacencyList[to-1];
            
            
            tmpline = getNextWord (tmpline, word);
			ac->capacity = (double) atof (word);
            ac->base_capacity = ac->capacity;
            
            initEdgeList[first].from = from;
            initEdgeList[first].to = to;
            initEdgeList[first].capacity = ac->capacity;

			++ first;

			++ ac->from->numAdjacent;
			++ ac->to->numAdjacent;

			break;

		case 'n':

			sscanf (line, "%c %d %c", &ch, &i, &ch1);

			if (ch1 == 's')
			{
				source = i;	
			}
			else if (ch1 == 't')
			{
				sink = i;	
			}
			else
			{
				printf ("Unrecognized character %c on line %d\n", ch1, numLines);
				exit (1);
			}

			break;
		}
	}

	for (i=0; i<numNodes; ++i) 
	{
		createOutOfTree (&adjacencyList[i]);
	}

	for (i=0; i<numArcs; i++) 
	{
		to = arcList[i].to->number;
		from = arcList[i].from->number;
		capacity = arcList[i].capacity;

		if (!((source == to) || (sink == from) || (from == to))) 
		{
			if ((source == from) && (to == sink)) 
			{
				arcList[i].flow = capacity;
			}
			else if (from == source)
			{
				addOutOfTreeNode (&adjacencyList[from-1], &arcList[i]);
			}
			else if (to == sink)
			{
				addOutOfTreeNode (&adjacencyList[to-1], &arcList[i]);
			}
			else
			{
				addOutOfTreeNode (&adjacencyList[from-1], &arcList[i]);
			}
		}
	}

	free (line);
	line = NULL;

	free (word);
	word = NULL;
    fclose(file);
}

static void
simpleInitialization (void) 
{
	int i, size;
	Arc *tempArc;

	size = adjacencyList[source-1].numOutOfTree;
	for (i=0; i<size; ++i) 
	{
		tempArc = adjacencyList[source-1].outOfTree[i];
		tempArc->flow = tempArc->capacity;
		tempArc->to->excess += tempArc->capacity;
	}

	size = adjacencyList[sink-1].numOutOfTree;
	for (i=0; i<size; ++i)
	{
		tempArc = adjacencyList[sink-1].outOfTree[i];
		tempArc->flow = tempArc->capacity;
		tempArc->from->excess -= tempArc->capacity;
	}

	adjacencyList[source-1].excess = 0;
	adjacencyList[sink-1].excess = 0;

	for (i=0; i<numNodes; ++i) 
	{
		if (adjacencyList[i].excess > 0)
        //if (adjacencyList[i].excess > EPSILON) 
		{
		    adjacencyList[i].label = 1;
			++ labelCount[1];

			addToStrongBucket (&adjacencyList[i], strongRoots[1].end);
		}
	}

	adjacencyList[source-1].label = numNodes;
	adjacencyList[sink-1].label = 0;
	labelCount[0] = (numNodes - 2) - labelCount[1];
}

static void
reInitialization_negative (double lambda, double initial_guess) 
{
	int i, size;
	Arc *tempArc, *ac;
    Node *nd;    
    
    for (i=0; i<numNodes; ++i)
	{
		initializeRoot (&strongRoots[i]);
        
        nd = &adjacencyList[i];    
        nd->label = 0;
        nd->excess = 0;
        nd->parent = NULL;
        nd->childList = NULL;
        nd->nextScan = NULL;
        nd->nextArc = 0;
        //nd->numOutOfTree = 0;
        nd->arcToParent = NULL;
        nd->next = NULL;
        nd->prev = NULL;
        nd->visited = 0; 
		
		labelCount[i] = 0;
	}

	for (i=0; i<numArcs; ++i)
	{
		ac = &arcList[i];    
        ac->flow = 0;
        ac->direction = 1;
	}
    
     
	size = adjacencyList[source-1].numOutOfTree;
	for (i=0; i<size; ++i) 
	{
		tempArc = adjacencyList[source-1].outOfTree[i];
        tempArc->capacity = lambda;
		tempArc->flow = tempArc->capacity;
		tempArc->to->excess += tempArc->capacity;
        printf ("c Reinit: source-adjacent arc (%d, %d): capacity %.3f.\n",
				tempArc->from->number,
				tempArc->to->number,
				tempArc->capacity);
	}

	size = adjacencyList[sink-1].numOutOfTree;
	for (i=0; i<size; ++i)
	{
		tempArc = adjacencyList[sink-1].outOfTree[i];
        tempArc->capacity = tempArc->base_capacity - initial_guess;
		tempArc->flow = tempArc->capacity;
		tempArc->from->excess -= tempArc->capacity;
        printf ("c Reinit: sink-adjacent arc (%d, %d): capacity %.3f.\n",
				tempArc->from->number,
				tempArc->to->number,
				tempArc->capacity);
	}

	adjacencyList[source-1].excess = 0;
	adjacencyList[sink-1].excess = 0;

	for (i=0; i<numNodes; ++i) 
	{
		if (adjacencyList[i].excess > 0) 
		{
		    adjacencyList[i].label = 1;
			++ labelCount[1];

			addToStrongBucket (&adjacencyList[i], strongRoots[1].end);
		}
	}

	adjacencyList[source-1].label = numNodes;
	adjacencyList[sink-1].label = 0;
	labelCount[0] = (numNodes - 2) - labelCount[1];
}

static void
reInitialization (double lambda, double initial_guess, int num_edges) 
{
	int i, size;
	Arc *tempArc, *ac;
    Node *nd;    
    
    for (i=0; i<numNodes; ++i)
	{
		initializeRoot (&strongRoots[i]);
        
        nd = &adjacencyList[i];    
        nd->label = 0;
        nd->excess = 0;
        nd->parent = NULL;
        nd->childList = NULL;
        nd->nextScan = NULL;
        nd->nextArc = 0;
        nd->arcToParent = NULL;
        nd->next = NULL;
        nd->prev = NULL;
        nd->visited = 0;
		
		labelCount[i] = 0;
	}

	for (i=0; i<numArcs; ++i)
	{
		ac = &arcList[i];    
        ac->flow = 0;
        ac->direction = 1;
	}
    
     
	size = adjacencyList[source-1].numOutOfTree;
	for (i=0; i<size; ++i) 
	{
		tempArc = adjacencyList[source-1].outOfTree[i];
        if (tempArc->to->number <= num_edges) {
            tempArc->capacity = initial_guess;
        } else {
            tempArc->capacity = lambda;
        }
		tempArc->flow = tempArc->capacity;
		tempArc->to->excess += tempArc->capacity;
        printf ("c Reinit: source-adjacent arc (%d, %d): capacity %.3f.\n",
				tempArc->from->number,
				tempArc->to->number,
				tempArc->capacity);
	}

	size = adjacencyList[sink-1].numOutOfTree;
	for (i=0; i<size; ++i)
	{
		tempArc = adjacencyList[sink-1].outOfTree[i];
		tempArc->flow = tempArc->capacity;
		tempArc->from->excess -= tempArc->capacity;
        printf ("c Reinit: sink-adjacent arc (%d, %d): capacity %.3f.\n",
				tempArc->from->number,
				tempArc->to->number,
				tempArc->capacity);
	}

	adjacencyList[source-1].excess = 0;
	adjacencyList[sink-1].excess = 0;

	for (i=0; i<numNodes; ++i) 
	{
		if (adjacencyList[i].excess > 0)
		{
		    adjacencyList[i].label = 1;
			++ labelCount[1];

			addToStrongBucket (&adjacencyList[i], strongRoots[1].end);
		}
	}

	adjacencyList[source-1].label = numNodes;
	adjacencyList[sink-1].label = 0;
	labelCount[0] = (numNodes - 2) - labelCount[1];
}

static inline int 
addRelationship (Node *newParent, Node *child) 
{
	child->parent = newParent;
	child->next = newParent->childList;
	newParent->childList = child;

	return 0;
}

static inline void
breakRelationship (Node *oldParent, Node *child) 
{
	Node *current;

	child->parent = NULL;

	if (oldParent->childList == child) 
	{
		oldParent->childList = child->next;
		child->next = NULL;
		return;
	}

	for (current = oldParent->childList; (current->next != child); current = current->next);

	current->next = child->next;
	child->next = NULL;
}

static void
merge (Node *parent, Node *child, Arc *newArc) 
{
	Arc *oldArc;
	Node *current = child, *oldParent, *newParent = parent;

#ifdef STATS
	++ numMergers;
#endif

	while (current->parent) 
	{
		oldArc = current->arcToParent;
		current->arcToParent = newArc;
		oldParent = current->parent;
		breakRelationship (oldParent, current);
		addRelationship (newParent, current);
		newParent = current;
		current = oldParent;
		newArc = oldArc;
		newArc->direction = 1 - newArc->direction;
	}

	current->arcToParent = newArc;
	addRelationship (newParent, current);
}


static inline void 
pushUpward (Arc *currentArc, Node *child, Node *parent, const double resCap) 
{
#ifdef STATS
	++ numPushes;
#endif

	if (resCap >= child->excess) 
	{
		parent->excess += child->excess;
		currentArc->flow += child->excess;
		child->excess = 0;
		return;
	}

	currentArc->direction = 0;
	parent->excess += resCap;
	child->excess -= resCap;
	currentArc->flow = currentArc->capacity;
	parent->outOfTree[parent->numOutOfTree] = currentArc;
	++ parent->numOutOfTree;
	breakRelationship (parent, child);

	addToStrongBucket (child, strongRoots[child->label].end);
}


static inline void
pushDownward (Arc *currentArc, Node *child, Node *parent, double flow) 
{
#ifdef STATS
	++ numPushes;
#endif

	if (flow >= child->excess) 
	{
		parent->excess += child->excess;
		currentArc->flow -= child->excess;
		child->excess = 0;
		return;
	}

	currentArc->direction = 1;
	child->excess -= flow;
	parent->excess += flow;
	currentArc->flow = 0;
	parent->outOfTree[parent->numOutOfTree] = currentArc;
	++ parent->numOutOfTree;
	breakRelationship (parent, child);

	addToStrongBucket (child, strongRoots[child->label].end);
}

static void
pushExcess (Node *strongRoot) 
{
	Node *current, *parent;
	Arc *arcToParent;

	for (current = strongRoot; (current->excess && current->parent); current = parent) 
	{
		parent = current->parent;
		arcToParent = current->arcToParent;
		if (arcToParent->direction)
		{
			pushUpward (arcToParent, current, parent, (arcToParent->capacity - arcToParent->flow)); 
		}
		else
		{
			pushDownward (arcToParent, current, parent, arcToParent->flow); 
		}
	}

	if (current->excess > 0)
    //if (current->excess > EPSILON) 
	{
		if (!current->next)
		{
			addToStrongBucket (current, strongRoots[current->label].end);
		}
	}
}


static Arc *
findWeakNode (Node *strongNode, Node **weakNode) 
{
	int i, size;
	Arc *out;

	size = strongNode->numOutOfTree;

	for (i=strongNode->nextArc; i<size; ++i) 
	{

#ifdef STATS
		++ numArcScans;
#endif

		if (strongNode->outOfTree[i]->to->label == (highestStrongLabel-1)) 
		{
			strongNode->nextArc = i;
			out = strongNode->outOfTree[i];
			(*weakNode) = out->to;
			-- strongNode->numOutOfTree;
			strongNode->outOfTree[i] = strongNode->outOfTree[strongNode->numOutOfTree];
			return (out);
		}
		else if (strongNode->outOfTree[i]->from->label == (highestStrongLabel-1)) 
		{
			strongNode->nextArc = i;
			out = strongNode->outOfTree[i];
			(*weakNode) = out->from;
			-- strongNode->numOutOfTree;
			strongNode->outOfTree[i] = strongNode->outOfTree[strongNode->numOutOfTree];
			return (out);
		}
	}

	strongNode->nextArc = strongNode->numOutOfTree;

	return NULL;
}


static void
checkChildren (Node *curNode) 
{
	for ( ; (curNode->nextScan); curNode->nextScan = curNode->nextScan->next)
	{
		if (curNode->nextScan->label == curNode->label)
		{
			return;
		}
		
	}	

	-- labelCount[curNode->label];
	++	curNode->label;
	++ labelCount[curNode->label];

#ifdef STATS
	++ numRelabels;
#endif

	curNode->nextArc = 0;
}

static void
processRoot (Node *strongRoot) 
{
	Node *temp, *strongNode = strongRoot, *weakNode;
	Arc *out;

	strongRoot->nextScan = strongRoot->childList;

	if ((out = findWeakNode (strongRoot, &weakNode)))
	{
		merge (weakNode, strongNode, out);
		pushExcess (strongRoot);
		return;
	}

	checkChildren (strongRoot);
	
	while (strongNode)
	{
		while (strongNode->nextScan) 
		{
			temp = strongNode->nextScan;
			strongNode->nextScan = strongNode->nextScan->next;
			strongNode = temp;
			strongNode->nextScan = strongNode->childList;

			if ((out = findWeakNode (strongNode, &weakNode)))
			{
				merge (weakNode, strongNode, out);
				pushExcess (strongRoot);
				return;
			}

			checkChildren (strongNode);
		}

		if ((strongNode = strongNode->parent))
		{
			checkChildren (strongNode);
		}
	}

	addToStrongBucket (strongRoot, strongRoots[strongRoot->label].end);

	++ highestStrongLabel;
}

static Node *
getHighestStrongRoot () 
{
	int i;
	Node *strongRoot;

	for (i=highestStrongLabel; i>0; --i) 
	{
		if (strongRoots[i].start->next != strongRoots[i].end)  
		{
			highestStrongLabel = i;
			if (labelCount[i-1]) 
			{
				strongRoot = strongRoots[i].start->next;
				strongRoot->next->prev = strongRoot->prev;
				strongRoot->prev->next = strongRoot->next;
				strongRoot->next = NULL;
				return strongRoot;				
			}

			while (strongRoots[i].start->next != strongRoots[i].end) 
			{

#ifdef STATS
				++ numGaps;
#endif
				strongRoot = strongRoots[i].start->next;
				strongRoot->next->prev = strongRoot->prev;
				strongRoot->prev->next = strongRoot->next;
				liftAll (strongRoot);
			}
		}
	}

	if (strongRoots[0].start->next == strongRoots[0].end) 
	{
		return NULL;
	}

	while (strongRoots[0].start->next != strongRoots[0].end) 
	{
		strongRoot = strongRoots[0].start->next;
		strongRoot->next->prev = strongRoot->prev;
		strongRoot->prev->next = strongRoot->next;

		strongRoot->label = 1;
		-- labelCount[0];
		++ labelCount[1];

#ifdef STATS
		++ numRelabels;
#endif

		addToStrongBucket (strongRoot, strongRoots[strongRoot->label].end);
	}	

	highestStrongLabel = 1;

	strongRoot = strongRoots[1].start->next;
	strongRoot->next->prev = strongRoot->prev;
	strongRoot->prev->next = strongRoot->next;
	strongRoot->next = NULL;

	return strongRoot;	
}

static double getNewCapacity(const double capParameter, const Arc *tempArc)
{
    double newcapacity = -1;
    if (tempArc->from->number == source)
    {
        newcapacity = capParameter;
    }
    else
    {
        if (tempArc->to->number == sink)
        {
            newcapacity = tempArc->capacity;
        }
    }
    
    return newcapacity;
}

static void
updateCapacities (const double capParameter)
{
	int i, size;
	double delta, newcapacity;
	Arc *tempArc;
    
	size = adjacencyList[source-1].numOutOfTree;
	for (i=0; i<size; ++i) 
	{
		
        tempArc = adjacencyList[source-1].outOfTree[i];
       
		newcapacity = getNewCapacity(capParameter, tempArc);
        delta = (newcapacity - tempArc->capacity);
               
		if (delta < 0)
		{
			printf ("c Error on source-adjacent arc (%d, %d): capacity decreases by %.3f.\n",
				tempArc->from->number,
				tempArc->to->number,
				(-delta));
			exit(0);
		}

		tempArc->capacity += delta;
		tempArc->flow += delta;
		tempArc->to->excess += delta;

		if ((tempArc->to->label < numNodes) && (tempArc->to->excess > 0))
		{
			pushExcess (tempArc->to);
		}
	}

	size = adjacencyList[sink-1].numOutOfTree;
	for (i=0; i<size; ++i)
	{
		tempArc = adjacencyList[sink-1].outOfTree[i];
        newcapacity = getNewCapacity(capParameter, tempArc);
        delta = (newcapacity - tempArc->capacity);       
       
		if (delta > 0)
		{
			printf ("c Error on sink-adjacent arc (%d, %d): capacity %.3f increases to %.3f.\n",
				tempArc->from->number,
				tempArc->to->number,
				tempArc->capacity,
				newcapacity);
			exit(0);
		}

		tempArc->capacity += delta;
		tempArc->flow += delta;
		tempArc->from->excess -= delta;

		if ((tempArc->from->label < numNodes) && (tempArc->from->excess > 0))
		{
			pushExcess (tempArc->from);
		}
	}

	highestStrongLabel = (numNodes-1);
}


static void
updateSrcCapacities (const double newcapacity, const int num_edges)
{
	int i, size;
	double delta;
	Arc *tempArc;
    
	size = adjacencyList[source-1].numOutOfTree;
	for (i=0; i<size; ++i) 
	{
		tempArc = adjacencyList[source-1].outOfTree[i];
        if (tempArc->to->number <= num_edges)
        {
            
            delta = (newcapacity - tempArc->capacity);
                   
            if (delta < -EPSILON)
            {
                printf ("c Error on source-adjacent arc (%d, %d): capacity decreases by %.13f.\n",
                    tempArc->from->number,
                    tempArc->to->number,
                    (-delta));
                exit(0);
            }

            tempArc->capacity += delta;
            tempArc->flow += delta;
            tempArc->to->excess += delta;

            if ((tempArc->to->label < numNodes) && (tempArc->to->excess > 0))
            {
                pushExcess(tempArc->to);
            }
        }
	}

	

	highestStrongLabel = (numNodes-1);
}

static void
updateSinkCapacities (const double delta)
{
	int i, size;	
	Arc *tempArc;

	size = adjacencyList[sink-1].numOutOfTree;
	for (i=0; i<size; ++i)
	{
		tempArc = adjacencyList[sink-1].outOfTree[i];       
		if (delta > 0)
		{
			printf ("c Error on sink-adjacent arc (%d, %d): capacity %.3f increases to %.3f.\n",
				tempArc->from->number,
				tempArc->to->number,
				tempArc->capacity,
				tempArc->capacity + delta);
			exit(0);
		}
        printf ("c sink-adjacent arc (%d, %d): capacity %.3f increases   to %.3f.\n",
				tempArc->from->number,
				tempArc->to->number,
				tempArc->capacity,
				tempArc->capacity + delta);
		tempArc->capacity += delta;
		tempArc->flow += delta;
		tempArc->from->excess -= delta;

		if ((tempArc->from->label < numNodes) && (tempArc->from->excess > 0))
		{
			pushExcess (tempArc->from);
		}
	}

	highestStrongLabel = (numNodes-1);
}




static double
computeMinCut (void)
{
	int i;
    double mincut = 0;

	for (i=0; i<numArcs; ++i) 
	{
		
        if ((arcList[i].from->label >= numNodes) && (arcList[i].to->label < numNodes))
		{
			
            mincut += arcList[i].capacity;
		}
	}
	return mincut;
}

static void
pseudoflowPhase1 (void) 
{
	Node *strongRoot;
	//int theparam = 0;
	double thetime;
    
	thetime = timer ();
    
	while ((strongRoot = getHighestStrongRoot()))  
	{ 
		processRoot (strongRoot);
	}
}

static void
checkOptimality (void) 
{
	int i, check = 1;
	double mincut = 0, *excess;
    double difference;

	excess = (double *) malloc (numNodes * sizeof (double));
	if (!excess)
	{
		printf ("%s Line %d: Out of memory\n", __FILE__, __LINE__);
		exit (1);
	}

	for (i=0; i<numNodes; ++i)
	{
		excess[i] = 0;
	}

	for (i=0; i<numArcs; ++i) 
	{
		if ((arcList[i].from->label >= numNodes) && (arcList[i].to->label < numNodes))
		{
			mincut += arcList[i].capacity;
		}

		if ((arcList[i].flow > arcList[i].capacity) || (arcList[i].flow < 0)) 
		{
			check = 0;
			printf("c Capacity constraint violated on arc (%d, %d)\n", 
				arcList[i].from->number,
				arcList[i].to->number);
		}
		excess[arcList[i].from->number - 1] -= arcList[i].flow;
		excess[arcList[i].to->number - 1] += arcList[i].flow;
	}
    
    printf("c Mincut: %.3f\n", mincut);

	for (i=0; i<numNodes; i++) 
	{
		if ((i != (source-1)) && (i != (sink-1))) 
		{
			if ((excess[i] > EPSILON) || (excess[i] < -EPSILON)) 
			{
				check = 0;
				printf ("c Flow balance constraint violated in node %d. Excess = %.15f\n", 
					i+1,
					excess[i]);
			}
		}
	}

	if (check)
	{
		printf ("c\nc Solution checks as feasible.\n");
	}

	check = 1;
    difference = excess[sink-1] - mincut;
	if ((difference > EPSILON) || (difference < -EPSILON)) 
	{
		check = 0;
		printf("c Flow is not optimal - max flow does not equal min cut! Difference: %.6f, %.6f\nc\n", excess[sink-1], mincut);
	}

	if (check) 
	{
		printf ("c\nc Solution checks as optimal.\nc \n");
		printf ("s Max Flow            : %.3f\n", mincut);
	}

	free (excess);
	excess = NULL;
}


static void
quickSort (Arc **arr, const int first, const int last)
{
	int i, j, left=first, right=last, mid, pivot;
    double x1, x2, x3, pivotval;
	Arc *swap;

	if ((right-left) <= 5)
	{// Bubble sort if 5 elements or less
		for (i=right; (i>left); --i)
		{
			swap = NULL;
			for (j=left; j<i; ++j)
			{
				if (arr[j]->flow < arr[j+1]->flow)
				{
					swap = arr[j];
					arr[j] = arr[j+1];
					arr[j+1] = swap;
				}
			}

			if (!swap)
			{
				return;
			}
		}

		return;
	}

	mid = (first+last)/2;

	x1 = arr[first]->flow; 
	x2 = arr[mid]->flow; 
	x3 = arr[last]->flow;

	pivot = mid;
	
	if (x1 <= x2)
	{
		if (x2 > x3)
		{
			pivot = left;

			if (x1 <= x3)
			{
				pivot = right;
			}
		}
	}
	else
	{
		if (x2 <= x3)
		{
			pivot = right;

			if (x1 <= x3)
			{
				pivot = left;
			}
		}
	}

	pivotval = arr[pivot]->flow;

	swap = arr[first];
	arr[first] = arr[pivot];
	arr[pivot] = swap;

	left = (first+1);

	while (left < right)
	{
		if (arr[left]->flow < pivotval)
		{
			swap = arr[left];
			arr[left] = arr[right];
			arr[right] = swap;
			-- right;
		}
		else
		{
			++ left;
		}
	}

	swap = arr[first];
	arr[first] = arr[left];
	arr[left] = swap;

	if (first < (left-1))
	{
		quickSort (arr, first, (left-1));
	}
	
	if ((left+1) < last)
	{
		quickSort (arr, (left+1), last);
	}
}

static void
sort (Node * current)
{
	if (current->numOutOfTree > 1)
	{
		quickSort (current->outOfTree, 0, (current->numOutOfTree-1));
	}
}

static void
minisort (Node *current) 
{
	Arc *temp = current->outOfTree[current->nextArc];
	int i, size = current->numOutOfTree;
    double tempflow = temp->flow;

	for(i=current->nextArc+1; ((i<size) && (tempflow < current->outOfTree[i]->flow)); ++i)
	{
		current->outOfTree[i-1] = current->outOfTree[i];
	}
	current->outOfTree[i-1] = temp;
}

static void
decompose (Node *excessNode, const int source, int *iteration) 
{
	Node *current = excessNode;
	Arc *tempArc;
	double bottleneck = excessNode->excess;

	for ( ;(current->number != source) && (current->visited < (*iteration)); 
				current = tempArc->from)
	{
		current->visited = (*iteration);
		tempArc = current->outOfTree[current->nextArc];

		if (tempArc->flow < bottleneck)
		{
			bottleneck = tempArc->flow;
		}
	}

	if (current->number == source) 
	{
		excessNode->excess -= bottleneck;
		current = excessNode;

		while (current->number != source) 
		{
			tempArc = current->outOfTree[current->nextArc];
			tempArc->flow -= bottleneck;

			if (tempArc->flow) 
			{
				minisort(current);
			}
			else 
			{
				++ current->nextArc;
			}
			current = tempArc->from;
		}
		return;
	}

	++ (*iteration);

	bottleneck = current->outOfTree[current->nextArc]->flow;

	while (current->visited < (*iteration))
	{
		current->visited = (*iteration);
		tempArc = current->outOfTree[current->nextArc];

		if (tempArc->flow < bottleneck)
		{
			bottleneck = tempArc->flow;
		}
		current = tempArc->from;
	}	
	
	++ (*iteration);

	while (current->visited < (*iteration))
	{
		current->visited = (*iteration);

		tempArc = current->outOfTree[current->nextArc];
		tempArc->flow -= bottleneck;

		if (tempArc->flow) 
		{
			minisort(current);
			current = tempArc->from;
		}
		else 
		{
			++ current->nextArc;
			current = tempArc->from;
		}
	}
}

static void
recoverFlow (void)
{
	int i, j, iteration = 1;
	Arc *tempArc;
	Node *tempNode;

	for (i=0; i<adjacencyList[sink-1].numOutOfTree; ++i) 
	{
		
        tempArc = adjacencyList[sink-1].outOfTree[i];
                
		if (tempArc->from->excess < 0) 
		{
			tempArc->flow -= (double) (-1*tempArc->from->excess);            
			tempArc->from->excess = 0;
		}        
	}
    
	for (i=0; i<adjacencyList[source-1].numOutOfTree; ++i) 
	{
		tempArc = adjacencyList[source-1].outOfTree[i];
		addOutOfTreeNode (tempArc->to, tempArc);
	}
    adjacencyList[source-1].excess = 0;
	adjacencyList[sink-1].excess = 0;
    
	for (i=0; i<numNodes; ++i) 
	{
		tempNode = &adjacencyList[i];

		if ((i == (source-1)) || (i == (sink-1)))
		{
			continue;
		}

		if (tempNode->label >= numNodes) 
		{
			tempNode->nextArc = 0;            
			if ((tempNode->parent) && (tempNode->arcToParent->flow))
            //if ((tempNode->parent) && ((tempNode->arcToParent->flow > EPSILON) || (tempNode->arcToParent->flow < -EPSILON)))
			{
				//printf ("c tempNode->arcToParent->flow: %.10f\n",  tempNode->arcToParent->flow);
                addOutOfTreeNode (tempNode->arcToParent->to, tempNode->arcToParent);
			}
            //printf ("aaaa");
			for (j=0; j<tempNode->numOutOfTree; ++j) 
			{
				if (!tempNode->outOfTree[j]->flow) 
                //if ((tempNode->arcToParent->flow < EPSILON) && (tempNode->arcToParent->flow > -EPSILON)) 
				{
					-- tempNode->numOutOfTree;
					tempNode->outOfTree[j] = tempNode->outOfTree[tempNode->numOutOfTree];
					-- j;
				}
			}

			sort(tempNode);
		}
	}
    for (i=0; i<numNodes; ++i) 
	{
		tempNode = &adjacencyList[i];
		while (tempNode->excess > EPSILON) 
		{
            ++ iteration;
			decompose(tempNode, source, &iteration);
		}
	}    
}



static void
freeMemory (void)
{
	int i;

	for (i=0; i<numNodes; ++i)
	{
		freeRoot (&strongRoots[i]);
	}

	free (strongRoots);

	for (i=0; i<numNodes; ++i)
	{
		if (adjacencyList[i].outOfTree)
		{
			free (adjacencyList[i].outOfTree);
		}
	}

	free (adjacencyList);

	free (labelCount);

	free (arcList);
}

extern void c_readDimacsFileCreateList(const char *filename)
{
    readDimacsFileCreateList(filename);
}
extern void c_simpleInitialization(){simpleInitialization();}
extern void c_reInitialization(const double lambda, const double initial_guess, const int num_edges){reInitialization(lambda, initial_guess, num_edges);}
extern void c_pseudoflowPhase1(){pseudoflowPhase1();}
extern void c_recoverFlow(){recoverFlow();}
extern void c_checkOptimality(){checkOptimality();}
extern double c_getMinCutValue(){
    double cutValue = computeMinCut();
    return cutValue;
}

extern void c_getMinCutEdgeSet(const int num_edges, int *sourceSet)
{
    
    int i = 0;
    
    for (i = 0; i < numNodes; ++i)
    {
        if (adjacencyList[i].number <= num_edges)
        {
            if (adjacencyList[i].label < numNodes)
            {
                sourceSet[adjacencyList[i].number-1] = 1;
            }
            else
            {
                sourceSet[adjacencyList[i].number-1] = 0;
            }
        }
    }   
}

extern void c_getSizeOfMinCutSet(const int num_edges, int *F_edges, int *F_nodes)
{
    
    int i;    
    *F_edges = 0;
    *F_nodes = 0;
    for (i = 0; i < numNodes; ++i)
    {
        
        if (adjacencyList[i].label < numNodes)
        {
            if (adjacencyList[i].number <= num_edges)
            {
                
                (*F_edges)++;
            }
            else
            {
                
                (*F_nodes)++;
            }           
        }
    }
    (*F_nodes)--;
    
}

extern void c_updateCapacities(const double capParameter){updateCapacities(capParameter);}
extern void c_updateSinkCapacities(const double delta){updateSinkCapacities(delta);}
extern void c_updateSrcCapacities(const double newcapacity, const int num_edges) { updateSrcCapacities(newcapacity, num_edges);}

extern void c_reCreateGraph(const int num_edges, const double lambda, const double c)
{ 
    freeMemory();
    highestStrongLabel = 1;
    reCreateGraph(num_edges, lambda, c); 
    simpleInitialization();
}

extern void c_finalfreeMemory()
{
    freeMemory();
    free(initEdgeList);
}

extern void c_freeMemory(){freeMemory();}
extern void c_print_pars()
{
	printf ("c Number of nodes     : %d\n", numNodes);
	printf ("c Number of arcs      : %d\n", numArcs);
#ifdef STATS
	printf ("c Number of arc scans : %lld\n", numArcScans);
	printf ("c Number of mergers   : %d\n", numMergers);
	printf ("c Number of pushes    : %lld\n", numPushes);
	printf ("c Number of relabels  : %d\n", numRelabels);
	printf ("c Number of gaps      : %d\n", numGaps);
#endif
}

int 
main(int argc, char ** argv) 
{    
	printf ("c Pseudoflow algorithm for parametric min cut (version 1.0)\n");
	readDimacsFileCreateList(argv[1]);

#ifdef PROGRESS
	printf ("c Finished reading file.\n"); fflush (stdout);
#endif

	simpleInitialization ();

#ifdef PROGRESS
	printf ("c Finished initialization.\n"); fflush (stdout);
#endif

	pseudoflowPhase1 ();
    

#ifdef PROGRESS
	printf ("c Finished phase 1.\n"); fflush (stdout);
#endif

#ifdef RECOVER_FLOW
	recoverFlow();
	checkOptimality ();
#endif

	printf ("c Number of nodes     : %d\n", numNodes);
	printf ("c Number of arcs      : %d\n", numArcs);
#ifdef STATS
	printf ("c Number of arc scans : %lld\n", numArcScans);
	printf ("c Number of mergers   : %d\n", numMergers);
	printf ("c Number of pushes    : %lld\n", numPushes);
	printf ("c Number of relabels  : %d\n", numRelabels);
	printf ("c Number of gaps      : %d\n", numGaps);
#endif


	freeMemory ();

	return 0;
}
