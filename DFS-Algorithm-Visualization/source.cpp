// Zaharia-Cristea Alina

// I implemented the following algorithms: DFS, Tarjan's algorithm, and Topological Sort on various directed graphs.
// I also implemented the `push` and `pop` functions for a stack, graph initialization and memory deallocation,
// edge addition between two nodes (based on their indices), and console output functions for each algorithm.

// The DFS (Depth-First Search) algorithm takes a graph as input.
// It initially marks all nodes as WHITE (unvisited), and sets a global `time` variable to 0
// to track both discovery and finishing times during traversal.
// Then, it iterates through all nodes; if a node is still WHITE, it calls `dfs_visit` on that node.
// Inside `dfs_visit`, the discovery time is updated, and for each adjacent node,
// if the node is WHITE (unvisited), the function is called recursively.
// Once all adjacent WHITE nodes are visited, the current node is marked BLACK (finished),
// the finish time is recorded, and the parent pointers are updated accordingly.
// The time complexity of DFS is O(V + E).

// Tarjan’s algorithm is used to find strongly connected components in a directed graph.
// The main function `Tarjan` initializes the traversal, and calls `strongConnected` on each node.
// During this traversal, each node is assigned an `index` (its discovery order) and a `lowLink`
// value (the smallest index reachable from the node, including via back-edges).
// Nodes are added to a stack to track the current path.
// At the end, each node gets a `comp` field representing the strongly connected component it belongs to.
// The time complexity of Tarjan's algorithm is O(V + E).

// The Topological Sort algorithm uses DFS to track finish times,
// and uses a stack to store the topological order of the nodes.
// The resulting stack is returned as the final topological sort.
// Note: Topological sorting only works for **DAGs** (Directed Acyclic Graphs),
// since cycles prevent a valid topological ordering.
// The time complexity of the Topological Sort algorithm is O(V + E).


#include<stdio.h>
#include<stdlib.h>
#include"Profiler.h";

Profiler pf("DFS");

typedef enum {
	WHITE = 0,
	GRAY,
	BLACK
}Color;

typedef struct node {
	int key;
	struct node* next;

	int adjSize;
	struct node** adj;
	Color color;
	struct node* parent;
	int discovery_time;
	int finishing_time;

	int index;
	int lowLink;
	boolean onStack;
	int comp;
}Node;

typedef struct {
	int nrNodes;
	Node** v;
}Graph;

void push(Node**stack, Node* node)
{
	if (*stack == NULL)
	{
		*stack = node;
	}
	else
	{
		node->next = *stack;
		*stack = node;
	}
}

Node* pop(Node** stack)
{
	if (*stack == NULL)
	{
		printf("EROARE!\n");
		return NULL;
	}
	else
	{
		Node* node = *stack;
		*stack = (*stack)->next;
		node->next = NULL;
		return node;
	}
}

Graph* initGraph(int nrNodes)
{
	Graph* graph = (Graph*)malloc(sizeof(Graph));

	graph->nrNodes = nrNodes;
	graph->v = (Node**)malloc(nrNodes * sizeof(Node*));

	for (int i = 0; i < nrNodes; i++)
	{
		graph->v[i] = (Node*)malloc(sizeof(Node));
		memset(graph->v[i], 0, sizeof(Node));

		graph->v[i]->key = i;
		graph->v[i]->adjSize = 0;
		graph->v[i]->adj = NULL;
		graph->v[i]->parent = NULL;
		graph->v[i]->discovery_time = -1;
		graph->v[i]->finishing_time = -1;
	}

	return graph;
}

void dfs_visit(Graph* graph, Node* v, int* time, Operation op)
{
	(*time)++;
	v->discovery_time = *time;
	v->color = GRAY;
	op.count(+2);
	
	for (int i = 0; i < v->adjSize; i++)
	{
		op.count();
		if (v->adj[i]->color == WHITE)
		{
			//printf("%d ", v->adj[i]->key);
			op.count();
			v->adj[i]->parent = v;
			dfs_visit(graph, v->adj[i], time, op);
		}
	}
	v->color = BLACK;
	(*time)++;
	v->finishing_time = *time;
	op.count(+2);
}

void dfs(Graph* graph, Operation op)
{
	for (int i = 0; i < graph->nrNodes; i++)
	{
		op.count(+2);
		graph->v[i]->color = WHITE;
		graph->v[i]->parent = NULL;
	}
	int time=0;

	for (int i = 0; i < graph->nrNodes; i++)
	{
		op.count();
		if (graph->v[i]->color == WHITE)
		{
			//printf("%d ", graph->v[i]->key);
			dfs_visit(graph, graph->v[i], &time, op);
		}
	}
}

void freeGraph(Graph* graph)
{
	if (graph->v != NULL) {
		for (int i = 0; i < graph->nrNodes; ++i) {
			if (graph->v[i] != NULL) {
				if (graph->v[i]->adj != NULL) {
					free(graph->v[i]->adj);
					graph->v[i]->adj = NULL;
				}
				graph->v[i]->adjSize = 0;
				free(graph->v[i]);
				graph->v[i] = NULL;
			}
		}
		free(graph->v);
		graph->v = NULL;
	}
	graph->nrNodes = 0;
}

void printGraph(Graph* graph)
{
	if (graph == NULL)
	{
		printf("\tGraful este gol!\n");
	}
	else
	{
		for (int i = 0; i < graph->nrNodes; i++)
		{
			printf("\t\tNode %d: ", graph->v[i]->key);

			for (int j = 0; j < graph->v[i]->adjSize; j++)
			{
				printf("%d ", graph->v[i]->adj[j]->key);
			}
			printf("\n");
		}
	}
}

void recursive_adj_print(Node* node, int lvl)
{

	if (node != NULL && node->color != BLACK)
	{
		printf("\t\t");
		for (int i = 0; i < 4 * lvl; i++)
		{
			printf(" ");
		}
		printf("%d, disc: %d, finish: %d\n\n", node->key, node->discovery_time, node->finishing_time);

		node->color = BLACK;
		for (int i = 0; i < node->adjSize; i++)
		{
			recursive_adj_print(node->adj[i], lvl + 1);
		}
	}
}

void pretty_printDFS(Graph* graph)
{
	if (graph == NULL)
	{
		printf("Graful este gol!\n");
	}
	else
	{
		for (int i = 0; i < graph->nrNodes; i++)
		{
			recursive_adj_print(graph->v[i], 0);
		}
	}
}

void paintAllNodesWhite(Graph* graph)
{
	for (int i = 0; i < graph->nrNodes; i++)
	{
		graph->v[i]->color = WHITE;
	}
}

void recursiveTopological(Node* node, Node** stack)
{
	if (node->color == GRAY) {
		printf("Graful are ciclu!\n");
		return;
	}

	push(stack, node);
	node->color = GRAY;
	for (int i = 0; i < node->adjSize; i++)
	{

		recursiveTopological(node->adj[i], stack);
	}
}

Node** topologicalSort(Graph* graph)
{
	Operation op = pf.createOperation("opDemoTopologicalSort", 12);
	dfs(graph, op);
	Node** stack = (Node**)malloc(graph->nrNodes * sizeof(Node*));
	for (int i = 0; i < graph->nrNodes; i++)
	{
		stack[i] = NULL;
	}

	for (int i = 0; i < graph->nrNodes; i++)
	{
		if (graph->v[i]->color == BLACK)
		{
			recursiveTopological(graph->v[i], stack);
		}
	}
	return stack;
}

void printList(Node** stack)
{
	Node* currentNode = *stack;
	while (currentNode != NULL)
	{
		printf("%d ", currentNode->key);
		currentNode = currentNode->next;
	}
}

void strongConnected(Graph* graph, Node* node, Node**stack, int index, int *nrComponents)
{
	node->index = node->lowLink = index;
	index++;
	push(stack, node);
	node->onStack = true;

	for (int i = 0; i < node->adjSize; i++)
	{
		if (node->adj[i]->index == -1)
		{
			strongConnected(graph, node->adj[i], stack, index, nrComponents);
			node->lowLink = min(node->lowLink, node->adj[i]->lowLink);
		}
		else if (node->adj[i]->onStack)
		{
			node->lowLink = min(node->lowLink, node->adj[i]->index);
		}
	}
	if (node->lowLink == node->index)
	{
		Node* v;
		(*nrComponents)++;
		do {
			v = pop(stack);
			v->onStack = false;
			v->comp = *nrComponents;
		} while (v != node);
	}
}

void tarjan(Graph* graph)
{
	int index = 0;
	int nrComponents = 0;
	Node** stack = (Node**)malloc(sizeof(Node*));

	for (int i = 0; i < graph->nrNodes; i++)
	{
		graph->v[i]->index = graph->v[i]->lowLink - 1;
		graph->v[i]->onStack = false;
		graph->v[i]->comp = 0;
	}
	for (int i = 0; i < graph->nrNodes; i++)
	{
		if (graph->v[i]->index == -1)
			strongConnected(graph, graph->v[i], stack, index, &nrComponents);
	}
}

void addEdge(Graph* graph, int u, int v)
{
	graph->v[u]->adjSize++;
	graph->v[u]->adj = (Node**)realloc(graph->v[u]->adj, (graph->v[u]->adjSize+1) * sizeof(Node*));
	graph->v[u]->adj[graph->v[u]->adjSize - 1] = graph->v[v];
}

void printTarjan(Graph* graph)
{
	for (int i = 0; i < graph->nrNodes; i++)
	{
		printf("\t\tNodul %d -> Componenta %d\n", graph->v[i]->key, graph->v[i]->comp);
	}
}

void demoDFS()
{
	Operation op = pf.createOperation("opDemoDFS", 12);
	Graph* graph;
	int nrNodes = 12;
	graph = initGraph(nrNodes);

	for (int i = 1; i < nrNodes; i++)
	{
		int randomIndex = rand() % i;
		addEdge(graph, randomIndex, i);
	}
	printf("\t* G R A F U L  I N I T I A L *\n\n");
	printGraph(graph);
	printf("\n\n");

	printf("\t\t* D F S *\n\n");
	dfs(graph, op);
	paintAllNodesWhite(graph);
	pretty_printDFS(graph);
	freeGraph(graph);
	printf("\n");
}

void demoSortareTopologica()
{
	Graph* graph;
	int nrNodes = 9;
	graph = initGraph(nrNodes);

	for (int i = 1; i < nrNodes; i++)
	{
		int randomIndex = rand() % i;
		addEdge(graph, randomIndex, i);
		/*graph->v[randomIndex]->adj = (Node**)realloc(graph->v[randomIndex]->adj, (graph->v[randomIndex]->adjSize + 1) * sizeof(Node*));
		graph->v[randomIndex]->adj[graph->v[randomIndex]->adjSize] = graph->v[i];
		(graph->v[randomIndex]->adjSize)++;*/
	}
	printf("\t* G R A F U L  I N I T I A L *\n\n");
	printGraph(graph);
	printf("\n\n\n");

	printf("\t* S O R T A R E A  T O P O L O G I C A *\n\n\t\t");
	Node**stack = topologicalSort(graph);
	printList(stack);
	printf("\n");
	free(stack);
	freeGraph(graph);
}

void demoTarjan()
{
	Graph* graph;
	int nrNodes = 8;
	graph = initGraph(nrNodes);

	addEdge(graph, 0, 2);
	addEdge(graph, 1, 0);
	addEdge(graph, 1, 3);
	addEdge(graph, 2, 1);
	addEdge(graph, 2, 3);
	addEdge(graph, 2, 4);
	addEdge(graph, 3, 5);
	addEdge(graph, 4, 5);
	addEdge(graph, 4, 6);
	addEdge(graph, 5, 3);
	addEdge(graph, 5, 7);
	addEdge(graph, 6, 4);
	addEdge(graph, 6, 7);

	printf("\t* G R A F U L  I N I T I A L *\n\n");
	printGraph(graph);
	printf("\n\n");

	printf("\t* COOMPONENTELE TARE CONEXE (TARJAN) *\n\n");
	tarjan(graph);
	printTarjan(graph);
	printf("\n");
	freeGraph(graph);
}

int edgeAlreadyExists(Graph* graph, int u, int v)
{
	for (int i = 0; i < graph->v[u]->adjSize; i++)
	{
		if (graph->v[u]->adj[i]->key == v)
			return 1;
	}
	return 0;
}

void performance()
{
	int n, i;

	for (n = 1000; n <= 4500; n += 100)
	{
		Operation op1 = pf.createOperation("DFS-edges", n);
		Graph* graph;
		int nrNodes = 100;
		graph = initGraph(nrNodes);

		for (int i = 0; i < n; i++)
		{
			int randomIndexU, randomIndexV;

			do {
				randomIndexU = rand() % nrNodes;
				randomIndexV = rand() % nrNodes;
			} while (randomIndexU == randomIndexV || edgeAlreadyExists(graph, randomIndexU, randomIndexV));

			addEdge(graph, randomIndexU, randomIndexV);
		}
		dfs(graph, op1);
		freeGraph(graph);
	}

	for (n = 100; n <= 200; n += 10)
	{	
		Operation op2 = pf.createOperation("DFS-vertices", n);
		Graph* graph;
		graph = initGraph(n);

		int edges = 4500;
		for (int i = 0; i < edges; i++)
		{
			int randomIndexU, randomIndexV;

			do {
				randomIndexU = rand() % n;
				randomIndexV = rand() % n;
			} while (randomIndexU == randomIndexV || edgeAlreadyExists(graph, randomIndexU, randomIndexV));

			addEdge(graph, randomIndexU, randomIndexV);
		}
		dfs(graph, op2);
		freeGraph(graph);
	}

	pf.showReport();
}


int main()
{
	demoDFS();
	demoSortareTopologica();
	demoTarjan();
	performance();
	return 0;
}