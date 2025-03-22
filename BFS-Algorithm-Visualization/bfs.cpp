//Zaharia Cristea Alina

// I implemented the 'get_neighb' function, which receives the position of a cell in the grid and colors its neighbors 
// (up, down, left, and right) if they exist and are within the grid boundaries.
//
// I also implemented the BFS algorithm, which starts traversal from a given start node.
// The traversal works as follows: initially, all nodes are colored WHITE, the 'dist' field (for distance) is 
// initialized to 0, and all node parents are set to NULL. During the BFS traversal, these fields are updated: 
// visited nodes are recolored to BLACK, and their 'dist' field is updated based on the distance from the start node, 
// with parents set accordingly.
//
// The complexity of the algorithm, as shown when calling the performance function (which generates two graphs: 
// one for vertex operations and another for edge operations), is: O(v + e).
//
// I also implemented a 'pretty_print' function for BFS, similar to the one used in lab 6 for displaying the k-ary tree 
// based on the parent vector.
//
// Additionally, I implemented the 'shortest_path' function, which determines the shortest path from a given start node 
// to a final 'end' node. If the 'end' node cannot be reached from the 'start' node, the function returns -1.
// This implementation uses the previously built BFS algorithm, but also computes the path length by 
// starting from the 'end' node and traversing backward through the parent nodes using a temporary node 'current'.


#include <stdlib.h>
#include <string.h>
#include "bfs.h"

#define CAPACITY 100

typedef struct _queue {
    Node* vec[CAPACITY];
    int size;
    int head;
    int tail;
} Queue;

void initQueue(Queue* my_queue)
{
    my_queue->size = 0;
    my_queue->head = 0;
    my_queue->tail = 0;
}

void enqueue(Queue* my_queue, Node* s)
{
    if (my_queue->size == CAPACITY - 1)
    {
        printf("Overflow\n");
        return;
    }

    my_queue->vec[my_queue->tail] = s;
    (my_queue->size)++;

    if (my_queue->tail == CAPACITY - 1)
        my_queue->tail = 0;
    else
        (my_queue->tail)++;
}

Node* dequeue(Queue* my_queue)
{
    if (my_queue->size == 0)
    {
        //printf("Underflow\n");
        return NULL;
    }

    Node* key = my_queue->vec[my_queue->head];
    if (my_queue->head == CAPACITY - 1)
        my_queue->head = 0;
    else
        (my_queue->head)++;
    (my_queue->size)--;
    return key;
}



int get_neighbors(const Grid* grid, Point p, Point neighb[])
{
    int counter = 0;
    if (grid->mat[p.row][p.col] == 0 && p.row > 0 && p.row < grid->rows && p.col > 0 && p.col < grid->cols)
    {
        if (grid->mat[p.row - 1][p.col] == 0)     // verific daca exista vecin sus
        {
            neighb[counter].row = p.row - 1;
            neighb[counter].col = p.col;
            counter++;
        }

        if (grid->mat[p.row + 1][p.col] == 0)   // verific daca exista vecin jos
        {
            neighb[counter].row = p.row + 1;
            neighb[counter].col = p.col;
            counter++;
        }

        if (grid->mat[p.row][p.col - 1] == 0)   // vecin stanga
        {
            neighb[counter].row = p.row;
            neighb[counter].col = p.col - 1;
            counter++;
        }

        if (grid->mat[p.row][p.col + 1] == 0)   // vecin dreapta
        {
            neighb[counter].row = p.row;
            neighb[counter].col = p.col + 1;
            counter++;
        }
    }

    // TODO: fill the array neighb with the neighbors of the point p and return the number of neighbors
    // the point p will have at most 4 neighbors (up, down, left, right)
    // avoid the neighbors that are outside the grid limits or fall into a wall
    // note: the size of the array neighb is guaranteed to be at least 4
    return counter;
}

void grid_to_graph(const Grid* grid, Graph* graph)
{
    //we need to keep the nodes in a matrix, so we can easily refer to a position in the grid
    Node* nodes[MAX_ROWS][MAX_COLS];
    int i, j, k;
    Point neighb[4];

    //compute how many nodes we have and allocate each node
    graph->nrNodes = 0;
    for (i = 0; i < grid->rows; ++i) {
        for (j = 0; j < grid->cols; ++j) {
            if (grid->mat[i][j] == 0) {
                nodes[i][j] = (Node*)malloc(sizeof(Node));
                memset(nodes[i][j], 0, sizeof(Node)); //initialize all fields with 0/NULL
                nodes[i][j]->position.row = i;
                nodes[i][j]->position.col = j;
                ++graph->nrNodes;
            }
            else {
                nodes[i][j] = NULL;
            }
        }
    }
    graph->v = (Node**)malloc(graph->nrNodes * sizeof(Node*));
    k = 0;
    for (i = 0; i < grid->rows; ++i) {
        for (j = 0; j < grid->cols; ++j) {
            if (nodes[i][j] != NULL) {
                graph->v[k++] = nodes[i][j];
            }
        }
    }

    //compute the adjacency list for each node
    for (i = 0; i < graph->nrNodes; ++i) {
        graph->v[i]->adjSize = get_neighbors(grid, graph->v[i]->position, neighb);
        if (graph->v[i]->adjSize != 0) {
            graph->v[i]->adj = (Node**)malloc(graph->v[i]->adjSize * sizeof(Node*));
            k = 0;
            for (j = 0; j < graph->v[i]->adjSize; ++j) {
                if (neighb[j].row >= 0 && neighb[j].row < grid->rows &&
                    neighb[j].col >= 0 && neighb[j].col < grid->cols &&
                    grid->mat[neighb[j].row][neighb[j].col] == 0) {
                    graph->v[i]->adj[k++] = nodes[neighb[j].row][neighb[j].col];
                }
            }
            if (k < graph->v[i]->adjSize) {
                //get_neighbors returned some invalid neighbors
                graph->v[i]->adjSize = k;
                graph->v[i]->adj = (Node**)realloc(graph->v[i]->adj, k * sizeof(Node*));
            }
        }
    }
}

void free_graph(Graph* graph)
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

void bfs(Graph* graph, Node* s, Operation* op)
{
    if (op != NULL)
    {
        op->count();
    }

    for (int i = 0; i < graph->nrNodes; i++)
    {
        if (op != NULL)
        {
            op->count(+3);
        }

        graph->v[i]->color = WHITE;
        graph->v[i]->dist = 0;
        graph->v[i]->parent = NULL;
    }

    s->color = GRAY;
    s->dist = 0;
    s->parent = NULL;
    if (op != NULL)
    {
        op->count(+3);
    }

    Queue my_queue;
    initQueue(&my_queue);
    enqueue(&my_queue, s);

    while (my_queue.size != 0)
    {
        Node* currentNodeU = dequeue(&my_queue);

        if (op != NULL)
        {
            op->count(+2);
        }
        
        for (int i = 0; i < currentNodeU->adjSize; i++)
        {
            Node* currentNodeV = currentNodeU->adj[i];
            if (op != NULL)
            {
                op->count();
            }

            if (currentNodeV->color == WHITE)
            {
                currentNodeV->color = GRAY;
                currentNodeV->dist = currentNodeU->dist + 1;
                currentNodeV->parent = currentNodeU;
                enqueue(&my_queue, currentNodeV);
                if (op != NULL)
                {
                    op->count(+3);
                }
            }
        }
        currentNodeU->color = BLACK;
        if (op != NULL)
        {
            op->count();
        }
    }
    // TOOD: implement the BFS algorithm on the graph, starting from the node s
    // at the end of the algorithm, every node reachable from s should have the color BLACK
    // for all the visited nodes, the minimum distance from s (dist) and the parent in the BFS tree should be set
    // for counting the number of operations, the optional op parameter is received
    // since op can be NULL (when we are calling the bfs for display purposes), you should check it before counting:
    // if(op != NULL) op->count();
}


void recursive_child_print(int p[], Point repr[], int length, int parent, int depth)
{
    for (int i = 0; i < length; i++)
    {
        if (p[i] == parent)
        {
            for (int j = 0; j < 8 * depth; j++)
            {
                printf(" ");
            }
            printf("(%d, %d)\n", repr[i].row, repr[i].col);

            recursive_child_print(p, repr, length, i, depth + 1);
        }
    }
}

void pretty_print(int p[], Point repr[], int length)
{
    for (int i = 0; i < length; i++)
    {
        if (p[i] == -1)
        {
            printf("(%d, %d)\n", repr[i].row, repr[i].col);
            recursive_child_print(p, repr, length, i, 1);
        }
    }
}


void print_bfs_tree(Graph* graph)
{
    //first, we will represent the BFS tree as a parent array
    int n = 0; //the number of nodes
    int* p = NULL; //the parent array
    Point* repr = NULL; //the representation for each element in p

    //some of the nodes in graph->v may not have been reached by BFS
    //p and repr will contain only the reachable nodes
    int* transf = (int*)malloc(graph->nrNodes * sizeof(int));
    for (int i = 0; i < graph->nrNodes; ++i) {
        if (graph->v[i]->color == BLACK) {
            transf[i] = n;
            ++n;
        }
        else {
            transf[i] = -1;
        }
    }
    if (n == 0) {
        //no BFS tree
        free(transf);
        return;
    }

    int err = 0;
    p = (int*)malloc(n * sizeof(int));
    repr = (Point*)malloc(n * sizeof(Node));
    for (int i = 0; i < graph->nrNodes && !err; ++i) {
        if (graph->v[i]->color == BLACK) {
            if (transf[i] < 0 || transf[i] >= n) {
                err = 1;
            }
            else {
                repr[transf[i]] = graph->v[i]->position;
                if (graph->v[i]->parent == NULL) {
                    p[transf[i]] = -1;
                }
                else {
                    err = 1;
                    for (int j = 0; j < graph->nrNodes; ++j) {
                        if (graph->v[i]->parent == graph->v[j]) {
                            if (transf[j] >= 0 && transf[j] < n) {
                                p[transf[i]] = transf[j];
                                err = 0;
                            }
                            break;
                        }
                    }
                }
            }
        }
    }
    free(transf);
    transf = NULL;

    if (!err) {
        pretty_print(p, repr, n);

        // TODO: pretty print the BFS tree
        // the parrent array is p (p[k] is the parent for node k or -1 if k is the root)
        // when printing the node k, print repr[k] (it contains the row and column for that point)
        // you can adapt the code for transforming and printing multi-way trees from the previous labs
    }

    if (p != NULL) {
        free(p);
        p = NULL;
    }
    if (repr != NULL) {
        free(repr);
        repr = NULL;
    }
}

int shortest_path(Graph* graph, Node* start, Node* end, Node* path[])
{
    if (graph == NULL || start == NULL || end == NULL)
    {
        return -1;
    }

    for (int i = 0; i < graph->nrNodes; i++)
    {
        graph->v[i]->color = WHITE;
        graph->v[i]->dist = 0;
        graph->v[i]->parent = NULL;
    }
    start->color = GRAY;
    start->dist = 0;
    start->parent = NULL;

    Queue my_queue;
    initQueue(&my_queue);
    enqueue(&my_queue, start);

    while (my_queue.size != 0)
    {
        Node* currentNodeU = dequeue(&my_queue);

        for (int i = 0; i < currentNodeU->adjSize; i++)
        {
            Node* currentNodeV = currentNodeU->adj[i];

            if (currentNodeV->color == WHITE)
            {
                currentNodeV->color = GRAY;
                currentNodeV->dist = currentNodeU->dist + 1;
                currentNodeV->parent = currentNodeU;
                enqueue(&my_queue, currentNodeV);
            }
        }
        currentNodeU->color = BLACK;

        if (currentNodeU == end)
        {
            int length = 0;
            Node* current = end;

            while (current != NULL)
            {
                path[length] = current;
                current = current->parent;
                (length)++;
            }
            return length;
        }
    }
    return -1;

    // TODO: 
    // compute the shortest path between the nodes start and end in the given graph
    // the nodes from the path, should be filled, in order, in the array path
    // the number of nodes filled in the path array should be returned
    // if end is not reachable from start, return -1
    // note: the size of the array path is guaranteed to be at least 1000
}


void performance()
{
    int n, i;
    Profiler p("bfs");

    // vary the number of edges
    for (n = 1000; n <= 4500; n += 100) 
    {
        Operation op1 = p.createOperation("bfs-edges", n);
        Graph graph;
        graph.nrNodes = 100;
        //initialize the nodes of the graph
        graph.v = (Node**)malloc(graph.nrNodes * sizeof(Node*));
        for (i = 0; i < graph.nrNodes; ++i) 
        {
            graph.v[i] = (Node*)malloc(sizeof(Node));
            memset(graph.v[i], 0, sizeof(Node));
        }
        graph.edges = (Edge*)malloc(n * sizeof(Edge));
        // TODO: generate n random edges
        // make sure the generated graph is connected
        int randomIndex = rand() % n;
        Node* currentNode = graph.v[randomIndex];
        op1.count();

        for (int i = n - 1; i > 0; i--)
        {
            randomIndex = rand() % i;
            graph.edges[randomIndex].u = currentNode;
            graph.edges[randomIndex].v = graph.v[randomIndex];
            currentNode = graph.v[randomIndex];
            op1.count(+3);
        }

        bfs(&graph, graph.v[0], &op1);
        free_graph(&graph);
    }

    // vary the number of vertices
    for (n = 100; n <= 200; n += 10) 
    {
        Operation op2 = p.createOperation("bfs-vertices", n);
        Graph graph;
        graph.nrNodes = n;
        //initialize the nodes of the graph
        graph.v = (Node**)malloc(graph.nrNodes * sizeof(Node*));
        for (i = 0; i < graph.nrNodes; ++i) 
        {
            graph.v[i] = (Node*)malloc(sizeof(Node));
            memset(graph.v[i], 0, sizeof(Node));
        }

        int randomIndex = rand() % n;
        Node* currentNode = graph.v[randomIndex];
        op2.count();

        for (int i = n - 1; i > 0; i--)
        {
            randomIndex = rand() % i;
            graph.edges[randomIndex].u = currentNode;
            graph.edges[randomIndex].v = graph.v[randomIndex];
            currentNode = graph.v[randomIndex];
            op2.count(+3);
        }

        for (int i = 4500 - n; i > 0; i-=2)
        {
            randomIndex = rand() % i;
            currentNode = graph.v[randomIndex];
            int randomIndex1 = rand() % i;
            graph.edges[randomIndex] = { currentNode, graph.v[randomIndex1] };
            op2.count(+3);
        }
        // TODO: generate 4500 random edges
        // make sure the generated graph is connected

        bfs(&graph, graph.v[0], &op2);
        free_graph(&graph);
    }
    p.showReport();
}
