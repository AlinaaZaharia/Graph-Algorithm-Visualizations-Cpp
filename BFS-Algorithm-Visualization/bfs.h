#ifndef __BFS_H__
#define __BFS_H__

#include "Profiler.h"

#define MAX_ROWS 100
#define MAX_COLS 100

typedef struct{
    int rows;
    int cols;
    int mat[MAX_ROWS][MAX_COLS];
}Grid;

typedef struct{
    int row;
    int col;
}Point;

typedef enum{
    WHITE = 0,
    GRAY,
    BLACK
} Color;

typedef struct _Node{
    Point position;
    int adjSize;
    struct _Node **adj;

    Color color;
    int dist;
    struct _Node *parent;
}Node;

typedef struct Edge {
    Node* u;
    Node* v;
}Edge;

typedef struct{
    int nrNodes;
    Node **v;
    Edge* edges;
}Graph;

int get_neighbors(const Grid *grid, Point p, Point neighb[]);
void grid_to_graph(const Grid *grid, Graph *graph);
void free_graph(Graph *graph);
void bfs(Graph *graph, Node *s, Operation *op=NULL);
void print_bfs_tree(Graph *graph);
int shortest_path(Graph *graph, Node *start, Node *end, Node *path[]);
void performance();

#endif