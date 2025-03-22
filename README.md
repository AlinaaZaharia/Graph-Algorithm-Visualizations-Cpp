# Graph Algorithm Visualizations in C++

A collection of C++ console applications demonstrating **graph traversal and analysis algorithms**, including:

- Breadth-First Search (BFS) on 2D grids with real-time visualization
- Depth-First Search (DFS), **Tarjan’s algorithm** (strongly connected components), and **Topological Sort** on directed graphs

Both projects include performance profiling and complexity analysis using custom benchmarks and graphical reports.


## 📁 Project Structure

### `bfs-grid/` – BFS Pathfinding in a 2D Grid

A console app that reads a customizable `grid.txt` file (`0 = free`, `1 = wall`) and allows BFS-based traversal and shortest path reconstruction.

#### 🔧 Implemented Functionality

- `get_neighbors()` – Finds up/down/left/right neighbors within grid bounds
- `bfs()` – Classic BFS with color marking, distance, and parent tracking
- `pretty_print()` – Tree-style display of BFS parent relationships
- `shortest_path()` – Computes the shortest path from source to destination
- `performance()` – Evaluates BFS complexity (O(V + E)) with profiler
- 🔲 Terminal grid display with direction arrows and color codes

#### 📊 Complexity Reports

Performance measured in terms of vertex/edge operations.

### `graph-algos/` – DFS, Tarjan, Topological Sort on Directed Graphs

A C++ app that generates and analyzes directed graphs using classic and advanced algorithms.

#### 🔧 Implemented Algorithms

- `DFS()` – Depth-First Search with discovery/finishing time and parent tracking
- `Tarjan()` – Strongly Connected Components detection (O(V + E)) using low-link values and a stack
- `TopologicalSort()` – Ordering nodes using finishing time from DFS
- `pretty_printDFS()` – Tree-style indentation of DFS traversal
- `performance()` – Benchmarks DFS on varying node/edge counts

#### 📊 Complexity Reports

Generated using custom profiler.


## 🚀 How to Run

1. Open project folder (`bfs-grid/` or `graph-algos/`) in your C++ IDE (e.g., Code::Blocks, CLion, VS Code).
2. Modify the `main()` to activate the desired demo (BFS traversal, pathfinding, Tarjan, etc.).
3. Build & run.
4. For grid-based traversal, edit `grid.txt` to define custom maps.


## 🛠️ Notes

Both projects were built on top of a provided C++ template which included base structures and utilities. All core algorithmic logic was implemented separately, including:

- BFS traversal
- Neighbor detection
- Shortest path finding
- DFS with time/parent management
- Tarjan’s SCC algorithm
- Topological sorting
