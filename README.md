# MineField

This is a C program for analyzing a graph representing a minefield. The program allows you to perform various operations on the graph, including adding and removing nodes (mines), adding and removing edges (connections between mines), performing depth-first searches (DFS) on the graph, calculating the effectiveness of a given mine, finding the most effective mines, and calculating the area covered by explosions from the most effective mine.

> My primary goals in developing this minefield project were to create an efficient graph representation for mines, implement a flexible DFS traversal, calculate and display mine effectiveness, identify the most effective mines, calculate area coverage by explosions, and enable graph visualization for enhanced minefield analysis.
> 
> This project was developed as the third university assignment for "Algorithms and Data Structures", University of Belgrade School of Electrical Engineering majoring Software Engineering. Please refer to the [instructions.pdf](instructions.pdf) file for detailed assignment instructions.


## Table of Contents

- [Features](#features)
- [Getting Started](#getting-started)
  - [Prerequisites](#prerequisites)
  - [Usage](#usage)
- [Menu Options](#menu-options)
- [File Format](#file-format)


## Features
- Create and manipulate a graph representing a minefield.
- Perform DFS traversal starting from any mine.
- Calculate the effectiveness of a given mine in terms of its explosion radius.
- Find and print the most effective mines in the graph.
- Calculate the area covered by explosions from the most effective mine.
- Visualize the minefield graph with coordinates in an output file.

## Getting Started

### Prerequisites

To compile and run the MineField simulation, you need a C compiler (e.g. GCC) installed on your system.

### Usage

1. Clone this repository to your local machine.
2. Open a terminal or command prompt.
3. Navigate to the project directory.
4. Compile the code using a C compiler (e.g., GCC)
5. Follow the on-screen instructions to create nodes, establish links, perform other operations on the graph and much more.

### Menu Options

1. **Dodaj cvor u graf**: Add a node (mine) to the graph.
2. **Ukloni cvor iz grafa**: Remove a node (mine) from the graph.
3. **Dodaj granu u graf**: Add an edge (connection) between two nodes (mines) in the graph.
4. **Ukloni granu iz grafa**: Remove an edge (connection) between two nodes (mines) in the graph.
5. **DFS ispis grafa**: Perform a depth-first search (DFS) traversal of the graph starting from a selected mine.
6. **Baci raketu**: Simulate launching a rocket at a specified mine location.
7. **Efikasnost zadate mine**: Calculate and print the effectiveness of a specified mine.
8. **Pronadji minu max efikasnosti**: Find and print the most effective mines in the graph.
9. **Povrsina mine sa max efikasnoscu**: Calculate and print the area covered by explosions from the most effective mine.
10. **Ispisi graf sa koordinatama**: Visualize the minefield graph with coordinates and save it to a file.
11. **Izlaz**: Exit the program.

## File format
When importing graph information from a file, the format should be as follows:
```
<number_of_nodes>
<x1> <y1> <radius1>
<x2> <y2> <radius2>
...
```
This project includes 3 examples, filenames `graf1.txt`, `graf2.txt` and `grafTest.txt` (this test was made to cover all edge cases of this program).
