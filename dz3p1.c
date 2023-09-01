#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#define MAX_NODES 1000
//<editor-fold desc="Typedefovi">
typedef struct {
    double x, y, r;
} Mina;

typedef struct Edge {
    Mina source, destination;
} Edge;

typedef struct Graph {
    int numNodes, numEdges;
    Edge *edges;
    Mina *nodes;
} Graph;

typedef struct Stack {
    int top;
    unsigned capacity;
    Mina *array;
} Stack;

typedef struct VisitedMina {
    double x, y, r;
    struct VisitedMina *next;
} VisitedMina;
//</editor-fold>

//<editor-fold desc="STEK funkcije">

Stack *createStack(int capacity) {
    Stack *stack = (Stack *) malloc(sizeof(Stack));
    stack->capacity = capacity;
    stack->top = -1;
    stack->array = (Mina *) malloc(stack->capacity * sizeof(Mina));
    return stack;
}

int isFull(Stack *stack) {
    return stack->top == stack->capacity - 1;
}

int isEmpty(Stack *stack) {
    return stack->top == -1;
}

void push(Stack *stack, Mina item) {
    if (isFull(stack))
        return;
    stack->array[++stack->top] = item;
}

Mina pop(Stack *stack) {
    if (isEmpty(stack)) {
        Mina empty_mina = {0, 0, 0};
        return empty_mina;
    }
    return stack->array[stack->top--];
}

int is_in_stack(Stack *stack, Mina mina) {
    for (int i = 0; i <= stack->top; i++) {
        if (stack->array[i].x == mina.x && stack->array[i].y == mina.y && stack->array[i].r == mina.r) {
            return 1;
        }
    }
    return 0;
}

int is_adjacent(Graph *graph, Mina current_mina, Mina mina) {
    for (int i = 0; i < graph->numEdges; i++) {
        if (graph->edges[i].source.x == current_mina.x && graph->edges[i].source.y == current_mina.y &&
            graph->edges[i].source.r == current_mina.r
            && graph->edges[i].destination.x == mina.x && graph->edges[i].destination.y == mina.y &&
            graph->edges[i].destination.r == mina.r) {
            return 1;
        }
    }
    return 0;
}

//</editor-fold>

//<editor-fold desc="Funkcije za rad sa grafom">
Edge createEdge(Mina source, Mina destination) {
    Edge newEdge;
    newEdge.source = source;
    newEdge.destination = destination;
    return newEdge;
}

Graph *createGraph(int numNodes) {
    Graph *graph = (Graph *) malloc(sizeof(Graph));
    graph->numNodes = numNodes;
    graph->numEdges = 0;
    graph->edges = NULL;
    graph->nodes = (Mina *) calloc(numNodes, sizeof(Mina));
    return graph;
}

//void addNode(Graph *graph, Mina mina, int index) { // stari addNode
////    if (index >= 0 && index < graph->numNodes) {
//    graph->nodes[index] = mina;
////    } else {
////        printf("Node cannot be added. It is outside the range.\n");
////    }
//}

void addNode(Graph *graph, Mina mina, int index) {
    if (index >= 0 && index <= graph->numNodes) {
        graph->nodes = realloc(graph->nodes, (graph->numNodes + 1) * sizeof(Mina));
        if (graph->nodes == NULL) {
            printf("Neuspesna alokacija memorije za novi Cvor\n");
            return;
        }
        for (int i = graph->numNodes; i > index; i--) {
            graph->nodes[i] = graph->nodes[i - 1];
        }
        graph->nodes[index] = mina;
    } else {
        printf("Nevalidan index. Cvor ne moze biti dodat.\n");
    }
}


void removeNode(Graph *graph, Mina mina) {
    int nodeIndex = -1;
    for (int i = 0; i < graph->numNodes; i++) {
        if (graph->nodes[i].x == mina.x && graph->nodes[i].y == mina.y && graph->nodes[i].r == mina.r) {
            nodeIndex = i;
            break;
        }
    }

    if (nodeIndex != -1) {
        for (int i = nodeIndex; i < graph->numNodes - 1; i++) {
            graph->nodes[i] = graph->nodes[i + 1];
        }
        graph->numNodes--;

        for (int i = 0; i < graph->numEdges; i++) {
            if (graph->edges[i].source.x == mina.x && graph->edges[i].source.y == mina.y &&
                graph->edges[i].source.r == mina.r
                || graph->edges[i].destination.x == mina.x && graph->edges[i].destination.y == mina.y &&
                   graph->edges[i].destination.r == mina.r) {

                for (int j = i; j < graph->numEdges - 1; j++) {
                    graph->edges[j] = graph->edges[j + 1];
                }
                graph->numEdges--;
                i--;
            }
        }
        printf("Uspesno uklonjen cvor!\n");
    } else {
        printf("Mina sa zadatim koordinatama ne postoji u grafu.\n");

    }
}

int edgeExists(Graph *graph, Mina source, Mina destination) {
    for (int i = 0; i < graph->numEdges; i++) {
        if (graph->edges[i].source.x == source.x && graph->edges[i].source.y == source.y &&
            graph->edges[i].source.r == source.r
            && graph->edges[i].destination.x == destination.x && graph->edges[i].destination.y == destination.y &&
            graph->edges[i].destination.r == destination.r) {
            return 1;
        }
    }
    return 0;
}

void addEdge(Graph *graph, Mina source, Mina destination) {
    if (!edgeExists(graph, source, destination)) {
        if (graph->numEdges % MAX_NODES == 0) {  // Whenever current numEdges is a multiple of MAX_NODES, realloc
            Edge *tempEdges = realloc(graph->edges, (graph->numEdges + MAX_NODES) * sizeof(Edge));
            if (tempEdges == NULL) {
                printf("Dodeljivanje memorije za novu granu nije uspelo.\n");
                return;
            }
            graph->edges = tempEdges;
        }
        graph->edges[graph->numEdges++] = createEdge(source, destination);
        printf("Uspesno dodata grana!\n");
    } else {
        printf("Grana vec postoji u grafu.\n");
    }
}

void removeEdge(Graph *graph, Mina source, Mina destination) {
    if (edgeExists(graph, source, destination)) {
        for (int i = 0; i < graph->numEdges; i++) {
            if (graph->edges[i].source.x == source.x && graph->edges[i].source.y == source.y &&
                graph->edges[i].source.r == source.r
                && graph->edges[i].destination.x == destination.x && graph->edges[i].destination.y == destination.y &&
                graph->edges[i].destination.r == destination.r) {

                for (int j = i; j < graph->numEdges - 1; j++) {
                    graph->edges[j] = graph->edges[j + 1];
                }
                graph->numEdges--;
                break;
            }
        }
    } else {
        printf("Grana ne postoji.\n");
    }
}

void printGraph(Graph *graph) {
    printf("Broj cvorova: %d\n", graph->numNodes);
    printf("Broj grana: %d\n", graph->numEdges);

    printf("\nCvorovi:\n");
    for (int i = 0; i < graph->numNodes; i++) {
        printf("Cvor %d: x = %.2lf, y = %.2lf, r = %.2lf\n", i, graph->nodes[i].x, graph->nodes[i].y,
               graph->nodes[i].r);
    }

    printf("\nGrane:\n");
    for (int i = 0; i < graph->numEdges; i++) {
        int sourceIndex = -1, destIndex = -1;

        for (int j = 0; j < graph->numNodes; j++) {
            if (graph->nodes[j].x == graph->edges[i].source.x &&
                graph->nodes[j].y == graph->edges[i].source.y &&
                graph->nodes[j].r == graph->edges[i].source.r) {
                sourceIndex = j;
            }
            if (graph->nodes[j].x == graph->edges[i].destination.x &&
                graph->nodes[j].y == graph->edges[i].destination.y &&
                graph->nodes[j].r == graph->edges[i].destination.r) {
                destIndex = j;
            }
        }

        printf("Grana %d: source (Cvor %d) = (%.2lf, %.2lf, %.2lf) \xC4> destination (Cvor %d) = (%.2lf, %.2lf, %.2lf)\n",
               i,
               sourceIndex,
               graph->edges[i].source.x, graph->edges[i].source.y, graph->edges[i].source.r,
               destIndex,
               graph->edges[i].destination.x, graph->edges[i].destination.y, graph->edges[i].destination.r);
    }
}


bool ubaciUGrafSaProverom(Graph *graph, double x, double y, double r, int index) {
    for (int i = 0; i < graph->numNodes; i++) {
        if (graph->nodes[i].x == x && graph->nodes[i].y == y && graph->nodes[i].r == r) {
            printf("Mina vec postoji na zadatim koordinatama.\n");
            return false;
        }
    }

    Mina mina = {x, y, r};
    printf("Ubacujem minu: %lf %lf %lf\n", x, y, r);
    addNode(graph, mina, index);
    printf("Uspesno ubacena mina!\n");
    return true;
}

void ubaciUGraf(Graph *graph, double x, double y, double r, int index) {
    for (int i = 0; i < graph->numNodes; i++) {
        if (graph->nodes[i].x == x && graph->nodes[i].y == y && graph->nodes[i].r == r) {
            printf("Mina vec postoji na zadatim koordinatama.\n");
            return;
        }
    }

    Mina mina = {x, y, r};
    printf("Ubacujem minu: %lf %lf %lf\n", x, y, r);
    addNode(graph, mina, index);
    printf("Uspesno ubacena mina!\n");
}
//</editor-fold>

//<editor-fold desc="Rad sa minama">
double distance(Mina a, Mina b) {
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    return sqrt(dx * dx + dy * dy);
}

int getMinaIndex(Graph *graph, Mina mina) {
    for (int i = 0; i < graph->numNodes; i++) {
        if (graph->nodes[i].x == mina.x && graph->nodes[i].y == mina.y && graph->nodes[i].r == mina.r) {
            return i;
        }
    }
    return -1;
}

//rekurzivni dfs, radi ali ne treba
//<editor-fold desc="DFS rekruzivni">
//void dfs(Graph* graph, Mina source, int* visited) {
//    for(int i = 0; i < graph->numEdges; i++) {
//        if(graph->edges[i].source.x == source.x && graph->edges[i].source.y == source.y && graph->edges[i].source.r == source.r) {
//            Mina destination = graph->edges[i].destination;
//            int destination_index = -1;
//            for (int j = 0; j < graph->numNodes; j++) {
//                if (graph->nodes[j].x == destination.x && graph->nodes[j].y == destination.y && graph->nodes[j].r == destination.r) {
//                    destination_index = j;
//                    break;
//                }
//            }
//            if (destination_index != -1 && !visited[destination_index]) {
//                visited[destination_index] = 1;
//                dfs(graph, destination, visited);
//            }
//        }
//    }
//}
//</editor-fold>

void dfs(Graph *graph, Mina source, int *visited) {
    Stack *stack = createStack(graph->numNodes);

    int source_index = getMinaIndex(graph, source);
    if (source_index == -1) {
        return;
    }

    push(stack, source);

    while (!isEmpty(stack)) {
        Mina mina = pop(stack);
        int mina_index = getMinaIndex(graph, mina);

        if (!visited[mina_index]) {
            visited[mina_index] = 1;

            for (int i = 0; i < graph->numNodes; i++) {
                if (!visited[i] && is_adjacent(graph, mina, graph->nodes[i]) && !is_in_stack(stack, graph->nodes[i])) {
                    push(stack, graph->nodes[i]);
                }
            }
        }
    }

    free(stack->array);
    free(stack);
}

void addTransitiveEdges(Graph *graph) {
    for (int i = 0; i < graph->numNodes; i++) {
        int visited[graph->numNodes];
        for (int j = 0; j < graph->numNodes; j++) {
            visited[j] = 0;
        }
        dfs(graph, graph->nodes[i], visited);
        for (int j = 0; j < graph->numNodes; j++) {
            if (i != j && visited[j] && !edgeExists(graph, graph->nodes[i], graph->nodes[j])) {
                addEdge(graph, graph->nodes[i], graph->nodes[j]);
            }
        }
    }
}

VisitedMina *createVisitedMina(double x, double y, double r) {
    VisitedMina *newMina = (VisitedMina *)malloc(sizeof(VisitedMina));
    newMina->x = x;
    newMina->y = y;
    newMina->r = r;
    newMina->next = NULL;
    return newMina;
}

void addVisitedMina(VisitedMina **head, double x, double y, double r) {
    VisitedMina *newMina = createVisitedMina(x, y, r);
    newMina->next = *head;
    *head = newMina;
}

bool isVisited(VisitedMina *head, Mina mina) {
    VisitedMina *tmp = head;
    while(tmp != NULL) {
        if(tmp->x == mina.x && tmp->y == mina.y && tmp->r == mina.r) {
            return true;
        }
        tmp = tmp->next;
    }
    return false;
}

//</editor-fold>

void ucitajNazivFajla(char *nazivFajla) {
    printf("-> ");
    scanf("%210s", nazivFajla);
}

Mina citanjeIzFajla(Graph *graph) {
    char file[210];
    printf("Unesite naziv fajla: \n");
    ucitajNazivFajla(file);
    FILE *fajl = fopen(file, "r");
    if (fajl == NULL) {
        printf("Greska pri otvaranju fajla!\n");
        exit(1);
    }

    int brMina;
    fscanf(fajl, "%d", &brMina);

    Mina *temp = realloc(graph->nodes, brMina * sizeof(Mina));
    if (temp == NULL) {
        printf("Greska pri realokaciji memorije!\n");
        exit(6);
    }
    graph->nodes = temp;

    graph->numNodes = brMina;

    double x, y, r;

    Mina start_mina;
    for (int i = 0; i < brMina; i++) {
        fscanf(fajl, "%lf %lf %lf", &x, &y, &r);
        if (i == 0) {
            start_mina.x = x;
            start_mina.y = y;
            start_mina.r = r;
        }
        ubaciUGraf(graph, x, y, r, i);
    }

    for (int i = 0; i < graph->numNodes; i++) {
        for (int j = i + 1; j < graph->numNodes; j++) {
            if (distance(graph->nodes[i], graph->nodes[j]) <= graph->nodes[i].r) {
                addEdge(graph, graph->nodes[i], graph->nodes[j]);
            }
            if (distance(graph->nodes[i], graph->nodes[j]) <= graph->nodes[j].r) {
                addEdge(graph, graph->nodes[j], graph->nodes[i]);
            }
        }
    }

    addTransitiveEdges(graph);

    fclose(fajl);

    return start_mina;
}

//<editor-fold desc="Konzola funkcije">

void dodajMinuKonzola(Graph *graph) {
    double x, y, r;
    printf("Unesite x koordinatu mine: \n");
    printf("->");
    scanf("%lf", &x);
    printf("Unesite y koordinatu mine: \n");
    printf("->");
    scanf("%lf", &y);
    printf("Unesite radijus mine: \n");
    printf("->");
    scanf("%lf", &r);
    if(ubaciUGrafSaProverom(graph, x, y, r, graph->numNodes)) {
        graph->numNodes++;
    }
}

void addEdgeKonzola(Graph *graph) {
    int sourceIndex, destinationIndex;
    printf("Unesite index pocetne mine: \n");
    printf("-> ");
    scanf("%d", &sourceIndex);
    printf("Unesite index krajnje mine: \n");
    printf("-> ");
    scanf("%d", &destinationIndex);

    if (sourceIndex < 0 || sourceIndex >= graph->numNodes || destinationIndex < 0 ||
        destinationIndex >= graph->numNodes) {
        printf("Nevalidan index.\n");
        return;
    }

    addEdge(graph, graph->nodes[sourceIndex], graph->nodes[destinationIndex]);
}

void removeEdgeKonzola(Graph *graph) {
    int sourceIndex, destinationIndex;
    printf("Unesite index pocetne mine: \n");
    printf("-> ");
    scanf("%d", &sourceIndex);
    printf("Unesite index krajnje mine: \n");
    printf("-> ");
    scanf("%d", &destinationIndex);

    if (sourceIndex < 0 || sourceIndex >= graph->numNodes || destinationIndex < 0 ||
        destinationIndex >= graph->numNodes) {
        printf("Nevalidan index mine.\n");
        return;
    }

    removeEdge(graph, graph->nodes[sourceIndex], graph->nodes[destinationIndex]);
}

void findAndPrintMostEffectiveMines(Graph *graph) {
    int maxEfficiency = 0;
    int countMaxEfficiency = 0;

    for (int i = 0; i < graph->numNodes; i++) {
        int efficiency = 0;
        for (int j = 0; j < graph->numEdges; j++) {
            if (graph->edges[j].source.x == graph->nodes[i].x && graph->edges[j].source.y == graph->nodes[i].y &&
                graph->edges[j].source.r == graph->nodes[i].r) {
                efficiency++;
            }
        }

        if (efficiency > maxEfficiency) {
            maxEfficiency = efficiency;
            countMaxEfficiency = 1;
        } else if (efficiency == maxEfficiency) {
            countMaxEfficiency++;
        }
    }

    for (int i = 0; i < graph->numNodes; i++) {
        int efficiency = 0;
        for (int j = 0; j < graph->numEdges; j++) {
            if (graph->edges[j].source.x == graph->nodes[i].x && graph->edges[j].source.y == graph->nodes[i].y &&
                graph->edges[j].source.r == graph->nodes[i].r) {
                efficiency++;
            }
        }

        if (efficiency == maxEfficiency) {
            if (countMaxEfficiency > 1) {
                if (maxEfficiency == 1) {
                    printf("Jedna od najefikasnijih mina je mina na koordinatama: (%.2lf, %.2lf) sa radijusom %.2lf. Ona ce detonirati %d minu.\n",
                           graph->nodes[i].x, graph->nodes[i].y, graph->nodes[i].r, maxEfficiency+1);
                } else {
                    printf("Jedna od najefikasnijih mina je mina na koordinatama: (%.2lf, %.2lf) sa radijusom %.2lf. Ona ce detonirati %d mine.\n",
                           graph->nodes[i].x, graph->nodes[i].y, graph->nodes[i].r, maxEfficiency+1);
                }
            } else {
                if (maxEfficiency == 1) {
                    printf("Najefikasnija mina se nalazi na koordinatama (%.2lf, %.2lf) sa radijusom %.2lf. Ona ce detonirati %d minu.\n",
                           graph->nodes[i].x, graph->nodes[i].y, graph->nodes[i].r, maxEfficiency+1);
                } else {
                    printf("Najefikasnija mina se nalazi na koordinatama (%.2lf, %.2lf) sa radijusom %.2lf. Ona ce detonirati %d mine.\n",
                           graph->nodes[i].x, graph->nodes[i].y, graph->nodes[i].r, maxEfficiency+1);
                }
            }
        }
    }
}

void findAndPrintMineEffectiveness(Graph *graph, double x, double y, double r) {
    int mineIndex = -1;
    for (int i = 0; i < graph->numNodes; i++) {
        if (graph->nodes[i].x == x && graph->nodes[i].y == y && graph->nodes[i].r == r) {
            mineIndex = i;
            break;
        }
    }

    if (mineIndex == -1) {
        printf("Mina sa zadatim koordinatama (%.2lf, %.2lf) i radijusom %.2lf ne postoji u grafu.\n", x, y, r);
        return;
    }

    int efficiency = 0;
    printf("Mina na poziciji (%.2lf, %.2lf) sa radijusom %.2lf ce detonirati sledece mine:\n", x, y, r);
    for (int i = 0; i < graph->numEdges; i++) {
        if (graph->edges[i].source.x == x && graph->edges[i].source.y == y && graph->edges[i].source.r == r) {
            printf("Minu na poziciji (%.2lf, %.2lf) sa radijusom %.2lf\n",
                   graph->edges[i].destination.x, graph->edges[i].destination.y, graph->edges[i].destination.r);
            efficiency++;
        }
    }
    printf("Ukupno ce biti detonirano %d mina.\n", efficiency+1);
}

void DFSdetonacija(Graph *graph, int index, bool *visited, int *detoniraneMine) {
    Stack *stek = createStack(graph->numNodes);
    push(stek, graph->nodes[index]);
    while (!isEmpty(stek)) {
        Mina current_mina = pop(stek);
        int current_index = getMinaIndex(graph, current_mina);
        if (!visited[current_index]) {
            visited[current_index] = true;
            (*detoniraneMine)++;
            printf("Minu na poziciji (%.2lf, %.2lf) sa radijusom %.2lf\n",
                   current_mina.x, current_mina.y, current_mina.r);
        }
        for (int i = 0; i < graph->numNodes; i++) {
            if (is_adjacent(graph, current_mina, graph->nodes[i]) && !visited[i]) {
                push(stek, graph->nodes[i]);
            }
        }
    }
    free(stek->array);
    free(stek);
}

void baciRaketu(Graph *graph, double x, double y, double r) {
    printf("Raketa je bacena na poziciju (%.2lf, %.2lf) sa radijusom %.2lf i detonirace sledece mine:\n", x, y, r);

    bool *visited = calloc(graph->numNodes, sizeof(bool));
    if (visited == NULL) {
        printf("Neuspesna alokacija memorije.\n");
        return;
    }

    int detoniraneMine = 0;
    for (int i = 0; i < graph->numNodes; i++) {
        double dx = graph->nodes[i].x - x;
        double dy = graph->nodes[i].y - y;
        double distance = sqrt(dx * dx + dy * dy);

        if (distance <= r && !visited[i]) {
            DFSdetonacija(graph, i, visited, &detoniraneMine);
        }
    }
    printf("Ukupno je detonirano %d mina.\n", detoniraneMine+1);
    free(visited);
}

//<editor-fold desc="Povrsina">
VisitedMina *add_to_visited(VisitedMina *head, Mina mina) {
    VisitedMina *new_node = (VisitedMina *) malloc(sizeof(VisitedMina));
    new_node->x = mina.x;
    new_node->y = mina.y;
    new_node->r = mina.r;
    new_node->next = head;
    return new_node;
}

int is_in_visited(VisitedMina *head, Mina mina) {
    VisitedMina *current = head;
    while (current != NULL) {
        if (current->x == mina.x && current->y == mina.y && current->r == mina.r) {
            return 1;
        }
        current = current->next;
    }
    return 0;
}

void dfsVisited(Mina *mine, Mina *start, Graph *graph, Stack *stack, double *minX, double *maxX, double *minY, double *maxY) {
    VisitedMina *visited = NULL;

    push(stack, *start);
    while (!isEmpty(stack)) {
        Mina current = pop(stack);
        if (!is_in_visited(visited, current)) {
            visited = add_to_visited(visited, current);

            if (*minX > current.x - current.r)
                *minX = current.x - current.r;
            if (*maxX < current.x + current.r)
                *maxX = current.x + current.r;
            if (*minY > current.y - current.r)
                *minY = current.y - current.r;
            if (*maxY < current.y + current.r)
                *maxY = current.y + current.r;

            for (int i = 0; i < graph->numNodes; i++) {
                if (is_adjacent(graph, current, graph->nodes[i]) && !is_in_visited(visited, graph->nodes[i])) {
                    push(stack, graph->nodes[i]);
                }
            }
        }
    }

    VisitedMina *current = visited, *next_node;
    while (current != NULL) {
        next_node = current->next;
        free(current);
        current = next_node;
    }
}

double povrsinaPokrivenaEksplozijom(Graph *graph) {
    int maxEfficiency = 0;
    Mina *mostEfficientMine = NULL;

    for (int i = 0; i < graph->numNodes; i++) {
        int efficiency = 0;
        for (int j = 0; j < graph->numEdges; j++) {
            if (graph->edges[j].source.x == graph->nodes[i].x && graph->edges[j].source.y == graph->nodes[i].y &&
                graph->edges[j].source.r == graph->nodes[i].r) {
                efficiency++;
            }
        }
        if (efficiency > maxEfficiency) {
            maxEfficiency = efficiency;
            mostEfficientMine = &graph->nodes[i];
        }
    }

    double minX = mostEfficientMine->x - mostEfficientMine->r;
    double maxX = mostEfficientMine->x + mostEfficientMine->r;
    double minY = mostEfficientMine->y - mostEfficientMine->r;
    double maxY = mostEfficientMine->y + mostEfficientMine->r;

    printf("Najefikasnija mina je na poziciji (%.2lf, %.2lf) sa radijusom %.2lf i efikasnoscu %d.\n",
           mostEfficientMine->x, mostEfficientMine->y, mostEfficientMine->r, maxEfficiency);

    Stack *stack = createStack(graph->numNodes);
    dfsVisited(mostEfficientMine, mostEfficientMine, graph, stack, &minX, &maxX, &minY, &maxY);
    free(stack->array);
    free(stack);

    int numPoints = 9000;
    int countInside = 0;
    for (int i = 0; i < numPoints; i++) {
        double x = minX + (maxX - minX) * (double) rand() / RAND_MAX;
        double y = minY + (maxY - minY) * (double) rand() / RAND_MAX;

        double dx = x - mostEfficientMine->x;
        double dy = y - mostEfficientMine->y;
        double distanceSquared = dx * dx + dy * dy;
        if (distanceSquared <= mostEfficientMine->r * mostEfficientMine->r) {
            countInside++;
        }
    }

    double totalArea = (maxX - minX) * (maxY - minY) * (double) countInside / numPoints+2;
    return totalArea;
}
//</editor-fold>

//</editor-fold>

//<editor-fold desc="Funkcije za crtanje">
void drawCircle(int x_centre, int y_centre, int r, char **matrix) {
    int x = r, y = 0;
    int P = 1 - r;
    while (x > y) {
        y++;

        if (P <= 0) {
            P = P + 2*y + 1;
        }
        else {
            x--;
            P = P + 2*y - 2*x + 1;
        }

        if (x < y) {
            break;
        }
        matrix[x_centre + x][y_centre + y] = '*';
        matrix[x_centre - x][y_centre + y] = '*';
        matrix[x_centre + x][y_centre - y] = '*';
        matrix[x_centre - x][y_centre - y] = '*';
        if (x != y) {
            matrix[x_centre + y][y_centre + x] = '*';
            matrix[x_centre - y][y_centre + x] = '*';
            matrix[x_centre + y][y_centre - x] = '*';
            matrix[x_centre - y][y_centre - x] = '*';
        }
    }
}

void printCoordinates(Graph *graph, char *filename) {
    FILE *fajl = fopen(filename, "w");
    if (fajl == NULL) {
        printf("Greska pri otvaranju fajla!\n");
        exit(5);
    }

    int minX = INT_MAX, maxX = INT_MIN, minY = INT_MAX, maxY = INT_MIN;
    int scaleFactor = 10;

    for (int i = 0; i < graph->numNodes; i++) {
        int x = (int) (scaleFactor * graph->nodes[i].x);
        int y = (int) (scaleFactor * graph->nodes[i].y);
        int r = (int) (scaleFactor * graph->nodes[i].r);

        if (x - r < minX) minX = x - r;
        if (x + r > maxX) maxX = x + r;
        if (y - r < minY) minY = y - r;
        if (y + r > maxY) maxY = y + r;
    }

    int width = maxX - minX + 1;
    int height = maxY - minY + 1;

    char **matrix = malloc(width * sizeof(char *));
    for (int i = 0; i < width; i++) {
        matrix[i] = malloc(height * sizeof(char));
        for (int j = 0; j < height; j++) {
            matrix[i][j] = ' ';
        }
    }

    for (int i = 0; i < graph->numNodes; i++) {
        int x = (int) (scaleFactor * graph->nodes[i].x) - minX;
        int y = (int) (scaleFactor * graph->nodes[i].y) - minY;
        int r = (int) (scaleFactor * graph->nodes[i].r);
        drawCircle(x, y, r, matrix);
    }

    for (int j = height - 1; j >= 0; j--) {
        for (int i = 0; i < width; i++) {
            fprintf(fajl, "%c", matrix[i][j]);
        }
        fprintf(fajl, "\n");
    }

    for (int i = 0; i < width; i++) {
        free(matrix[i]);
    }
    free(matrix);

    fclose(fajl);
}
//</editor-fold>

//<editor-fold desc="DFS ispis">

void DFSIspis(Graph *graph, Mina start_mina) {
    Stack *stack = createStack(graph->numNodes);
    VisitedMina *visited = NULL;

    int startIndex = getMinaIndex(graph, start_mina);
    if (startIndex == -1) {
        printf("Pocetna mina (start_mina) nije pronadjena u grafu.\n");
        return;
    }
    push(stack, start_mina);

    printf("DFS prolazak pocevsi od pocetne mine: (%.2lf, %.2lf) sa radijusom %.2lf:\n", start_mina.x,
           start_mina.y, start_mina.r);

    while (!isEmpty(stack)) {
        Mina current_mina = pop(stack);

        if (!isVisited(visited, current_mina)) {
            addVisitedMina(&visited, current_mina.x, current_mina.y, current_mina.r);
            printf("Posecena mina na poziciji: (%.2lf, %.2lf) sa radijusom %.2lf\n", current_mina.x, current_mina.y,
                   current_mina.r);

            for (int i = 0; i < graph->numNodes; i++) {
                if (!isVisited(visited, graph->nodes[i]) && is_adjacent(graph, current_mina, graph->nodes[i])) {
                    push(stack, graph->nodes[i]);
                }
            }
        }
    }

    for (int i = 0; i < graph->numNodes; i++) {
        if (!isVisited(visited, graph->nodes[i])) {
            push(stack, graph->nodes[i]);
            printf("DFS prolazak pocevsi od pocetne mine: (%.2lf, %.2lf) sa radijusom %.2lf:\n", graph->nodes[i].x,
                   graph->nodes[i].y, graph->nodes[i].r);

            while (!isEmpty(stack)) {
                Mina current_mina = pop(stack);

                if (!isVisited(visited, current_mina)) {
                    addVisitedMina(&visited, current_mina.x, current_mina.y, current_mina.r);
                    printf("Posecena mina na poziciji: (%.2lf, %.2lf) sa radijusom %.2lf\n", current_mina.x, current_mina.y,
                           current_mina.r);

                    for (int j = 0; j < graph->numNodes; j++) {
                        if (!isVisited(visited, graph->nodes[j]) && is_adjacent(graph, current_mina, graph->nodes[j])) {
                            push(stack, graph->nodes[j]);
                        }
                    }
                }
            }
        }
    }

    free(stack->array);
    free(stack);
    while(visited != NULL) {
        VisitedMina *tmp = visited;
        visited = visited->next;
        free(tmp);
    }
}




//</editor-fold>

void izlaz() {
    printf("\n\xDA\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xBF\n");
    printf("\xB3         Izlaz iz programa...         \xB3\n");
    printf("\xC0\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xD9\n");
    printf("----------------------------------------");
    printf("\n\xDA\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xBF\n");
    printf("\xB3       Program uspesno zavrsen!       \xB3\n");
    printf("\xC0\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xD9\n");
}

void deleteGraph(Graph *graph) {

    free(graph->nodes);
    free(graph->edges);
    free(graph);
}

int main() {

    Graph *graph = createGraph(0);
    Mina start_mina = citanjeIzFajla(graph);
    puts("Fajl je uspesno ucitan!\n");

    printGraph(graph);

//<editor-fold desc = "MENI"
    int izbor;
    do {
        printf("\n\xDA\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xBF\n");
        printf("\xB3                 MENI                 \xB3\n");
        printf("\xC3\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xB4\n");
        printf("\xB3  1. DODAJ CVOR U GRAF                \xB3\n");
        printf("\xC3\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xB4\n");
        printf("\xB3  2. UKLONI CVOR IZ GRAFA             \xB3\n");
        printf("\xC3\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xB4\n");
        printf("\xB3  3. DODAJ GRANU U GRAF               \xB3\n");
        printf("\xC3\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xB4\n");
        printf("\xB3  4. UKLONI GRANU IZ GRAFA            \xB3\n");
        printf("\xC3\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xB4\n");
        printf("\xB3  5. DFS ISPIS GRAFA                  \xB3\n");
        printf("\xC3\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xB4\n");
        printf("\xB3  6. BACI RAKETU                      \xB3\n");
        printf("\xC3\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xB4\n");
        printf("\xB3  7. EFIKASNOST ZADATE MINE           \xB3\n");
        printf("\xC3\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xB4\n");
        printf("\xB3  8. PRONADJI MINU MAX EFIKASNOSTI    \xB3\n");
        printf("\xC3\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xB4\n");
        printf("\xB3  9. POVRSINA MINE SA MAX EFIKASNOSCU \xB3\n");
        printf("\xC3\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xB4\n");
        printf("\xB3  10. ISPISI GRAF SA KOORDINATAMA     \xB3\n");
        printf("\xC3\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xB4\n");
        printf("\xB3  11. Izlaz                           \xB3\n");
        printf("\xC0\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xC4\xD9\n");
        printf("Unesite broj za odabir zeljene opcije: \n");
        printf("-> ");
        scanf("%d", &izbor);

        switch (izbor) {
            case 1:
                printf("Izabrali ste Opciju \"DODAJ CVOR U GRAF\"\n");
                dodajMinuKonzola(graph);
                break;
            case 2:
                printf("Izabrali ste Opciju \"UKLONI CVOR IZ GRAFA\"\n");
                double x, y, r;
                printf("Unesite x koordinatu mine: \n");
                printf("-> ");
                scanf("%lf", &x);
                printf("Unesite y koordinatu mine: \n");
                printf("-> ");
                scanf("%lf", &y);
                printf("Unesite radius mine: \n");
                printf("-> ");
                scanf("%lf", &r);
                Mina mina = {x, y, r};
                removeNode(graph, mina);
                break;
            case 3:
                printf("Izabrali ste Opciju \"DODAJ GRANU U GRAF\"\n");
                addEdgeKonzola(graph);
                break;
            case 4:
                printf("Izabrali ste Opciju \"UKLONI GRANU IZ GRAFA\"\n");
                removeEdgeKonzola(graph);
                printf("Uspesno ste uklonili granu iz grafa!\n");
                break;
            case 5:
                printf("Izabrali ste Opciju \"DFS ISPIS GRAFA\"\n");
                DFSIspis(graph, start_mina);
                printGraph(graph);
                break;
            case 6:
                printf("Izabrali ste Opciju \"BACI RAKETU\"\n");
                double x2, y2, r2;
                printf("Unesite x koordinatu mine: \n");
                printf("-> ");
                scanf("%lf", &x2);
                printf("Unesite y koordinatu mine: \n");
                printf("-> ");
                scanf("%lf", &y2);
                printf("Unesite radius mine: \n");
                printf("-> ");
                scanf("%lf", &r2);
                baciRaketu(graph, x2, y2, r2);
                break;
            case 7:
                printf("Izabrali ste Opciju \"EFIKASNOST ZADATE MINE\"\n");
                double x1, y1, r1;
                printf("Unesite x koordinatu mine: \n");
                printf("-> ");
                scanf("%lf", &x1);
                printf("Unesite y koordinatu mine: \n");
                printf("-> ");
                scanf("%lf", &y1);
                printf("Unesite radius mine: \n");
                printf("-> ");
                scanf("%lf", &r1);
                findAndPrintMineEffectiveness(graph, x1, y1, r1);
                break;
            case 8:
                printf("Izabrali ste Opciju \"PRONADJI MINU MAX EFIKASNOSTI\"\n");
                findAndPrintMostEffectiveMines(graph);
                break;
            case 9:
                printf("Izabrali ste Opciju \"POVRSINA MINE SA MAX EFIKASNOSCU\"\n");
                printf("Povrsina eksplozije te mine je: %.4lf\n", povrsinaPokrivenaEksplozijom(graph));
                break;
            case 10:
                printf("Izabrali ste Opciju \"ISPISI GRAF SA KOORDINATAMA\"\n");
                printf("Unesite naziv fajla u kojem zelite da se iscrta (pozeljno je da ne koristite fajlove koji se vec koriste): \n");
                char nazivFajla[210];
                ucitajNazivFajla(nazivFajla);
                printCoordinates(graph, nazivFajla);
                printf("Uspesno ste ispisali graf sa koordinatama u fajl \"%s\"!\n", nazivFajla);
                break;
            case 11:
                izlaz();
                break;
            default:
                printf("Pogresan izbor, molimo unesite broj 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 ili 11.\n");
        }

    } while (izbor != 11);
//</editor-fold>

    deleteGraph(graph);

    return 0;
}