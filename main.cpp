#include <iostream>
#include <list>
#include <vector>
#include <cstdlib>
#include <limits>
#include <cctype>
#include <iomanip>

using Vertex = unsigned int;
using uint = unsigned int;

template<typename T>
class PriorityQueue {
private:
    std::vector<T> heap;
    std::vector<int> &dist; //referencia para o vetor de distancias

    int parent(int i) { return (i - 1) / 2; }
    int left(int i) { return 2 * i + 1; }
    int right(int i) { return 2 * i + 2; }

    void sift_up(int i);
    void sift_down(int i);

public:
    // Construtor que recebe o vetor de distancias
    PriorityQueue(std::vector<int> &dist) : dist(dist) {}
    bool empty() const { return heap.empty();}
    int size() const { return heap.size(); }
    void insert(T element);
    T extract_min();

};

template<typename T>
void PriorityQueue<T>::sift_up(int i) {
    while (i > 0 && dist[heap[i]] < dist[heap[parent(i)]]) {
        std::swap(heap[i], heap[parent(i)]);
        i = parent(i);
    }
}

template<typename T>
void PriorityQueue<T>::sift_down(int i) {
    int n = heap.size();
    while (true) {
        int smallest = i;
        int l = left(i);
        int r = right(i);

        if (l < n && dist[heap[l]] < dist[heap[smallest]]) {
            smallest = l;
        }
        if (r < n && dist[heap[r]] < dist[heap[smallest]]) {
            smallest = r;
        }
        if (smallest != i) {
            T temp = heap[i];
            heap[i] = heap[smallest];
            heap[smallest] = temp;

            i = smallest;
        } else {
            break;
        }
    }
}

template<typename T>
void PriorityQueue<T>::insert(T element) {
    heap.push_back(element);
    sift_up(heap.size() - 1);

}

template<typename T>
T PriorityQueue<T>::extract_min() {
    if (heap.empty()) {
        throw std::runtime_error("Heap vazio");
    }
    T min = heap[0];
    heap[0] = heap.back();
    heap.pop_back();

    if (!heap.empty()) {
        sift_down(0);
    }
    return min;
}

class Graph_Board {
private:
    uint num_vertices;
    uint N;

    std::vector<std::vector<std::pair<Vertex, int>>> adj; //par do adj e peso
    std::vector<std::pair<int, int>> moviments = {
        {1, 2}, {1, -2}, {-1, 2}, {-1, -2}, {2, 1}, {2, -1}, {-2, 1}, {-2, -1}
    };


    void createGraph();
    void insertionSort(std::vector<std::pair<Vertex, int> > &vec);

public:
    Graph_Board(uint N = 8);
    ~Graph_Board() = default;
    void add_edge(Vertex u, Vertex v);
    int calculateWeight(int u, int v, int N);
    const std::vector<std::pair<Vertex, int>> &get_adj(Vertex u) const { return adj[u]; }
    uint get_vertices() { return num_vertices; }
    uint get_size() const { return N; }
};

void Graph_Board::add_edge(Vertex u, Vertex v) {
    for (const auto &neighbor : adj[u]) {
        if (neighbor.first == v) {
            return;
        }
    }
    int weight = calculateWeight(u, v, N);
    adj[u].push_back({v, weight});
    adj[v].push_back({u, weight});
}

Graph_Board::Graph_Board(uint N) : num_vertices(N * N), N(N) {
    adj.resize(num_vertices);
    createGraph();
}

void Graph_Board::insertionSort(std::vector<std::pair<Vertex, int>> &vec) {
    for (size_t i = 1; i < vec.size(); i++) {
        auto key = vec[i];
        int j = i - 1;
        while (j >= 0 && vec[j].first > key.first) {
            vec[j + 1] = vec[j];
            j--;
        }
        vec[j + 1] = key;
    }
}

//static aqui pois o hackerhank da aviso de compilacao
void Graph_Board::createGraph() {
    for (uint line = 0; line < N; ++line) {
        for (uint column = 0; column < N; ++column) {
            uint vertex_origin = line * N + column;
            for (const auto &mov: moviments) {
                int new_line = line + mov.first;
                int new_column = column + mov.second;
                if (new_line >= 0 && new_line < static_cast<int>(N)  &&
                    new_column >= 0 && new_column < static_cast<int>(N) ) {
                    uint vertex_destiny = new_line * static_cast<int>(N)  + new_column;
                    if (vertex_origin < vertex_destiny) {
                        add_edge(vertex_origin, vertex_destiny);
                    }
                    }
            }
        }
    }
    for (uint u = 0; u < num_vertices; ++u) {
        if (!adj[u].empty()) {
            insertionSort(adj[u]);
        }
    }
}

int Graph_Board::calculateWeight(int u, int v, int N) {
    int line_u = u / N;
    int column_u = u % N;

    char alpha_u = 'a' + column_u; //au
    int beta_u = line_u + 1; //bu

    int line_v = v / N;
    int column_v = v % N;

    char alpha_v = 'a' + column_v; //av
    int beta_v = line_v + 1; //bv
    return ((int(alpha_u) * beta_u + int(alpha_v) * beta_v) % 19);
};

// Estruturas auxiliaresr
struct resultPath {
    int dist;
    std::vector<int> path;
};

struct Army {
    std::string color;
    std::string position;
    std::vector<std::string> enemies;
};

struct resultArmy {
    std::string color;
    int moviments;
    int weight;
};

struct EstadoExercito {
    std::string color;
    std::vector<int> path;
    int pos_idx = 0;
    int moviments = 0;
    int total_weight = 0;
    bool arrived = false;
    bool stopped = false;
    std::vector<std::string> enemies;
    std::vector<std::string> allys;
    int position = -1;
    int arrival_round = -1;
};

class ArmysAttack {
private:
    std::vector<Army> armys;
    std::string castelPosition;
    std::vector<std::string> storms;

    std::pair<int, int> ChessNotToPos(const std::string &position);
    bool isStorm(int vertex, int N);
    std::vector<int> reconstructPath(int destiny, const std::vector<int> &parent);
    int chooseDestiny(const std::vector<int> &list_threats, const std::vector<int> &dist, int &less_distancy);
    resultPath bestPathToCastle(Graph_Board &game_board, int at_beggining, const std::vector<int> &list_threats);
    std::pair<std::vector<int>, std::vector<int>> Dijkstra(Graph_Board &game_board, int at_beggining);
    std::vector<EstadoExercito> startArmys(Graph_Board &game_graph, int N, int castel_vertex);
    void move(EstadoExercito &state_army, int next_position, int weight, int castel_vertex, int N, int round);
    bool tryEnemy(EstadoExercito &state, int next_position, std::vector<EstadoExercito> &stateArmys);
    void executeRound(std::vector<EstadoExercito> &state_army, Graph_Board &game_graph, int round, int castle_vertex);
    std::vector<resultArmy> processResult(const std::vector<EstadoExercito> &state_army);
    void insertionSortArmy(std::vector<resultArmy> &vec);
    bool canBeAlly(const EstadoExercito &a, const EstadoExercito &b);
    bool contains(const std::vector<std::string>& vec, const std::string& value);

public:
    ArmysAttack(std::vector<Army> listArmy, std::string castelPosition, std::vector<std::string> storms);
    void solve_game(Graph_Board &game_board);
};


bool ArmysAttack::canBeAlly(const EstadoExercito &a, const EstadoExercito &b) {
    if (contains(a.enemies, b.color)) {
        return false;
    }
    if (contains(b.enemies, a.color)) {
        return false;
    }
    return true;
}

ArmysAttack::ArmysAttack(std::vector<Army> la, std::string cp, std::vector<std::string> s) : armys(la), castelPosition(cp), storms(s) {}

std::pair<int, int> ArmysAttack::ChessNotToPos(const std::string &pos) {
    return {std::stoi(pos.substr(1)) - 1, pos[0] - 'a'};
}

bool ArmysAttack::isStorm(int vertex, int N) {
    for (auto &s : storms) {
        auto [line, row] = ChessNotToPos(s);
        if (vertex == line * N + row) return true;
    }
    return false;
}

std::vector<int> ArmysAttack::reconstructPath(int d, const std::vector<int> &p) {
    std::vector<int> path;
    for (int v = d; v != -1; v = p[v]) {
        path.insert(path.begin(), v);
    }
    return path;
}

int ArmysAttack::chooseDestiny(const std::vector<int> &list_threats, const std::vector<int> &dist, int &less_distancy) {
    int final_destiny = -1;
    less_distancy = 1e9;
    for (int target_pos : list_threats) {
        if (dist[target_pos] < less_distancy) {
            less_distancy = dist[target_pos];
            final_destiny = target_pos;
        }
    }
    return final_destiny;
}

resultPath ArmysAttack::bestPathToCastle(Graph_Board &game_board, int at_beggining, const std::vector<int> &list_threats) {
    auto [dist, parent] = Dijkstra(game_board, at_beggining);
    int less_distance;
    int final_destiny = chooseDestiny(list_threats, dist, less_distance);

    if (final_destiny != -1) {
        return {less_distance, reconstructPath(final_destiny, parent)};
    }
    return {-1, {}};
}

std::pair<std::vector<int>, std::vector<int>> ArmysAttack::Dijkstra(Graph_Board &game_board, int at_beggining) {
    int V = game_board.get_vertices();
    std::vector<int> dist(V, 1e9);
    std::vector<int> parent(V, -1);
    std::vector<int> moves(V, 1e9);
    PriorityQueue<int> pq(dist);

    dist[at_beggining] = 0;
    moves[at_beggining] = 0;
    pq.insert(at_beggining);

    while (!pq.empty()) {
        int u = pq.extract_min();
        if (dist[u] == 1e9) {
            continue;
        }

        for (const auto &[neighbor, weight] : game_board.get_adj(u)) {
            if (dist[u] + weight < dist[neighbor]) {
                dist[neighbor] = dist[u] + weight;
                parent[neighbor] = u;
                moves[neighbor] = moves[u] + 1;
                pq.insert(neighbor);
            } else if (dist[u] + weight == dist[neighbor]) {
                if (moves[u] + 1 < moves[neighbor]) {
                    parent[neighbor] = u;
                    moves[neighbor] = moves[u] + 1;
                }
            }
        }
    }
    return {dist, parent};
}

std::vector<EstadoExercito> ArmysAttack::startArmys(Graph_Board &game_graph, int N, int castel_vertex) {
    std::vector<EstadoExercito> estados;
    for (auto &exercito : armys) {
        auto [line, row] = ChessNotToPos(exercito.position);
        resultPath result = bestPathToCastle(game_graph, line * N + row, {castel_vertex});

        EstadoExercito estado;
        estado.color = exercito.color;
        estado.path = result.path;
        estado.position = line * N + row;
        estado.enemies = exercito.enemies;
        estados.push_back(estado);
    }
    return estados;
}

void ArmysAttack::move(EstadoExercito &state_army, int next_position, int weight, int castel_vertex, int N, int round) {
    state_army.moviments++;
    state_army.total_weight += weight;
    state_army.pos_idx++;
    state_army.position = next_position;

    if (next_position == castel_vertex && !state_army.arrived) {
        state_army.arrived = true;
        state_army.arrival_round = round;
    }
}

bool ArmysAttack::contains(const std::vector<std::string>& vec, const std::string& value) {
    for (size_t i = 0; i < vec.size(); ++i) {
        if (vec[i] == value) {
            return true;
        }
    }
    return false;
}

bool ArmysAttack::tryEnemy(EstadoExercito &state, int next_position, std::vector<EstadoExercito> &stateArmys) {
    for (auto &other : stateArmys) {
        if (&other != &state && other.position == next_position && !other.arrived) {
            if (!canBeAlly(state, other)) {
                state.stopped = true;
                return true;
            } else {
                if (!contains(state.allys, other.color)) {
                    state.allys.push_back(other.color);
                }
                if (!contains(other.allys, state.color)) {
                    other.allys.push_back(state.color);
                }
            }
        }
    }
    return false;
}

void ArmysAttack::executeRound(std::vector<EstadoExercito> &state_army, Graph_Board &game_graph, int round, int castle_vertex) {
    int N = game_graph.get_size();
    for (auto &estado : state_army) {
        if (estado.arrived || estado.stopped) continue;
        if (static_cast<size_t>(estado.pos_idx + 1) >= estado.path.size()) continue;

        int next_pos = estado.path[estado.pos_idx + 1];

        if (isStorm(next_pos, N)) {
            bool is_ally= false;
            for(const auto& outro : state_army) {
                if(&outro != &estado && outro.position == estado.position && canBeAlly(estado, outro)) {
                    is_ally = true;
                    break;
                }
            }
            if(!is_ally){
                 estado.stopped = true;
                 continue;
            }
        }

        if (tryEnemy(estado, next_pos, state_army)) continue;

        move(estado, next_pos, game_graph.calculateWeight(estado.position, next_pos, N), castle_vertex, N, round);
    }

    for(auto& estado : state_army) if(estado.stopped) estado.stopped = false;
}

std::vector<resultArmy> ArmysAttack::processResult(const std::vector<EstadoExercito> &state_army) {
    std::vector<resultArmy> results;
    std::vector<const EstadoExercito*> arrived;
    for(const auto& e : state_army) if(e.arrived) arrived.push_back(&e);

    if (arrived.empty()) return results;

    int min_round = 1e9;
    for(const auto* v : arrived) {
        if(v->arrival_round < min_round) {
            min_round = v->arrival_round;
        }
    }
    for(const auto* final_state : arrived) {
        if(final_state->arrival_round == min_round) {
            results.push_back({final_state->color, final_state->moviments, final_state->total_weight});
        }
    }
    insertionSortArmy(results);
    return results;
}

void ArmysAttack::insertionSortArmy(std::vector<resultArmy> &vec) {
    for (size_t i = 1; i < vec.size(); i++) {
        auto key = vec[i];
        int j = i - 1;
        while (j >= 0 && vec[j].color > key.color) {
            vec[j + 1] = vec[j];
            j--;
        }
        vec[j + 1] = key;
    }
}

void ArmysAttack::solve_game(Graph_Board &game_board) {
    int N = game_board.get_size();
    auto [l, c] = ChessNotToPos(castelPosition);
    int castel_vertex = l * N + c;

    auto states = startArmys(game_board, N, castel_vertex);

    for (int rounds = 1; rounds <= 100; ++rounds) {
        executeRound(states, game_board, rounds, castel_vertex);
        bool all_arrived_or_stopped = true;
        for(const auto& estado : states) {
            if (!estado.arrived) {
                all_arrived_or_stopped = false;
                break;
            }
        }
        if(all_arrived_or_stopped) break;
    }

    auto results = processResult(states);

    for (size_t i = 0; i < results.size(); ++i) {
        std::cout << results[i].color << " " << results[i].moviments << " " << results[i].weight << (i == results.size() - 1 ? "" : " ");
    }
}

int main() {
    int line;
    std::cin >> line;

    if (line < 8 || line > 15) return 1;

    int armyNumber;
    std::cin >> armyNumber;
    std::cin.ignore();

    std::vector<Army> listArmys;
    for (int i = 0; i < armyNumber; i++) {
        std::string ln;
        std::getline(std::cin, ln);
        std::istringstream iss(ln);
        Army army;
        iss >> army.color >> army.position;
        std::string enemy;

        while (iss >> enemy) {
            army.enemies.push_back(enemy);
        }
        listArmys.push_back(army);
    }

    std::string castelPosition;
    std::cin >> castelPosition;

    int stormsNumber;
    std::cin >> stormsNumber;
    std::vector<std::string> stormsPositions(stormsNumber);
    for (int i = 0; i < stormsNumber; ++i) {
        std::cin >> stormsPositions[i];
    }

    Graph_Board tab(line);
    ArmysAttack jogo(listArmys, castelPosition, stormsPositions);
    jogo.solve_game(tab);

    return 0;
}