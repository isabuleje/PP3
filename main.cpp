#include <iostream>
#include <list>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <sstream>


using Vertex = unsigned int;
using uint = unsigned int;

template<typename T>
class PriorityQueue {
private:
    std::vector<T> heap;
    std::vector<int> &dist; // Referencia para o vetor de distancias

    int parent(int i) { return (i - 1) / 2; }
    int left(int i) { return 2 * i + 1; }
    int right(int i) { return 2 * i + 2; }

    void sift_up(int i) {
        while (i > 0 && dist[heap[i]] < dist[heap[parent(i)]]) {
            std::swap(heap[i], heap[parent(i)]);
            i = parent(i);
        }
    }

    void sift_down(int i) {
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
                std::swap(heap[i], heap[smallest]);
                i = smallest;
            } else {
                break;
            }
        }
    }

public:
    // Construtor que recebe o vetor de distancias
    PriorityQueue(std::vector<int> &distancias) : dist(distancias) {
    }

    bool empty() const { return heap.empty(); }
    int size() const { return heap.size(); }

    void insert(T element);

    T minimum();

    T extract_min();

    void decrease_key(T element);

    bool contains(T element) const {
        return std::find(heap.begin(), heap.end(), element) != heap.end();
    }
};

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

template<typename T>
T PriorityQueue<T>::minimum() {
    if (heap.empty()) {
        throw std::runtime_error("Heap vazio");
    }
    return heap[0];
}

template<typename T>
//Essa funcaor e mt lentar, dps me lembro de trocar
// Metodo para atualizar a prioridade (quando a distancia muda)
void PriorityQueue<T>::decrease_key(T element) {
    for (size_t i = 0; i < heap.size(); i++) {
        if (heap[i] == element) {
            sift_up(i);
            break;
        }
    }
}

// -x-x-x-x- Classe Graph_Board baseada num grafo lista de adjacencias  -x-x-x-x-
class Graph_Board {
private:
    uint num_vertices;
    uint num_edges;
    uint N; //dimensao do tabuleiro

    //tirei o ponteiro pq tava foda
    std::vector<std::vector<std::pair<Vertex, int> > > adj; //aq virou um par pois guarda o adj e o seu weight
    std::vector<std::pair<int, int> > moviments = {
        {1, 2}, {1, -2}, {-1, 2}, {-1, -2}, {2, 1}, {2, -1}, {-2, 1}, {-2, -1}
    };

    void criarGrafo();

    void insertionSort(std::vector<std::pair<Vertex, int> > &vec);

public:
    Graph_Board(uint N = 8);

    ~Graph_Board() = default;

    void add_edge(Vertex u, Vertex v);

    uint get_edges();

    const std::vector<std::pair<Vertex, int> > &get_adj(Vertex u) const;

    uint get_vertices();

    uint get_size() const { return N; };

    int calculateWeight(int u, int v, int N);
};

//construtor do Grafo tabuleiro
Graph_Board::Graph_Board(uint N) : num_vertices(N*N), num_edges(0), N(N)  {
    if (num_vertices == 0) {
        throw std::invalid_argument("Numero de vertices nao pode ser zero");
    }

    adj.resize(num_vertices);
    criarGrafo(); //cria o grafo assim que inicializa
}


//calculadora de weight q o prof pediu
//w(u,v)=(ascii(au)_Bu+ascii(av)_Bv)mod19 (formula)
int Graph_Board::calculateWeight(int u, int v, int N) {
    int linha_u = u / N;
    int coluna_u = u % N;

    char alpha_u = 'a' + coluna_u; //au
    int beta_u = linha_u + 1; //bu

    int linha_v = v / N;
    int coluna_v = v % N;

    char alpha_v = 'a' + coluna_v; //av
    int beta_v = linha_v + 1; //bv

    //esse int(letra) = ascii(letra) com mod 19
    return ((int(alpha_u) * beta_u + int(alpha_v) * beta_v) % 19);
}

void Graph_Board::add_edge(Vertex u, Vertex v) {
    if (u >= num_vertices || v >= num_vertices || u == v) {
        throw std::invalid_argument("Valores invalidos no add_edge");
    }

    //verificar se o vizinho ja existe
    for (const auto &neighbor: adj[u]) {
        if (neighbor.first == v) {
            return;
        }
    }

    int weight = calculateWeight(u, v, N);

    //aq o adj cria uma nova aresta
    adj[u].push_back({v, weight});
    adj[v].push_back({u, weight});
    num_edges += 1;
}

const std::vector<std::pair<Vertex, int> > &Graph_Board::get_adj(Vertex u) const {
    if (u >= num_vertices) {
        throw std::invalid_argument("Vertice invalido no get_adj: ");
    }
    return adj[u];
}

uint Graph_Board::get_edges() {
    return num_edges;
}

uint Graph_Board::get_vertices() {
    return num_vertices;
}

//como agr o tamanho do tabuleiro pode ser maior que 8 e melhor usar o N
void Graph_Board::criarGrafo() {
    num_edges = 0;

    for (uint line = 0; line < N; ++line) {
        for (uint column = 0; column < N; ++column) {
            uint vertex_origin = line * N + column;
            for (const auto &mov: moviments) {
                uint new_line = line + mov.first;
                uint new_column = column + mov.second;
                if (new_line >= 0 && new_line < N &&
                    new_column >= 0 && new_column < N) {
                    uint vertex_destiny = new_line * N + new_column;
                    if (vertex_origin < vertex_destiny) {
                        add_edge(vertex_origin, vertex_destiny);
                    }
                }
            }
        }
    }

    //aq usa o sort feito manualmente
    for (uint u = 0; u < num_vertices; ++u) {
        if (!adj[u].empty()) {
            insertionSort(adj[u]);
        }
    }
}

void Graph_Board::insertionSort(std::vector<std::pair<Vertex, int> > &vec) {
    if (vec.empty()) {
        return;
    }

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

// ----------------- Estrutura auxiliar para guardar o path -----------------
// ta aqui a struct que guarda a distancia e o path percorrido
struct resultPath {
    int distancia;
    std::vector<int> caminho;
};

//struct que guarda as informacoes do exercito
struct Army {
    std::string color;
    std::string position;
    std::vector<std::string> friends;

    std::vector<int> path;
    int idx_path = 0;
    bool arrived = false;
    bool stopped = false;
};

//struct que guarda as informacoes para printar no final
struct resultArmy {
    std::string color;
    int moviments;
    int weight;
};

//struct que guarda as condicoes do exercito
//dps n esquecer de unir o Army e o EstadoExercito
struct EstadoExercito {
    std::string color;
    std::vector<int> path;
    int pos_idx = 0;
    int moviments = 0;
    int total_weight = 0;
    bool arrived = false;
    bool stopped = false;
    std::vector<std::string> allys;
    int position = -1;
};


// -x-x-x-x- Classe AtaqueDosExercitos  -x-x-x-x-
class ArmysAttack {
private:
    std::vector<Army> armys;
    std::string castelPosition;
    std::vector<std::string> storms;

    //metodos (ta gigante D: pq eu modularizei tudo q conseguia)
    std::vector<int> castleChess(const Graph_Board &game_board);

    bool isStorm(int v, int N);

    resultPath bestPathToCastle(Graph_Board &game_board, int at_beggining,
                                                   const std::vector<int> &list_threats);

    std::pair<std::vector<int>, std::vector<int> > Dijkstra(Graph_Board &game_board, int at_beggining);

    int chooseDestiny(const std::vector<int> &list_threats, const std::vector<int> &dist,
                                   int &less_distancy);

    std::vector<int> reconstructPath(int destiny, const std::vector<int> &parent);

    bool tryStorm(EstadoExercito &estado, int next_position, int N);

    std::vector<EstadoExercito> startArmys(Graph_Board &game_graph, int N, int castel_vertex);

    bool tryEnemy(EstadoExercito &estado, int next_position, std::vector<EstadoExercito> &stateArmys);

    void move(EstadoExercito &estado, int next_position, int peso, int castel_vertex, int N);

    bool executeRound(std::vector<EstadoExercito> &estados, Graph_Board &game_graph, int round, int castle_vertex,
                      int &less_moviments_arrival);

    std::vector<resultArmy>
    processResult(const std::vector<EstadoExercito> &estados, int less_moviments_arrival);

    std::pair<int, int> ChessNotToPos(const std::string &position);

public:
    ArmysAttack(std::vector<Army> listArmy, std::string castelPosition, std::vector<std::string> storms);

    void solve_game(Graph_Board &game_board);
};

ArmysAttack::ArmysAttack(std::vector<Army> listArmy, std::string castelPosition, std::vector<std::string> storms) {
    this->castelPosition = castelPosition;
    this->armys = listArmy;
    this->storms = storms;
}

//converte notacao tipo "a1" em coordenadas (linha, coluna)
std::pair<int, int> ArmysAttack::ChessNotToPos(const std::string &position) {
    if (position.size() < 2) {
        throw std::invalid_argument("Posicao invalida");
    }

    //aqui precisa do stoi pq se for a10, a11 etc da errador
    int coluna = position[0] - 'a';
    int linha = std::stoi(position.substr(1)) - 1;

    return {linha, coluna};
}

std::vector<EstadoExercito> ArmysAttack::startArmys(Graph_Board &game_graph, int N, int castel_vertex) {
    std::vector<EstadoExercito> estados;
    for (auto &exercito: armys) {
        auto [l, c] = ChessNotToPos(exercito.position);
        int v = l * N + c;
        resultPath res = bestPathToCastle(game_graph, v, {castel_vertex});

        EstadoExercito estado;
        estado.color = exercito.color;
        estado.path = res.caminho;
        estado.position = v;
        estado.allys = exercito.friends; // copia lista de allys

        estados.push_back(estado);
    }
    return estados;
}

//determina como a tempestade do proximo movimento sera tratada
//se nao estiver junto de um aliado ele para, mas se estiver destroe a tempestade na hora
bool ArmysAttack::tryStorm(EstadoExercito &estado, int next_position, int N) {
    if (isStorm(next_position, N)) {
        if (estado.allys.empty()) {
            estado.stopped = true;
            return true; // perde rodada
        } else {
            // remove tempestade
            storms.erase(std::remove_if(storms.begin(), storms.end(),
                                        [&](const std::string &s) {
                                            auto [l, c] = ChessNotToPos(s);
                                            return l == next_position / N && c == next_position % N;
                                        }), storms.end());
        }
    }
    return false;
}

//gera as posicoes que sao vizinhas do castelo (um cavalo pode alcancar)
std::vector<int> ArmysAttack::castleChess(const Graph_Board &game_board) {
    //moviments dos exercitos em formato de L (fazuelli)
    const std::vector<std::pair<int, int> > movimentos = {
        {2, 1}, {2, -1}, {-2, 1}, {-2, -1}, {1, 2}, {1, -2}, {-1, 2}, {-1, -2}
    };

    int column = castelPosition[0] - 'a';
    int line = std::stoi(castelPosition.substr(1)) - 1;

    //a listSus e uma possivel lista de ameacas ao castelo
    std::vector<int> listSus;

    int board_size = static_cast<int>(game_board.get_size());

    for (auto mov: movimentos) {
        int new_line = line + mov.first;
        int new_column = column + mov.second;

        if (new_line >= 0 && new_line < board_size  && new_column >= 0 && new_column < board_size) {
            listSus.push_back(new_line * board_size  + new_column);
        }
    }

    return listSus;
}

//verificacao se tem Storm na posicao do tabuleiro
bool ArmysAttack::isStorm(int v, int N) {
    int linha = v / N;
    int coluna = v % N;

    for (auto &s: storms) {
        auto [l, c] = ChessNotToPos(s);
        if (l == linha && c == coluna) return true;
    }
    return false;
}

//incrementa os atributos do estado do exercito
//verifica se o exercito chega no castelo no proximo movimento
void ArmysAttack::move(EstadoExercito &estado, int next_position, int peso, int castel_vertex, int N) {
    estado.moviments++;
    estado.total_weight += peso;
    estado.pos_idx++;
    estado.position = next_position;

    if (next_position == castel_vertex) {
        estado.arrived = true;
    }
}

//verifica se o exercito e amiguinho ou nao
bool isAlly(const EstadoExercito &estado, const std::string &color) {
    return std::find(estado.allys.begin(), estado.allys.end(), color) != estado.allys.end();
}

//funcao principal que chama as outras e verifica os estados etcr
bool ArmysAttack::executeRound(std::vector<EstadoExercito> &estados, Graph_Board &game_graph,
                               int round, int castle_vertex, int &less_moviments_arrival) {
    int N = game_graph.get_size();
    bool moviment_happened = false;
    std::vector<std::string> arrived_this_round;

    std::cout << "RODADA " << round << " --------\n";

    // Processa todos os exercitos nesta rodada
    for (auto &estado: estados) {
        std::cout << "Exercito: " << estado.color
                  << " | Posicao atual: " << estado.position
                  << " | Movimentos: " << estado.moviments
                  << " | Status: " << (estado.arrived ? "CHEGOU" : estado.stopped ? "PARADO" : "EM MOVIMENTO")
                  << "\n";

        if (estado.arrived) continue;

        if (estado.stopped) {
            estado.stopped = false; // libera o stop
            moviment_happened = true;
            continue;
        }

        if (estado.pos_idx + 1 >= (int)estado.path.size()) continue;

        int next_position = estado.path[estado.pos_idx + 1];

        if (tryStorm(estado, next_position, N)) {
            std::cout << "  -> Encontrou tempestade, parada obrigatoria.\n";
            continue;
        }

        if (tryEnemy(estado, next_position, estados)) {
            std::cout << "  -> Encontrou inimigo, parada obrigatoria ou formacao de alianca.\n";
            continue;
        }

        int weight = game_graph.calculateWeight(estado.position, next_position, N);
        move(estado, next_position, weight, castle_vertex, N);

        std::cout << "  -> Movido para " << next_position
                  << ", peso do movimento: " << weight
                  << "peso total:" << estado.total_weight
                  << ", status agora: " << (estado.arrived ? "CHEGOU" : "EM MOVIMENTO")
                  << "\n";

        if (estado.arrived) {
            arrived_this_round.push_back(estado.color);
            less_moviments_arrival = std::min(less_moviments_arrival, estado.moviments);
        }

        moviment_happened = true;
    }

    // Verifica se ainda existe algum exercito que pode se mover
    bool any_movable_left = false;
    for (auto &estado: estados) {
        if (!estado.arrived && !estado.stopped && estado.pos_idx + 1 < (int)estado.path.size()) {
            any_movable_left = true;
            break;
        }
    }

    std::cout << "Fim da RODADA " << round << "\n\n";

    return any_movable_left;
}


std::vector<resultArmy> ArmysAttack::processResult(const std::vector<EstadoExercito> &estados,
                                                          int less_moviments_arrival) {
    std::vector<resultArmy> results;
    for (auto &estado: estados) {
        if (estado.arrived && estado.moviments == less_moviments_arrival) {
            results.push_back({estado.color, estado.moviments, estado.total_weight});
        }
    }
    std::sort(results.begin(), results.end(),
              [](const resultArmy &a, const resultArmy &b) { return a.color < b.color; });
    return results;
}

std::vector<int> ArmysAttack::reconstructPath(int destiny, const std::vector<int> &parent) {
    std::vector<int> path;
    for (int v = destiny; v != -1; v = parent[v]) {
        path.push_back(v);
    }
    std::reverse(path.begin(), path.end());
    return path;
}

int ArmysAttack::chooseDestiny(const std::vector<int> &list_threats,
                               const std::vector<int> &dist,
                               int &less_distancy) {
    int final_destiny = -1;
    less_distancy = 1e9;
    for (int pos_alvo: list_threats) {
        if (dist[pos_alvo] < less_distancy) {
            less_distancy = dist[pos_alvo];
            final_destiny = pos_alvo;
        }
    }
    return final_destiny;
}

// Dijkstra
std::pair<std::vector<int>, std::vector<int>> ArmysAttack::Dijkstra(Graph_Board &game_board, int at_beggining) {
    int V = game_board.get_vertices();
    int N = game_board.get_size();

    std::vector<int> dist(V, 1e9);
    std::vector<int> parent(V, -1);
    std::vector<bool> processed(V, false);
    std::vector<bool> inHeap(V, false);

    PriorityQueue<int> pq(dist);

    dist[at_beggining] = 0;
    pq.insert(at_beggining);
    inHeap[at_beggining] = true;

    while (!pq.empty()) {
        int u = pq.extract_min();
        inHeap[u] = false;

        if (processed[u]) continue;
        processed[u] = true;

        try {
            const auto &neighbors = game_board.get_adj(u);
            for (const auto &[neighbor, weight]: neighbors) {
                if (isStorm(neighbor, N)) continue; // pular tempestades

                if (!processed[neighbor] && dist[u] + weight < dist[neighbor]) {
                    dist[neighbor] = dist[u] + weight;
                    parent[neighbor] = u;

                    if (inHeap[neighbor]) {
                        pq.decrease_key(neighbor);
                    } else {
                        pq.insert(neighbor);
                        inHeap[neighbor] = true;
                    }
                }
            }
        } catch (const std::exception &e) {
            std::cout << "Erro ao acessar vizinhos do vertice " << std::endl;
        }
    }
    return {dist, parent};
}

//verifica se o exercito vai encontra um aliado ou um inimigo
bool ArmysAttack::tryEnemy(EstadoExercito &estado, int next_position, std::vector<EstadoExercito> &stateArmys) {
    for (auto &other: stateArmys) {
        if (other.position == next_position && !other.arrived) {
            //se ele encontrar ele perde a rodada
            if (!isAlly(estado, other.color)) {
                estado.stopped = true;
                return true; // encontrou inimigo
            } else {
                //se for amiguinho forma uma alianca
                estado.allys.push_back(other.color);
                other.allys.push_back(estado.color);
            }
        }
    }
    return false;
}

resultPath ArmysAttack::bestPathToCastle(Graph_Board &game_board,
                                               int at_beggining,
                                               const std::vector<int> &list_threats) {
    auto [dist, parent] = Dijkstra(game_board, at_beggining);

    int less_distant;
    int destino_final = chooseDestiny(list_threats, dist, less_distant);

    if (destino_final != -1) {
        auto caminho = reconstructPath(destino_final, parent);

        return {less_distant, caminho};
    }

    return {-1, {}};
}

//solucionar executa o calculo para todos os exercitos
void ArmysAttack::solve_game(Graph_Board &game_board) {
    int N = game_board.get_size();
    int castel_vertex = ChessNotToPos(castelPosition).first * N + ChessNotToPos(castelPosition).second;

    auto estados = startArmys(game_board, N, castel_vertex);

    int rodada = 0;
    int menor_movimentos_chegada = 1e9;

    while (rodada < 50) {
        rodada++;
        bool continuar = executeRound(estados, game_board, rodada, castel_vertex, menor_movimentos_chegada);
        if (!continuar) break;
    }

    auto results = processResult(estados, menor_movimentos_chegada);

    //imprime a lista de resultados finais
    for (auto &results: results) {
        std::cout << results.color << " " << results.moviments << " " << results.weight << " ";
    }
}

int main() {
    int line = 0;
    std::cin >> line; //linhas e colunas sao a mesma quantidade D:

    //um "catch" para caso a quantidade de linhas/colunas nao respeitem as regras
    if (line < 8 || line > 15) {
        throw std::invalid_argument("O numero de colunas linhas deve ser maior que 8 e menor que 15");
    }

    int armyNumber;
    std::cin >> armyNumber;
    std::cin.ignore();

    std::vector<Army> listArmys;

    //loop para adicionar as cores do exercito, inimigos, local de inicio etc no listaExercitos
    for (int i = 0; i < armyNumber; i++) {
        std::string line;
        std::getline(std::cin, line); // le a linha inteira

        std::istringstream iss(line);
        Army army;
        iss >> army.color >> army.position;

        std::string aliado;
        while (iss >> aliado) {
            army.friends.push_back(aliado);
        }

        listArmys.push_back(army);
    }

    std::string castelPosition;
    std::cin >> castelPosition;

    //adiciona as tormentas a uma lista
    int stormsNumber;
    std::cin >> stormsNumber;
    std::vector<std::string> stormsPositions(stormsNumber);
    for (int i = 0; i < stormsNumber; ++i) {
        std::cin >> stormsPositions[i];
    }

    //cria tabuleiro e executa a simulacao
    Graph_Board tab(line);
    ArmysAttack jogo(listArmys, castelPosition, stormsPositions);
    jogo.solve_game(tab);


    return 0;
}
/*
                                                     .. ..---:..... .....
                                                    ....=%@@@@@*-........
                                                     ..-%@@@%@@@@#+*##*-... .
                                                    ...=%@@*:-*@@@@@@@@@%=...
                                                    ....*@@#-::=%@@@##@@@@+..
                                                    ....:*@@%-::-**+-:=%@@@=.
                                                    .....-#@@#:::::::::=@@@%-
                                                        ..=@@@#:::::::::*@@@=. ..
                                                        ...+@@@+::::::::-%@@#:. .
            .. ..... ............                        ..:#@@@=::::::::+@@@=....
            .:+#%#+:..-+##%%%#=.                         ...=@@@*-:::::::-@@@+....
        ...:*@@@@@@%#@@@@@@@@@@=.                       . ...*@@%-::::::::%@@#:...
         .=@@@@**@@@@@@%%#*%@@%-.                       . . .-@@@+::::::::%@@%-...
        .-%@@@#:-*@@@%=::=#@@@=..                         ....%@@#-::::::-%@@%-...
        .*@@@#-:::-=-:::+%@@%+...                         ....#@@@=:::::::%@@%-...
       .:%@@@=:::::::::=%@@%=...                          ....*@@@*:::::::*@@%-...
     ...-@@@#-::::::::-#@@%=.   .                         ....+@@@#-::::::+%@%-...
    ....-@@@+:::::::::+@@@*... .                          ....+@@@#-::::::*@@#:...
    ....-%@@+:::::::::*@@%=.                              ....*@@@#-:::::=%@@*: .
      ..:#@@#-:::::::-%@@#-..                             ...:#@@@*-:::::+@@@+...
        .*@@@+:::::::-%@@#:...                          . .. -@@@@+:::::-#@@%-..
        .=@@@#-::::::-#@@%-...                           .. .*@@@*-:::::+%@@*...
       ...#@@@*:::::::*@@%=.                            ....=@@@%=:::::-#@@@=...
       ...-%@@%=::::::+@@@*..                           ...:#@@@+::::::=@@@%-...
       ....=@@@#=:::::-%@@%-......  ..... ... .... . ......*@@@#-:::::-#@@@+....
       .....+@@@#=:::::+@@@#-::-=++++++++++===--:::.......=@@@%=::::::+@@@#:.
        ....:#@@@%=:::::*@@@@@@@@@@@@@@@@@@@@@@@@@@@%#*=-+%@@@+::::::=%@@#:..
            .:#@@@@*-::::-*#%%%@@@%%%%%%%@@@@@@@@@@@@@@@@@@@@#-::::-*@@@*:...
            ...*@@@@+::::::::::--::::::::::::-----==++*#%%@@@@@#=:=#@@@*.....
            ..:#@@@#=:::::::::::::::::::::::::::::::::::::-=#@@@@%@@@@*..
            .-#@@@*-:::::::::::::::::::::::::::::::::::::::::=#@@@@@%=...
        ....=%@@@*:::::::::::::::::::::::::::::::::::::::::::::+@@@@%=...
          .-%@@@*:::::::::::::::::::::::::::::::::::::::::::::::-#@@@@=..
        ..-%@@@*::::::::::::::::::::::::::::::::::::::::::::::::::*@@@%-.....
        .-%@@@*::::::::::::::::::::::::::::::::::::::::::::::::::::*@@@%:...
    .. .-%@@@*-:::::::::::::::::::::::::::::::::::::::::::::::::::::*@@@*....
     ..-#@@@*-::::::::::::----::::::::::::::::::::::::-=======-:::::-#@@@+...
    ..:*@@@%-::::::=*%@@@@@@@@@@%##+:::::::::::::::+%@@@@@@@@@@@@%=::+@@@@-..
    ..+@@@@+::::::*@@@@@@@@@@@@@@@@%=:::::::::::::-%@@@@@@@@@@@@@@%=.-*@@@#:
    .-%@@@#-:::::=@@@@@@@@@@@@@@@@@@*:::::::::::::-%@@@@@@@@@@@@@@%-::-%@@@=.
....:*@@@@=::::::=@@@@@@@@@@@@@@@@@@+:::::::::::::-#@@@@@@@@@@@@@@+::::=@@@#: ...
. ..=%@@@#::::::::+@@@@@@@@@@@@@@@@#-::::::::::::::+@@@@@@@@@@@@@#-:::::+@@@+....
 ...*@@@%=:::::::::=%@@@@@@@@@@@@@%=::::::::::::::::+@@@@@@@@@@@+-::::::=%@@%- ..
...-%@@@#:::::::-==::=%@@@@@@@@@@#-::::::::::::::::::=%@@@@%#+==+*+::::::+@@@+...
. .*@@@@+:::::::-#@@@@@@@@@@@@@%=::::::::::::::::::::::=*%@@@@@@@%+::::::=%@@#-.
..:#@@@%=:::::::::=#%@@@@@@@%*=:::::::::::::::::::::::::::-=++*+-:::::::::#@@@=..
..-%@@@#-:::::::::::::---:::::::::::::::::::::::::::::::::::::::::::::::::*@@@+..
..=%@@@#::::::::::::::::::::::::::::=+**##%%%###**+-::::::::::::::::::::::=@@@+..
..=@@@@*::::::::::::::::::::::::::+%@@@@@@@@@@@@@@@@@+-:::::::::::::::::::=%@@%-.
..-%@@@#:::::::::::::::::::::::::+%@@@@%%#*++++*#%@@@@*-::::::::::::::::::=@@@@@%-..
...+@@@%+::::::::::::::::::::::::#@@@#-:::::::::::-%@@@+::::::::::::::::::=@@@@@@@=.
...:#@@@@*-:::::::::::::::::::::-%@@%=:::::::::::::=%@@#=::::::::::::::::-*@@@@@@@%:
. ..:+%@@@@*=:::::::::::::::::::=@@@*-::::::::::::::*@@@+::::::::::::::-=#@@@@@@@@*.
    ...+%@@@@@%##**+===------::-+@@@*-:::::::-:::-:-*@@@*------====++*%@@@@@@@@@*-..
     ...:+#%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*-...
          ....:-=+**##%%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@%#**%@@@%+:. .  .
            ..................=%@@%##%@@####%@@%##%%%%@@@@@@@%@@%*++*#@@@@#-.....
             .................*@@@%#%@@%#*##%@@%***###@@@@%@@@@#+++*%@@@@%:......
                            .:*@@@%%@@@##*#%@@@#*#*###@@@%%@@@#++#%@##%@@@*:.
                            .-%@@%%@@@%**##@@@%#*####%@@@%#@@@%%@@%*+*%@@@@=.
                         ....+@@@@@@@%#**#%@@@#***###@@@@##%@@@@%#*#%@@@@#-..
                        . ..-%@@@@@@%####%@@@%##*#*#%@@@%#*#@@@@##@@@@@*:.  .
                        . .:*@@@@@@%#*##@@@@%#*#**#%@@@@%#*#%@@@@@@@%+:..
                        ...+@@@@@@###*#@@@@%#*##*##@@@@@@@#*#@@@@%*=....
                     .....*@@@@@%#*##%@@@@%#*#####@@@@%@@@%*#%@@@#: .
                    ... .*@@@@@##**#%@@@@%#*##**#@@@@@#%@@@#*#@@@@=..
                    ...:#@@@@##***#@@@@@###**###%@@@@#*#%@@%*#%@@@*.....
                    ..=@@@@%#*###%@@@@%#*#**#*#@@@@@###*#@@@%##@@@#:....
                    :*@@@@#*#**#@@@@##*##**##%@@@@@##***#%@@%#*@@@@-....
                .. :*@@@%#####%@@@%#*##**#*#%@@@@@@##*#**#%@@##%@@@=....
                ...+@@@%#*#*#@@@@%######*##@@@@@@@@%**##*#%@@%##@@@=....
                ..:@@@%#####@@@@##*##**###@@@@@##@@@#*#***%@@@#*@@@=.
                ..-@@@#***%@@@%#**#*#####@@@@%#*#%@@%####*%@@@#*%@@+.
                .:#@@##*#%@@@%#*#**##*##@@@@##***#@@@%###*#@@@%#%@@*.
                .=@@@###%@@@%#*#****##@@@@@#*######@@@@##*#%@@@#%@@*:
                .*@@@##@@@@%#*****#%@@@@@@@%#*##*#*#@@@@%#*#@@@@@@@+..
                :#@@@#%@@@%#####%@@@@@@@@@@@@%#**####@@@@%#*#@@@@@%-..
                :*@@@%@@@%#*#%@@@@@@@@@@@@@@@@@%#**###%@@@@%#%@@@@*...
               ..=@@@@@@%%%@@@@@@@@@@@@@@@@@@@@@@@%%##*#%@@@@%%@@@#:.
               ..:%@@@@@@@@@@@@@@@@@#+====+++%@@@@@@@@@%%%%@@@@@@@@*:
                .:*@@@@@@@@%%@@@@@@+..... ....*@@@@@@@@@@@@@@@@@@@@@*....
                ..=@@@@%*-:.:%@@@@*:...  .....*@@@@@%%@@@@@@@@@@@@@#=....
                ...=+-:..   .-%@@#- .        .*@@@%=...-=+*****==-:..
                     ....   . ......        ...... ..................    */