#include <iostream>
#include <list>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <sstream>

using Vertex = unsigned int;
using uint = unsigned int;

template <typename T>
class PrioridadeMinima
{
private:
    std::vector<T> heap;

    int parent(int i) { return (i - 1) / 2; }
    int left(int i)   { return 2 * i + 1; }
    int right(int i)  { return 2 * i + 2; }

    //"subir" na árvore do heap se não atender aos requesitos da árvore
    //caso o filho de cima seja menor que o de baixo etc
    void sift_up(int i)
    {
        while (i > 0 && heap[i] < heap[parent(i)]) {
            T temp = heap[i];
            heap[i] = heap[parent(i)];
            heap[parent(i)] = temp;
            i = parent(i);
        }
    }

    //quando o minimo e removido, fica um buraco no árvore do heap
    //o sift_down "desce" o menor número e faz dnv a comparação para verificar novamente
    void sift_down(int i)
    {
        int n = heap.size();
        while (true) {
            int smallest = i;
            if (left(i) < n && heap[left(i)] < heap[smallest]) {
                smallest = left(i);
            }
            if (right(i) < n && heap[right(i)] < heap[smallest]) {
                smallest = right(i);
            }
            if (smallest != i) {
                T temp = heap[i];
                heap[i] = heap[smallest];
                heap[smallest] = temp;
                i = smallest;
            } else break;
        }
    }

public:
    bool empty() const { return heap.empty(); }
    int size() const { return heap.size(); }
    T minimum(std::vector<T> heap) const;
    void insert(T element);
    T extract_min();
};

//retorna mínimo
template<typename T>
T PrioridadeMinima<T>::minimum(std::vector<T> heap) const{
    if (heap.empty()) {
        throw std::runtime_error("Heap vazior");
    }
    return heap[0];
}

//insere novo elemento
template<typename T>
void PrioridadeMinima<T>::insert(T element) {
    heap.push_back(element);
    sift_up(heap.size() - 1);
}

// remove e retorna o mínimo
template<typename T>
T PrioridadeMinima<T>::extract_min() {
    if (heap.empty()) {
        throw std::runtime_error("Heap vazior");
    }
    T min = heap[0];
    heap[0] = heap.back();
    heap.pop_back();
    if (!heap.empty()) {
        sift_down(0);
    }

    return min;
}


//tem que ajeitar as pouco essa classe pro PP3
// -x-x-x-x- Classe Tabuleiro baseada num grafo lista de adjacencias  -x-x-x-x-
class Tabuleiro
{
private:
    uint num_vertices;
    uint num_edges;
    uint N; //dimensão do tabuleiro

    std::vector<std::pair<Vertex,int>> *adj; //aq virou um par pois guarda o adk e o seu peso
    std::vector<std::pair<int, int>> movimentos = {
        {1, 2}, {1, -2}, {-1, 2}, {-1, -2}, {2, 1}, {2, -1}, {-2, 1}, {-2, -1}};

    void criarGrafo();
    void insertionSort(std::vector<std::pair<Vertex, int>> &vec);


public:
    Tabuleiro(uint N = 8);
    ~Tabuleiro();

    void add_edge(Vertex u, Vertex v);
    uint get_edges();
    const std::vector<std::pair<Vertex,int>>& get_adj(Vertex u) const;
    uint get_vertices();
    uint get_size() const { return N;};
    int calcularPeso(int u, int v, int N);
};

//construtor do Grafo tabuleiro
Tabuleiro::Tabuleiro(uint N)
{
    this->N = N;
    num_vertices = N * N;
    num_edges = 0;
    adj = new std::vector<std::pair<Vertex,int>>[num_vertices];
    std::cout << "O tabuleiro foi criado com sucessor!" << std::endl;

    criarGrafo();
}

//destrutor do Grafo tabuleiro
Tabuleiro::~Tabuleiro()
{
    delete[] adj;
    adj = nullptr;
}

//calculadora de peso q o prof pediu
//w(u,v)=(ascii(αu)⋅βu+ascii(αv)⋅βv)mod19 (fórmula)
int Tabuleiro::calcularPeso(int u, int v, int N) {
    int linha_u = u / N;
    int coluna_u = u % N;

    char alpha_u = 'a' + coluna_u; //au
    int beta_u = linha_u + 1; //bu

    int linha_v = v / N;
    int coluna_v = v % N;

    char alpha_v = 'a' + coluna_v; //av
    int beta_v = linha_v + 1; //bv

    std::cout << "O calcular peso ta funcionando eu achor" << std::endl;
    //esse int(letra) = ascii(letra)
    return ((int(alpha_u) * beta_u + int(alpha_v) * beta_v) % 19);
}

void Tabuleiro::add_edge(Vertex u, Vertex v)
{
    if (u >= num_vertices || v >= num_vertices || u == v)
    {
        throw std::invalid_argument("Valores invalidos no add_edge");
    }

    int weight = calcularPeso(u, v, N);

    //aq o adj cria uma nova aresta
    adj[u].push_back({v,weight});
    adj[v].push_back({u, weight});
    num_edges += 1;

    std::cout << "O add_edge ta de boas tbm eu acho" << std::endl;
}

const std::vector<std::pair<Vertex,int>>& Tabuleiro::get_adj(Vertex u) const
{
    if (u >= num_vertices)
    {
        throw std::invalid_argument("Valores invalidos no get_adj");
    }
    return adj[u];
}

uint Tabuleiro::get_edges()
{
    return num_edges;
}

uint Tabuleiro::get_vertices()
{
    return num_vertices;
}

//como agr o tamanho do tabuleiro pode mudar de tamanho coloca-se N ao inves do 8 anterior
void Tabuleiro::criarGrafo()
{
    std::cout << "O grafo esta sendo criador no inicior!" << std::endl;
    for (int linha = 0; linha < N; ++linha)
    {
        for (int coluna = 0; coluna < N; ++coluna)
        {
            int vertice_origem = linha * N + coluna;
            for (const auto &mov : movimentos)
            {
                int nova_linha = linha + mov.first;
                int nova_coluna = coluna + mov.second;
                if (nova_linha >= 0 && nova_linha < N && nova_coluna >= 0 && nova_coluna < N)
                {
                    int vertice_destino = nova_linha * N + nova_coluna;

                    if (vertice_origem < vertice_destino)
                    {
                        add_edge(vertice_origem, vertice_destino);
                    }
                }
            }
        }
    }

    for (uint u = 0; u < num_vertices; ++u)
    {
        std::cout << "O vertice " << u << " esta sendo adicionador" << std::endl;
        insertionSort(adj[u]);
    }
}

void Tabuleiro::insertionSort(std::vector<std::pair<Vertex,int>> &vec)
{
    if (vec.size() < 2)
    {
        return;
    }

    for (int i = 1; i < vec.size(); i++)
    {
        auto key = vec[i];
        int j = i - 1;

        while (j >= 0 && vec[j].first > key.first)
        {
            vec[j + 1] = vec[j];
            j = j - 1;
        }
        vec[j + 1] = key;
    }
}

// ----------------- Estrutura auxiliar para guardar o caminho -----------------
// ta aqui a struct que guarda a distância e o caminho percorrido
struct ResultadoCaminho {
    int distancia;
    std::vector<int> caminho;

    //construtor pq da erro la na frente por causa do pair
    ResultadoCaminho(int d, std::vector<int> c)
        : distancia(d), caminho(std::move(c)) {}
};

struct Army {
    std::string color;
    std::string position;
    std::vector<std::string> friends;
};

// -x-x-x-x- Classe AtaqueDosExercitos  -x-x-x-x-
class AtaqueDosExercitos
{
private:
    std::vector<Army> armys;
    std::string castelPosition;
    std::vector<std::string> storms;

    std::vector<int> castleChess(Tabuleiro tabuleiro);
    bool isStorm(int v, int N);
    ResultadoCaminho melhorCaminhoAoRei(Tabuleiro &tabubu, int no_inicio, const std::vector<int> &listaAmeacas);

    std::pair<int, int> ChessNotToPos(const std::string &position);

public:
    AtaqueDosExercitos(std::vector<Army> listArmy, std::string castelPosition, std::vector<std::string> storms);
    void solucionar( Tabuleiro &tabubu);
};

AtaqueDosExercitos::AtaqueDosExercitos(std::vector<Army> listArmy, std::string castelPosition, std::vector<std::string> storms)
{
    this->castelPosition = castelPosition;
    this->armys = listArmy;
    this->storms = storms;
    std::cout << "Construtor dos exércitos criado" << std::endl;
}

//converte notação tipo "a1" em coordenadas (linha, coluna)
std::pair<int, int> AtaqueDosExercitos::ChessNotToPos(const std::string &position)
{
    if (position.size() < 2) {
        throw std::invalid_argument("Posição inválida");
    }

    //aqui precisa do stoi pq se for a10, a11 etc da errador
    int coluna = position[0] - 'a';
    int linha = std::stoi(position.substr(1)) - 1;
    return {linha, coluna};
}

//gera as posições que são vizinhas do castelo (um cavalo pode alcançar)
std::vector<int> AtaqueDosExercitos::castleChess(Tabuleiro tabuleiro)
{
    //movimentos dos exércitos em formato de L (fazuelli)
    const std::vector<std::pair<int, int>> movimentos = {
        {2, 1}, {2, -1}, {-2, 1}, {-2, -1}, {1, 2}, {1, -2}, {-1, 2}, {-1, -2}};

    int coluna = castelPosition[0] - 'a';
    int linha = castelPosition[1] - '1';

    //a listSus é uma lista de possíveis chegadas ao castelo
    std::vector<int> listSus;

    for (auto mov : movimentos)
    {
        int nova_linha = linha + mov.first;
        int nova_coluna = coluna + mov.second;
        if (nova_linha >= 0 && nova_linha < tabuleiro.get_size() && nova_coluna >= 0 && nova_coluna < 8)
        {
            listSus.push_back(nova_linha * tabuleiro.get_size() + nova_coluna);
        }
    }

    return listSus;
}

//verificação se tem Storm na posição do tabuleiro
bool AtaqueDosExercitos::isStorm(int v, int N) {
    int linha = v / N;
    int coluna = v % N;
    for (auto &s : storms) {
        int l = s[1] - '1';
        int c = s[0] - 'a';
        if (l == linha && c == coluna) return true;
    }
    return false;
}

// ----------------- Implementação do Dijkstra + reconstrução de caminho -----------------
ResultadoCaminho AtaqueDosExercitos::melhorCaminhoAoRei(Tabuleiro &tabubu, int no_inicio, const std::vector<int> &listaAmeacas) {
    int V = tabubu.get_vertices();
    std::vector<int> dist(V, 1e9); //1e9 é o infinito
    std::vector<int> parent(V, -1);

    // pair{distancia acumulada, vértice}
    PrioridadeMinima<std::pair<int,int>> pm;

    // Inicialização
    dist[no_inicio] = 0;
    pm.insert({0, no_inicio});

    // Enquanto houver vértices a processar
    while (!pm.empty()) {
        auto par = pm.extract_min();
        int d = par.first;
        int u = par.second;

        if (d > dist[u]) continue;

        if (d > dist[u]) continue;

        // Se o vértice atual é uma posição válida do castelo -> encontramos caminho
        if (std::find(listaAmeacas.begin(), listaAmeacas.end(), u) != listaAmeacas.end()) {
            // Reconstruir o caminho até 'u'
            std::vector<int> caminho;
            for (int v = u; v != -1; v = parent[v]) {
                caminho.push_back(v);
            }
            std::reverse(caminho.begin(), caminho.end());

            return {dist[u], caminho};
        }

        // Relaxamento das arestas
        for (auto [v, peso] : tabubu.get_adj(u)) {
            // Aqui você pode checar se 'v' é tormenta e ignorar
            if (dist[v] > dist[u] + peso) {
                dist[v] = dist[u] + peso;
                parent[v] = u;
                pm.insert({dist[v], v});
            }
        }
    }

    // Se não encontrou caminho até o castelo
    return {-1, {}};
}

//solucionar executa o cálculo para todos os exércitos
void AtaqueDosExercitos::solucionar(Tabuleiro &tabubu)
{
    std::cout << std::endl << "Solucionar está sendo chamado" << std::endl;

    std::vector<int> posicao_alvo = castleChess(tabubu);

    for (const auto &exercito : this->armys)
    {
        std::string cor = exercito.color;
        std::string posicao = exercito.position;

        std::pair<int, int> cava_pos = ChessNotToPos(posicao);
        Vertex no_inicio = cava_pos.first * tabubu.get_size() + cava_pos.second;

        ResultadoCaminho resultado = melhorCaminhoAoRei(tabubu, no_inicio, posicao_alvo);

        if (resultado.distancia != -1) {
            std::cout << "Exército " << cor << " em " << posicao
                  << " chega em " << resultado.distancia << " movimentos.\n";

            std::cout << "Caminho: ";
            for (auto v : resultado.caminho) {
                int l = v / tabubu.get_size() ;
                int c = v % tabubu.get_size();
                char colChar = 'a' + c;
                char rowChar = '1' + l;
                std::cout << colChar << rowChar << " ";
            }
            std::cout << "\n";
        } else {
            std::cout << "Exército " << cor << " em " << posicao
                  << " não encontrou caminho.\n";
        }
    }
}



int main()
{
    int line = 0;
    std::cin >> line; //linhas e colunas são a mesma quantidade D:

    //um "catch" para caso a quantidade de linhas/colunas não respeitem as regras
    if (line < 8 || line > 15) {
        throw std::invalid_argument("O número de colunas linhas deve ser maior que 8 e menor que 15");
    }

    int armyNumber;
    std::cin >> armyNumber;
    std::cin.ignore();

    std::vector<Army> listArmys;

    //loop para adicionar as cores do exercito, inimigos, local de início etc no listaExercitos
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



    //cria tabuleiro e executa a simulação
    Tabuleiro tab(line);
    AtaqueDosExercitos jogo(listArmys, castelPosition, stormsPositions);

    std::cout << "Tabuleiro " << line << "x" << line << std::endl;
    for (auto &e : listArmys) {
        std::cout << "Exercito: " << e.color << " posicao: " << e.position;
        if (!e.friends.empty()) {
            std::cout << " aliados: ";
            for (auto &a : e.friends)
                std::cout << a << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "Castelo em: " << castelPosition << std::endl;
    for (auto &s : stormsPositions)
        std::cout << "Storm em: " << s << std::endl;

    jogo.solucionar(tab);


    return 0;
}