#include <iostream>
#include <list>
#include <vector>


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
    void sift_up(int i) {
        while (i > 0 && heap[i] < heap[parent(i)]) {
            int temp = heap[i];
            heap[i] = heap[parent(i)];
            heap[parent(i)] = temp;
            i = parent(i);
        }
    }

    //quando o minimo e removido, fica um buraco no árvore do heap
    //o sift_down "desce" o menor número e faz dnv a comparação para verificar novamente
    void sift_down(int i) {
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
                int temp = heap[i];
                heap[i] = heap[smallest];
                heap[smallest] = temp;
                i = smallest;
            } else break;
        }
    }

public:
    bool empty() const { return heap.empty(); }
    int size() const { return heap.size(); }
    int minimum(std::vector<T> heap) const;
    void insert(T element);
    int extract_min();
};

//retorna mínimo
template<typename T>
int PrioridadeMinima<T>::minimum(std::vector<T> heap) const{
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
int PrioridadeMinima<T>::extract_min() {
    if (heap.empty()) {
        throw std::runtime_error("Heap vazior");
    }
    int min = heap[0];
    heap[0] = heap.back();
    heap.pop_back();
    if (!heap.empty()) {
        sift_down(0);
    }

    return min;
}

//tem que ajeitar as pouco essa classe pro PP3
// -x-x-x-x- Classe Tabuleiro baseada num grafo lista de adjacencias  -x-x-x-x-
/*class Tabuleiro
{
private:
    uint num_vertices;
    uint num_edges;
    std::vector<Vertex> *adj;
    std::vector<std::pair<int, int>> movimentos = {
        {1, 2}, {1, -2}, {-1, 2}, {-1, -2}, {2, 1}, {2, -1}, {-2, 1}, {-2, -1}};

    void criar_Grafo();
    void insertionSort(std::vector<Vertex> &vec);

public:
    Tabuleiro();
    ~Tabuleiro();
    void add_edge(Vertex u, Vertex v);
    std::vector<Vertex> get_adj(Vertex u);

    uint get_edges();
    uint get_vertices();
};

Tabuleiro::Tabuleiro()
{
    num_vertices = 64;
    num_edges = 0;
    adj = new std::vector<Vertex>[num_vertices];
    criar_Grafo();
}

Tabuleiro::~Tabuleiro()
{
    delete[] adj;
    adj = nullptr;
}

void Tabuleiro::add_edge(Vertex u, Vertex v)
{
    if (u >= num_vertices || v >= num_vertices || u == v)
    {
        throw std::invalid_argument("Valores invalidos");
    }

    adj[u].push_back(v);
    adj[v].push_back(u);
    num_edges += 1;
}

std::vector<Vertex> Tabuleiro::get_adj(Vertex u)
{
    if (u >= num_vertices)
    {
        throw std::invalid_argument("Valores invalidos");
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

void Tabuleiro::criar_Grafo()
{
    for (int linha = 0; linha < 8; ++linha)
    {
        for (int coluna = 0; coluna < 8; ++coluna)
        {
            int vertice_origem = linha * 8 + coluna;
            for (const auto &mov : movimentos)
            {
                int nova_linha = linha + mov.first;
                int nova_coluna = coluna + mov.second;
                if (nova_linha >= 0 && nova_linha < 8 && nova_coluna >= 0 && nova_coluna < 8)
                {

                    int vertice_destino = nova_linha * 8 + nova_coluna;


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
        insertionSort(adj[u]);
    }
}

void Tabuleiro::insertionSort(std::vector<Vertex> &vec)
{
    if (vec.size() < 2)
    {
        return;
    }

    for (Vertex i = 1; i < vec.size(); i++)
    {
        Vertex key = vec[i];
        int j = i - 1;

        while (j >= 0 && vec[j] > key)
        {
            vec[j + 1] = vec[j];
            j = j - 1;
        }
        vec[j + 1] = key;
    }
}

class AtaqueDosCavaleiros
{
private:
    std::vector<std::string> cavaleiros;
    std::string reiPosicao;
    std::vector<int> ameacasRei();
    int melhorCaminhoAoRei(Tabuleiro &tabubu, int no_inicio, const std::vector<int> &listaAmeacas);

    std::pair<int, int> ChessNotToPos(const std::string &position);

public:
    AtaqueDosCavaleiros(std::vector<std::string> knights, std::string king);
    void solucionar( Tabuleiro &tabubu);
};

AtaqueDosCavaleiros::AtaqueDosCavaleiros(std::vector<std::string> knights, std::string king)
{
    this->reiPosicao = king;
    this->cavaleiros = knights;
}

std::pair<int, int> AtaqueDosCavaleiros::ChessNotToPos(const std::string &position)
{
    return {position[1] - '1', position[0] - 'a'};
}

std::vector<int> AtaqueDosCavaleiros::ameacasRei()
{
    const std::vector<std::pair<int, int>> movimentos = {
        {2, 1}, {2, -1}, {-2, 1}, {-2, -1}, {1, 2}, {1, -2}, {-1, 2}, {-1, -2}};

    int coluna = reiPosicao[0] - 'a';
    int linha = reiPosicao[1] - '1';

    std::vector<int> listaAmeacas;
    for (auto mov : movimentos)
    {
        int nova_linha = linha + mov.first;
        int nova_coluna = coluna + mov.second;
        if (nova_linha >= 0 && nova_linha < 8 && nova_coluna >= 0 && nova_coluna < 8)
        {
            listaAmeacas.push_back(nova_linha * 8 + nova_coluna);
        }
    }

    return listaAmeacas;
}



void AtaqueDosCavaleiros::solucionar(Tabuleiro &tabubu)
{
    std::vector<int> posicao_alvo = ameacasRei();
    std::vector<int> resultados;

    int minimo_mov = 3000;

    for (const auto &cavalo : this->cavaleiros)
    {
        std::pair<int, int> cava_pos = ChessNotToPos(cavalo);
        Vertex no_inicio = cava_pos.first * 8 + cava_pos.second;
        int move = melhorCaminhoAoRei(tabubu, no_inicio, posicao_alvo);
        resultados.push_back(move);
        if (move != -1 && move < minimo_mov)
        {
            minimo_mov = move;
        }
    }


    bool first = true;
    for (int move_count : resultados)
    {
        if (move_count == minimo_mov)
        {
            if (!first)
                std::cout << " ";
            std::cout << move_count;
            first = false;
        }
    }
    std::cout << std::endl;
}*/

int main()
{
    int line = 0;
    std::string army;
    std::vector <std::string> listArmys;

    //linhas e colunas são a mesma quantidade D:
    std::cin >> line;

    //um "catch" para caso a quantidade de linhas/colunas não respeitem as regras
    if (line < 8 || line > 15) {
        throw std::invalid_argument("O número de colunas linhas deve ser maior que 8 e menor que 15");
    }

    //loop para adicionar as cores do exercito, inimigos, local de início etc no listaExercitos
    for (int i=0; i < 5; i++) {
        std::cin >> army;
        listArmys.push_back(army);
    }

    return 0;
}