#include <iostream>
#include <math.h>
#include <map>
#include <list>
#include <omp.h>
#include <fstream>
#include <queue>
#include <time.h>
#include <windows.h>

int THREADS = 4;
int GRAINSIZE = 128;
int BAGSIZE = 256;
using namespace std;
//nda pennant
class Node{
	friend class Bag;
	private:
		int value = -1; //representa o número do vertice
		int distance = -1; //representa a distancia do vertice
		Node* left; //ponteiro para o filho a esquerda
		Node* right;//ponteiro para o filho a direita
		list<Node*> adjacent; //lista de adjacencia do vertice
		int pennantSize(Node* root, int tamanho);//metodo que retorna o tamanho da pennant 
		void pennantToVector(list<Node*>& vector,Node* pennant, int position);//metodo que converte a pennant em uma lista
	public:
		Node(int value);//construtor
		Node();//construtor
		void insertAdjacent(map<int, Node*> graph, int i);//metodo responsavel por inserir o vizinho no vetor de adjacencia
		void printValue();//metodo auxiliar para imprimir a posicao do vetor
		void printGraph(map<int, Node*> graph);//metodo para imprimir o grafo
};

//construtor do no com dado
Node::Node(int data) {
	value = data;
	left = NULL;
	right = NULL;
}

//construtor de no vazio
Node::Node() {
	value = -1;
	left = NULL;
	right = NULL;
}

//coloca os vertices da pennant em ordem em uma lista para percorrer na funcao processPennant
void Node::pennantToVector(list<Node*> &pennant,Node* inPennant, int position){
	if (inPennant != NULL) {
		pennant.push_back(inPennant);
	pennantToVector(pennant, inPennant->left, position++);
	pennantToVector(pennant, inPennant->right, position++);
	}
	
}

//insere nos adjacentes, insere lista de adjacencia(grafo)
void Node::insertAdjacent(map<int, Node*> graph, int i){
	this->adjacent.push_back(graph[i]);
	this->adjacent.unique();
	}	



//imprime grafo
void Node::printGraph(map<int, Node*> graph) {
	cout << "GRAPH" << endl;
	for (int i = 0; i < graph.size(); i++){
		cout << i;
		for (Node* x : graph[i]->adjacent) {
			cout << "  -> " << x->value << "(" << x->distance << ")" << " ";
		}
		cout << endl;
	}
}


//metodo aux para impressao(Debug)
void Node::printValue() {
	cout << this->value << endl;
}

//metodo recursivo para contar o tamanho da pennant
int Node::pennantSize(Node* node, int tamanho){
	if (node != NULL) {
		pennantSize(node->left, tamanho + 1);
		pennantSize(node->right, tamanho + 1);
		return tamanho + 1;
	}
}

//classe bag
class Bag {
	private:
		Node** vector; //vetor de nos 
		void processLayer(Bag* inBag, Bag* outBag, int d); //metodo que processa cada camda do grafo
		void processPenant(Node* inPennant, Bag* outBag,int d); //metodo que processa cada pennant
		Node* pennantUnion(Node* x, Node* y); //realiza a uniao de pennants
		Node* pennantSplit(Node* x); //separa duas pennants
		int bagSize();//retorna o tamanho da bag 
		int elementsInBag = 0;
	public:
		int size;
		Bag(int graphSize); //construtor
		void insertBag(Node* x);//insere na bag
		void debug(); //percorre a bag e imprime os valores
		void percorre(Node* node);//percorre a pennant em ordem
		void PBFS(map<int, Node*> graph); //busca em largura com bags
		void serialBFS(map<int, Node*> graph); //busca em largura com fila
};

//metodo para fazer busca em largura
void Bag::PBFS(map<int, Node*> graph){
	//calcula tempo de execucao
	LARGE_INTEGER clockFrequency;
	QueryPerformanceFrequency(&clockFrequency);
	LARGE_INTEGER start_time;
	LARGE_INTEGER end_time;
	QueryPerformanceCounter(&start_time);

	graph[0]->distance = 0; // coloca a distancia do vertice inical em 0
	int d = 0; 
	Bag* V0 = new Bag(BAGSIZE);	//cria nova bag(primeira inbag)
	V0->insertBag(graph[0]);
	map<int, Bag*> vectorBags; // cria map de bags 
	vectorBags.insert(pair<int, Bag*>(0, V0));	//insere a primeira bag
	while (vectorBags[d]->elementsInBag > 0) {
		/*
		cout << "Nivel: " << d << endl;
		vectorBags[d]->debug();
		cout << endl;
		*/
		vectorBags.insert(pair<int, Bag*>(d + 1, new Bag(BAGSIZE)));
		processLayer(vectorBags[d], vectorBags[d + 1], d);
		d++;
	}
	QueryPerformanceCounter(&end_time);
	LARGE_INTEGER delta;
	delta.QuadPart = end_time.QuadPart - start_time.QuadPart;
	double deltaSeconds = ((double)delta.QuadPart) / clockFrequency.QuadPart;
	cout << "Tempo de execucao Busca em Largura paralela: " << deltaSeconds << endl;
}


//metodo process layer da busca
void Bag::processLayer(Bag* inBag, Bag* outBag, int d) {
	#pragma omp parallel for num_threads(THREADS)
	for (int k = 0; k < 9; k++) {
		if (inBag->vector[k] != NULL) {
			processPenant(inBag->vector[k], outBag, d);
		}
		//#pragma omp critical
	//	cout << "NumeroThread: " << omp_get_thread_num() << endl;
	}
}


void Bag::processPenant(Node* inPennant, Bag* outBag, int d){
	list<Node*>pennant;
	pennant.push_back(inPennant);
	inPennant->pennantToVector(pennant, inPennant->left, 0);
	if (pennant.size() < GRAINSIZE+1){		
		for (int i = pennant.size(); i > 0 ; i--)
		{
			//Pega primeiro elemento da lista(pennant)
			inPennant = pennant.front();

			//Remove primeiro elemento da lista(pennant)
			pennant.pop_front();

			//Para cada adjacente insere na outBag
			#pragma omp parallel for num_threads(THREADS)
			for(int j = 0; j < inPennant->adjacent.size(); j++){
				Node* x = inPennant->adjacent.front();
				inPennant->adjacent.pop_front();
				inPennant->adjacent.push_back(x);
				#pragma omp critical
				if (x->distance == -1) {
					x->distance = d+1;
					outBag->insertBag(x);
				}
		//#pragma omp critical
		//	cout << "NumeroThread: " << omp_get_thread_num() << endl;
			
			}
		}
	}
	else {
		Node* newPennant = pennantSplit(inPennant);
		#pragma omp task
		processPenant(newPennant, outBag, d);
		processPenant(inPennant, outBag, d);
		#pragma omp taskwait
	}
}

//retorna o tamanho da bag
int Bag::bagSize() {
	int maior = 0;
	for (int i = 0; i < 4; i++) {
		if (this->vector[i] != NULL) {
			maior = i;
		}
	}
	return maior + 1;
}


//construtor da bag
Bag::Bag(int graphSize){
	if (graphSize > 0){
		size = log2(graphSize)+1;
		vector = new Node * [size];
		for (int i = 0; i < size; i++){
			vector[i] = NULL;
		}
	}
	else {
		cerr << "Grafo deve ter pelo menos um n�" << endl;
	}	
}

//insere na bag
void Bag::insertBag(Node* x){
	int k = 0;
	while (vector[k] != NULL){
		x = pennantUnion(vector[k], x);
		vector[k] = NULL;
		k++;
	}
	vector[k] = x;
	this->elementsInBag++;
}

//split pennant
Node* Bag::pennantSplit(Node* x){
	Node* y = new Node();
	y = x->left;
	x->left = y->right;
	y->right = NULL;
	return y;
}

//pennant union
Node* Bag::pennantUnion(Node* x, Node* y) {
	y->right = x->left;
	x->left = y;
	return x;
}

//debug
void Bag::debug() {
	for (int i = 0; i < size; i++)
	{
		if (this->vector[i] == NULL){
			cout << i << ":null" << endl;
		}
		else {
			cout << i << ": " << vector[i]->value << " ";
			percorre(vector[i]->left);
			cout << endl;
		}
	}
}


//percorre a bag em ordem
void Bag::percorre(Node* node) {
	if(node != NULL){
		cout << node->value << " ";
		percorre(node->left);
		percorre(node->right);
	}
	
}

//busca em largura serial com fila
void Bag::serialBFS(map<int, Node*>graph){
	LARGE_INTEGER clockFrequency;
	QueryPerformanceFrequency(&clockFrequency);
	LARGE_INTEGER start_time;
	LARGE_INTEGER end_time;
	QueryPerformanceCounter(&start_time);
	queue<Node*> fila;
	graph[0]->distance = 0;
	fila.push(graph[0]);
	while(fila.size()>0){
		Node* vertice = fila.front();
		fila.pop();
		for (Node* x : vertice->adjacent) {
			if (x->distance == -1) {
				x->distance = vertice->distance + 1;
				fila.push(x);
			}
		}
	}
	QueryPerformanceCounter(&end_time);
	LARGE_INTEGER delta;
	delta.QuadPart = end_time.QuadPart - start_time.QuadPart;
	double deltaSeconds = ((double)delta.QuadPart) / clockFrequency.QuadPart;
	cout << "Tempo de execucao Busca em Largura com fila: " << deltaSeconds << endl;
}

int main() {
	int graphSize;
	int edges;
	Bag* bag = new Bag(BAGSIZE);
	ifstream grafFile;
	grafFile.open("Teste.txt");
	grafFile >> graphSize >> edges;
	map<int, Node*> graph;
	for (int i = 0; i < graphSize; i++){
		graph.insert(pair<int, Node*>(i, new Node(i)));
	}

	for (int i = 0; i < edges; i++){
		int insert, adjacent;
		grafFile >> insert >> adjacent;
		Node* node = graph[insert-1];
		node->insertAdjacent(graph,adjacent-1);
		node = graph[adjacent - 1];
		node->insertAdjacent(graph, insert - 1);
	}

	//graph[0]->printGraph(graph);
	//bag->serialBFS(graph);
	bag->PBFS(graph);
	//graph[0]->printGraph(graph);
	
}