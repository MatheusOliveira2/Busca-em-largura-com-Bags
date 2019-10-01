#include <iostream>
#include <math.h>
#include <map>
#include <list>
#include <omp.h>
using namespace std;


//nó da pennant
class Node{
	friend class Bag;
	private:
		int value = -1;
		int distance = -1; 
		Node* left;
		Node* right;
		list<Node*> adjacent;
		int pennantSize(Node* root, int tamanho);
		void pennantToVector(list<Node*> vector,Node* pennant, int position);
	public:
		Node(int value);
		Node();
		void insertAdjacent(map<int, Node*> graph, int i);
		void printValue();
		void printGraph(map<int, Node*> graph);
};

//construtor do nó com dado
Node::Node(int data) {
	value = data;
	left = NULL;
	right = NULL;
}

//construtor de nó vazio
Node::Node() {
	value = -1;
	left = NULL;
	right = NULL;
}

//coloca os vértices da pennant em ordem em uma lista para percorrer na função processPennant
void Node::pennantToVector(list<Node*> pennant,Node* inPennant, int position){
	if (inPennant != NULL) {
		pennant.push_back(inPennant);
	pennantToVector(pennant, inPennant->left, position++);
	pennantToVector(pennant, inPennant->right, position++);
	}
	
}

//insere nós adjacentes, insere lista de adjacencia(grafo)
void Node::insertAdjacent(map<int, Node*> graph, int i){
	int vertex = 0;
	cout << "Quais vertices ligam no vertice " << i << endl;
	while (vertex > -1) {
		cin >> vertex;
		if(vertex < graph.size())
			this->adjacent.push_back(graph[vertex]);
	}	
}


//imprime grafo
void Node::printGraph(map<int, Node*> graph) {
	cout << "GRAPH" << endl;
	for (int i = 0; i < graph.size(); i++){
		cout << i;
		for (Node* x : graph[i]->adjacent) {
			cout << "  -> " << x->value << " ";
		}
		cout << endl;
	}
}


//metodo aux para impressão(Debug)
void Node::printValue() {
	cout << this->value << endl;
}

//método recursivo para contar o tamanho da pennant
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
		Node** vector;
		void processLayer(Bag* inBag, Bag* outBag, int d);
		void processPenant(Node* inPennant, Bag* outBag,int d);
		Node* pennantUnion(Node* x, Node* y);
		Node* pennantSplit(Node* x);
		int bagSize();
		int elementsInBag = 0;
	public:
		int size;
		Bag(int graphSize);
		void insertBag(Node* x);
		void debug();
		void percorre(Node* node);
		void PBFS(map<int, Node*> graph);
};

//método para fazer busca em largura
void Bag::PBFS(map<int, Node*> graph){
	graph[0]->distance = 0;
	//cout << graph[0]->value << endl;
	int d = 0;
	Bag* V0 = new Bag(128); // definir constante gransize
	V0->insertBag(graph[0]);
	//V0->debug();
	cout << endl;
	map<int, Bag*> vectorBags;
	vectorBags.insert(pair<int, Bag*>(0, V0));
	while (vectorBags[d]->elementsInBag > 0) {
		cout << "Nível: " << d << endl;
		vectorBags[d]->debug();
		cout << endl;
		vectorBags.insert(pair<int, Bag*>(d + 1, new Bag(128)));
		processLayer(vectorBags[d], vectorBags[d + 1], d);
		d++;
	}
}

//retorna o tamanho da bag
int Bag::bagSize() {
	int maior = 0;
	for (int i = 0; i < 8; i++) {
		if (this->vector[i] != NULL) {
			maior = i;
		}
	}
	return maior+1;
}


//método process layer da busca
void Bag::processLayer(Bag* inBag, Bag* outBag, int d) {
	#pragma omp parallel for num_threads(4)
	for (int k = 0; k < inBag->bagSize(); k++) {
		if (inBag->vector[k] != NULL) {
			processPenant(inBag->vector[k], outBag, d);
			cout <<	"ochikubo: " << inBag->vector[k]->value << endl;
		}
			
	}
}


void Bag::processPenant(Node* inPennant, Bag* outBag, int d){
	list<Node*>pennant;
	pennant.push_back(inPennant);
	inPennant->pennantToVector(pennant, inPennant->left, 0);
	if (inPennant->pennantSize(inPennant->left,0) < 128/*verificar grainsize*/){
		cout << "tamanho do vetor: " << inPennant->pennantSize(inPennant->left, 0) << endl;
		for (int i = 0; i < pennant.size() ; i++)
		{
			inPennant = pennant.front();
			pennant.pop_front();
			#pragma omp parallel default(none)
			for(Node* x : inPennant->adjacent){
				if (x->distance == -1) {
					x->distance = d+1;
					outBag->insertBag(x);
					//cout << omp_get_thread_num() << endl;
				}
			}
		}
	}
	else {
		Node* newPennant = pennantSplit(inPennant);
		#pragma omp task default(none)
		processPenant(newPennant, outBag, d);
		processPenant(inPennant, outBag, d);
		#pragma omp taskwait default(none)
	}
	

}


//construtor da bag
Bag::Bag(int graphSize){
	if (graphSize > 0){
		size = log2(graphSize) + 1;
		vector = new Node * [size];
		for (int i = 0; i < size; i++){
			vector[i] = NULL;
		}
	}
	else {
		cerr << "Grafo deve ter pelo menos um nó" << endl;
	}	
}

//insere na bag
void Bag::insertBag(Node* x){
	int k = 0;
	while (vector[k] != NULL){
		x = pennantUnion(vector[k], x);
		vector[k] = NULL;
		k++;
		//cout << k <<endl;
	}
	vector[k] = x;
	this->elementsInBag++;
}

//split pennant
Node* Bag::pennantSplit(Node* x) {
	Node* y = new Node();
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
			cout << i << " :null" << endl;
		}
		else {
			cout << vector[i]->value << " ";
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

int main() {
	int graphSize;
	
	cin >> graphSize;
	Bag* bag = new Bag(128); // tratar exception caso 0
	map<int, Node*> graph;
	for (int i = 0; i < graphSize; i++){
		graph.insert(pair<int, Node*>(i, new Node(i)));
	}
	for (int i = 0; i < graphSize; i++){
		Node* node = graph[i];
		node->insertAdjacent(graph,i);
	}
	graph[0]->printGraph(graph);

	bag->PBFS(graph);

	//bag->debug();
}


