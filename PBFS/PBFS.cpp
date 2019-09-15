#include <iostream>
#include <math.h>
#include <map>
#include <list>
using namespace std;

class Node {
	friend class Bag;
	private:
		int value = -1;
		int distance = -1; 
		Node* left;
		Node* right;
		list<Node*> adjacent;
		int pennantSize(Node* root, int tamanho);
	public:
		Node(int value);
		Node();
		void insertAdjacent(map<int, Node*> graph);
		void printValue();
		void printGraph(map<int, Node*> graph);
};

Node::Node(int data) {
	value = data;
	left = NULL;
	right = NULL;
}

Node::Node() {
	value = -1;
	left = NULL;
	right = NULL;
}

void Node::insertAdjacent(map<int, Node*> graph){
	bool adj;
	for (int i = 0; i < graph.size(); i++) {
		cout << "O vertice " << i << " e adjacente do vertice " << this->value << "? " << endl;
		cin >> adj;
		if (adj)
			this->adjacent.push_back(graph[i]);
	}
	
}

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

void Node::printValue() {
	cout << this->value << endl;
}

int Node::pennantSize(Node* node, int tamanho) {
	pennantSize(node->left, tamanho + 1);
	pennantSize(node->right, tamanho + 1);
	return tamanho + 1;

}


class Bag {
	private:
		Node** vector;
		void processLayer(Bag* inBag, Bag* outBag, int d);
		void processPenant(Node* inPennant, Bag* outBag,int d);
		Node* pennantUnion(Node* x, Node* y);
		Node* pennantSplit(Node* x);
		int bagSize();
	public:
		int size;
		Bag(int graphSize);
		void insertBag(Node* x);
		void debug();
		void percorre(Node* node);
		void PBFS(map<int, Node*> graph);
		


};
void Bag::PBFS(map<int, Node*> graph){
	graph[0]->distance = 0;
	int d = 0;
	Bag* V0 = new Bag(128); // definir constante gransize
	V0->insertBag(graph[0]);
	map<int, Bag*> vectorBags;
	vectorBags.insert(pair<int, Bag*>(0, V0));
	while (vectorBags[d] != NULL) {
		vectorBags.insert(pair<int, Bag*>(d + 1, new Bag(128)));
		processLayer(vectorBags[d], vectorBags[d + 1], d);
		d++;
	}
}

int Bag::bagSize() {
	int maior = 0;
	for (int i = 0; i < 8; i++) {
		if (this->vector[i] != NULL) {
			maior = i;
		}
	}
	return maior+1;
}

void Bag::processLayer(Bag* inBag, Bag* outBag, int d) {
	for (int k = 0; k < inBag->bagSize(); k++) {
		if (inBag->vector[k] != NULL)
			processPenant(inBag->vector[k], outBag, d);
	}
}


void Bag::processPenant(Node* inPennant, Bag* outBag, int d){
	if (inPennant->pennantSize(inPennant->left, 0)) {
		for (size_t i = 0; i < length; i++)
		{
			for (Node* x : inPennant->adjacent) {

			}
		}
		
	}
	
}

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

void Bag::insertBag(Node* x){
	int k = 0;
	while (vector[k] != NULL){
		x = pennantUnion(vector[k], x);
		vector[k] = NULL;
		k++;
		//cout << k <<endl;
	}
	vector[k] = x;
}

Node* Bag::pennantSplit(Node* x) {
	Node* y = new Node();
	x->left = y->right;
	y->right = NULL;
	return y;
}

Node* Bag::pennantUnion(Node* x, Node* y) {
	y->right = x->left;
	x->left = y;
	return x;
}

void Bag::debug() {
	for (int i = 0; i < size; i++)
	{
		if (this->vector[i] == NULL){
			cout << i << " :null" << endl;
		}
		else {
			//cout << vector[i]->value << endl;
			//percorre(vector[i]->left);
		}
	}
}

void Bag::percorre(Node* node) {
	if(node != NULL){
		cout << node->value << endl;
		percorre(node->left);
		percorre(node->right);
	}
}

int main() {
	int graphSize;
	cin >> graphSize;
	Bag* bag = new Bag(graphSize); // tratar exception caso 0
	map<int, Node*> graph;
	for (int i = 0; i < graphSize; i++){
		graph.insert(pair<int, Node*>(i, new Node(i)));
	}
	for (int i = 0; i < graphSize; i++){
		Node* node = graph[i];
		node->insertAdjacent(graph);
	}
	graph[0]->printGraph(graph);
	
	
	//bag->debug();
}


