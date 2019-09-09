#include <iostream>

using namespace std;

class Node{
	private:
		int value = 0;
		Node* left;
		Node* right;
	public:
		Node(int value);
		void printValue();
};

Node::Node(int data){
		value = data;
		left = nullptr;
		right = nullptr;
}

void Node:: printValue(){
	cout << this->value;
	}

int main(){
	Node* newNode = new Node(5);
	newNode->printValue();
}
