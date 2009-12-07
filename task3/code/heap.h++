using namespace std;

template <class T>
class HeapNode {
public:
	// only internal nodes:
	T small;
	HeapNode<T> *children[3];

	// only leafs:
	int key_count;
	T keys[3];
};



template <class T>
class Heap {
public:
	HeapNode<T>* search(T key);
	void insert(T key);
	void delete_key(HeapNode<T> *node);
	HeapNode<T>* minimum();
	void decrease_key(HeapNode<T> *node, T new_key);
	HeapNode<T>* extract_min();
	Heap<T>* merge(Heap<T> *heap);
	void print(HeapNode<T> *node);
private:
	HeapNode<T>* root;
	HeadNode<T>* minimum;
};
