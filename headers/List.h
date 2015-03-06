#ifndef LIST_H
#define LIST_H
#define NODE_EMPTY -1

class ListNode {
private:
	ListNode* next;
	int 	data;
	ListNode* prev;
public:
	ListNode();
	ListNode(int data);
	int getData() const;
	ListNode* getNext();
	ListNode* getPrev();
	void setNext(ListNode* next);
	void setPrev(ListNode* prev);
};

class List {
private:
	// front points to the 1st element 
	// of the list 
	ListNode* front;
	// back points to the last element
	ListNode* back;
	int size;
public:
	List ();
	~List();
	void clear();
	int pop_front();
	int pop_back();
	void push_back(int data);
	ListNode* remove(ListNode* node);
	void print();
	ListNode* begin();
	ListNode* end();
	int get_front();
	int get_back();
	int get_size();
	bool empty();
	void sort();
private:
	ListNode* partition(ListNode* head);
	ListNode* merge(ListNode* low, ListNode* mid);
	void recursivelyMergeSort(ListNode* &subList);
};

#endif
