#include "List.h"
#include "iostream"

using namespace std;

ListNode::ListNode() {
	this->data = NODE_EMPTY;
	this->next = NULL;
	this->prev = NULL;
}

ListNode::ListNode(int data) {
	this->data = data;
	this->next = NULL;
	this->prev = NULL;
}

int ListNode::getData() const {
	return this->data;
}

ListNode* ListNode::getNext() {
	return this->next;
}

ListNode* ListNode::getPrev() {
	return this->prev;
}

void ListNode::setNext(ListNode* next) {
	this->next = next;
}

void ListNode::setPrev(ListNode* prev) {
	this->prev = prev;
}

List::List () {
	this->front = NULL;
	this->back = NULL;
	this->size = 0;
}

List::~List() {
	ListNode* toDelete;
	while(this->back != NULL) {
		toDelete = this->back;
		this->back = this->back->getPrev();
		delete toDelete;
	}
}

void List::clear() {
	this->size = 0;
	this->front = NULL;
	ListNode* toDelete;
	while(this->back != NULL) {
		toDelete = this->back;
		this->back = this->back->getPrev();
		delete toDelete;
	}
}

int List::pop_front() {
	ListNode* toDelete = this->front;
	int data;
	if(this->size == 0) return -1;
	else {
		this->size --;
		if(this->size == 0) {
			this->front = NULL;
			this->back = NULL;
		} else {
			toDelete->getNext()->setPrev(NULL);
			this->front = toDelete->getNext();
		}
		data = toDelete->getData();
		delete toDelete;
	}
	return data;
}

int List::pop_back() {
	ListNode* toDelete = this->back;
	int data;
	if(size == 0) return -1;
	else {
		this->size --;
		if(size == 0) {
			this->front = NULL;
			this->back = NULL;
		}
		else {
			toDelete->getPrev()->setNext(NULL);
			this->back = toDelete->getPrev();
		}
		data = toDelete->getData();
		delete toDelete;
	}
	return data;
}

void List::push_back(int data) {
	ListNode* newNode = new ListNode(data);
	if(size == 0) {
		this->front = newNode;
		this->back = newNode;
	}
	else {
		this->back->setNext(newNode);
		newNode->setPrev(this->back);
		this->back = newNode;		
	}
	size ++;
}
	
ListNode* List::remove(ListNode* node) {
	if(node == NULL) return NULL;
	auto toReturn = node->getNext();
	if(this->size == 1) {
		this->front = NULL;
		this->back = NULL;
	} else {
		if(node == this->front)
			this->front = this->front->getNext();
		else if(node == this->back)
			this->back = this->back->getPrev();
		if(node->getPrev() != NULL)
			node->getPrev()->setNext(toReturn);
		if(toReturn != NULL) 
			node->getNext()->setPrev(node->getPrev());
	}
	delete node;
	this->size --;	
	return toReturn;
}		
	
void List::print() {
	for(auto node = this->front; 
		node != NULL;
		node = node->getNext())
		cout << node->getData();
	cout << endl;
}	
	
ListNode* List::begin() {
	return this->front;	
}	
	
ListNode* List::end() {
	return NULL;
}
	
int List::get_front() {
	return this->front->getData();
}

int List::get_back() {
	return this->back->getData();
}

int List::get_size() {
	return this->size;
}

bool List::empty() {
	return size == 0;
}

ListNode* List::partition(ListNode* head) {
	ListNode* fast = head->getNext(), 
			*slow = head,
			*mid;
	while(fast != NULL) {
		fast = fast->getNext();
		if(fast != NULL) {
			fast = fast->getNext();
			slow = slow->getNext();
		}
	}
	mid = slow->getNext();
	slow->setNext(NULL);
	return mid;
}

void List::sort() {
	if(size <= 1) return;
	recursivelyMergeSort(front);
	// after sort, we need to restore the
	// relationship between nodes
	ListNode* temp = this->front;
	temp->setPrev(NULL);
	for(ListNode* it = temp->getNext();
		it != NULL; it = it->getNext()) {
		it->setPrev(temp); 
	}
}

ListNode* List::merge(ListNode* low, ListNode* mid) {
	ListNode* p = low,
			*q = mid;
	ListNode dummy;
	ListNode* combined = &dummy;
	while(p != NULL && q != NULL) {
		if(p->getData() < q->getData())	{		
			combined->setNext(p);
			p = p->getNext();
		}
		else {
			combined->setNext(q);
			q = q->getNext();
		}
		combined = combined->getNext();
	}		
	// concatenate the rest to buffer
	if(p == NULL) // part1 exausted 
		combined->setNext(q);
	else combined->setNext(p); // part2 exausted
	return dummy.getNext();
}

void List::recursivelyMergeSort(ListNode* &subList) {
	if(subList != NULL && subList->getNext() != NULL) {
		ListNode* mid = this->partition(subList);
		this->recursivelyMergeSort(subList);
		this->recursivelyMergeSort(mid);
		subList = this->merge(subList, mid);
	}
}
