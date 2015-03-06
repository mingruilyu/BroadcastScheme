#include "Schedule.h"
#include "iostream"

using namespace std;
Schedule::Schedule(int count) {
	this->count = count;
	this->schedule = new int[count];
	this->parent = new int[count];
	this->modifyFlag = false;
}

Schedule::~Schedule(){
	delete schedule;
	delete parent;
}

void Schedule::modify(int index,int value) {
	if (index < 0 || index >= this->count) return;
	this->previousValue = this->schedule[index];
	this->schedule[index] = value;
	this->lastIndex = index;
	this->modifyFlag = true;
}

void Schedule::rollback() {
	if (!modifyFlag) return;
	this->schedule[this->lastIndex] = this->previousValue;
	modifyFlag = false;
}

long Schedule::getPower(int index) const {
	return this->schedule[index];
}
int Schedule::getParent(int index) const {
	return this->parent[index];
}
void Schedule::setPower(int index, long power) {
	this->schedule[index] = power;
}

void Schedule::setParent(int index, int parent) {
	this->parent[index] = parent;
}

long Schedule::calculateTotalPower() const {
	long power = 0;
	for(int i = 0; i < count; i ++)
		power += this->schedule[i];
	return power;
}

void Schedule::printSchedule() const {
	for(int i = 0; i < count; i ++)
		cout << this->schedule[i] << '\t';
	cout << endl;
}
