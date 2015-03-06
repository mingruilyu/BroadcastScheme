#include "Station.h"
#include "iostream"
#include "math.h"
using namespace std;

int Station::stationCount = 0;

Station::Station(){
	this->index = stationCount ++;
	this->mark = STATION_UNMARK;
	this->relayPower = 0;
}
	
void Station::setCoordinate(int x, int y) {
	this->x = x;
	this->y = y;
}

long Station::getPower(const Station& station) const {
	return (this->x - station.x) * (this->x - station.x) +
			(this->y - station.y) * (this->y - station.y);
}

long Station::getXDiff(const Station& station) const {
	return (this->x - station.x) * (this->x - station.x);
}

void Station::setMark(bool mark) {
	this->mark = mark;
}

bool Station::isMarked() const {
	return this->mark;
}

int Station::getIndex() const {
	return this->index;
}

void Station::addRelayTarget(int target) {
	this->relayTarget.push_back(target);
}

long Station::getRelayPower() const {
	return this->relayPower;
}

void Station::setRelayPower(long power) {
	this->relayPower = power;
}

void Station::addReachable(int index) {
	this->reachable.push_back(index);
}

int Station::popReachable() {
	return this->reachable.pop_back();
}

void Station::commit(int index) {
	this->reachable.push_back(index);
}

ListNode* Station::removeReachable(ListNode* toRemove) {
	return this->reachable.remove(toRemove);
}

void Station::clearReachable() {
	this->reachable.clear();
}

ListNode* Station::getReachable() {
	return this->reachable.begin();
}

List* Station::getRelayTarget() {
	return &this->relayTarget;
}

void Station::sortRelayTarget() {
	this->relayTarget.sort();
}

int Station::getX() const {
	return this->x;
}

int Station::getY() const {
	return this->y;
}

StationSet::StationSet(int count) {
	this->capacity = count;
	this->set = new int[count];
}
void StationSet::setStations(Station* stations) {
	this->stations = stations;
}

StationSet::~StationSet() {
	delete[] set;
}

void StationSet::unmark(){
	for(int i = 0; i < size; i ++) 
		stations[this->set[i]].setMark(STATION_UNMARK);
}

int StationSet::getSize() const {
	return size;
}

void StationSet::insert(int index) {
	// add station index to set
	this->set[size ++] = index;
	this->modifyFlag = true;
}

bool StationSet::full() const {
	return size == capacity;
}

int StationSet::get(const int pos) {
	return set[pos];
}

void StationSet::printSet() const {
	for(int i = 0; i < size; i ++)
		cout << this->set[i] << '\t';
	cout << endl;
}

void StationSet::rollback() {
	if(modifyFlag) this->size --;
	else cerr << "Trying to rollback an uncommited operation";
	this->modifyFlag = false;
}

