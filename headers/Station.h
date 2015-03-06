#ifndef STATION_H
#define STATION_H

#include "memory"
#include "List.h"

#define STATION_UNMARK	false
#define STATION_MARK 	true
#define STATION_SOURCE	-1
#define STATION_IN_SET	-1
#define STATION_UNDETERMINED	-1
using namespace std;

class Station;
class StationSet;

class Station{
private:
	int x;
	int y;
	int index;
	bool mark;
	static int stationCount;
	// relayTarget keep all the other stations
	// that are using this station as relay
	List relayTarget;
	// reachabel records all the stations in
	// the set that lie within the radius of
	// the current transmission power. Each
	// time the trans mission power is changed,
	// this list should be updated.
	List reachable;
	// the current transmission power.
	long relayPower;
public:
	Station();
	void setCoordinate(const int x, const int y);
	void setMark(const bool mark);
	// get the power required to transmit from the
	// station to the target
	long getPower(const Station& target) const;
	long getXDiff(const Station& target) const;
	int getIndex() const;
	static bool powerComparator(const Station& a, const Station& b);
	bool isMarked() const;
	// add the target to relayTarget list
	// so this station will become a relay
	// station of the target station
	void addRelayTarget(int target);
	// get the relay power of this station
	long getRelayPower() const;
	void setRelayPower(long power);
	// add the index to reachable. this means
	// that the index station is currently
	// sitting within the relay power radius
	void addReachable(const int index);
	// pop the last element of the reachable list
	int popReachable();
	// remove the node pointed by toRemove
	// from the reachable list
	ListNode* removeReachable(ListNode* toRemove);
	void clearReachable();
	ListNode* getReachable();
	List* getRelayTarget();
	void sortRelayTarget();
	void commit(int index);
	int getX() const;
	int getY() const;
};

class StationSet{
private:
	int capacity;
	int size;
	bool modifyFlag;
	Station* stations;
	int* set;
public:
	StationSet(const int count);
	~StationSet();
	void unmark();
	void setStations(Station* stations);
	int getSize() const;
	void insert(int index);
	/* get the pos th element in the set. Note that the
	 * pos th element is not the index of that station*/
	int get(const int pos);
	void printSet() const;
	// roll back the last change to the set
	void rollback();
	// return true if the set is full
	bool full() const;
};

#endif
