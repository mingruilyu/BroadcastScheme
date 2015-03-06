#ifndef BROADCAST_SCHEME_H
#define BROADCAST_SCHEME_H

#include "Station.h"
#include "Schedule.h"
#define METHOD_ACSF		1
#define METHOD_ACSET	2
#define METHOD_ACF		3
#define METHOD_DC		4
#define METHOD_OWN		5
#define METHOD_DP_RECOMP	6
#define METHOD_DP_NO_RECOMP 7
#define METHOD_DP_ITER	8
#define GREATEST_POWER	100000000000
#define GREATEST_INDEX  10000
#define COMMIT			true
#define ATTEMP			false
#define SORT_BY_DIST2SOURCE true
#define SORT_BY_X_AXIS		false
#define SOURCE_INPUT	true
#define SOURCE_LEFTMOST false
using namespace std;

class BroadcastScheme {
private:
	int source;
	int method;
	int stationCount;
	int leftMost;
	/* priorityQueue is used to maintain a queue of
	 * stations sorted by the distance square(power)
	 * to the source station 						*/
	Station* stations;
	Schedule schedule;
	StationSet set;
private:
	/* initialize coordinate of stations and priority queue */
	void initializeStations();
	/* generate transmission tree from schedule and station
	 * set. This procedure does not change schedule or
	 * station.											*/
	long generateTree(const bool add);
	void unmarkStations();
	// find the station that takes the lowest
	// power to transmit to the target station
	int  findRelayStation(int target);
	// use the relay station to transmit to target
	void addRelayStation(int relay, int target);
	void AddClosest2SourceFirst();
	void AddClosest2SetFirst();
	void AddCheapestFirst();
	void DivideConquer();
	void ownMethod();
	void dynamicProgRecomp();
	void dynamicProgNoRecomp();
	void dynamicProgIter();
	long recursivelyDPRecomp(int startPos, int endPos,
							  int* sequence,
							  int maxTransmission, int* pivot, int* firstHalf);
	long recursivelyDPNoRecomp(int startPos, int endPos,
								int* sequence,
								int m, long*** val,
								int*** keyP, int*** keyQ);
	// return a list of station indice that lie within
	// relay's power radius. The indice are sorted from
	// small to big
	inline	List* getEncompassed(const int relay, const long power, const bool add);
	void recursivelyMergeSort(int low, int high, int* array, int* buffer, bool standard) const;
	void merge(int low, int high, int* array, int* buffer, bool standard) const;
	void recursivelyPrintNodes(int root) const;
	// sort the array using mergesort
	void sort(int* array, bool standard) const;
	/* update the schedule associated*/
	void updateSchedule();
	void constructScheduleRecomp(int startPos, int endPos, 
						   		int* sequence, int maxTransmission);
	void constructScheduleNoRecomp(int startPos, int endPos, int* sequence, int m, long*** val, int*** keyP, int*** keyQ);
	void reorderGroup(int* group, int lstGrpStart, int lstGrpEnd, 
					  int newGrpStart, int newGrpEnd);
	int findRelayFromGroup(int target, int grpStart, int grpEnd, int* group);	
public:
	BroadcastScheme(const int method, const int stationCount, const int source);
	~BroadcastScheme();
	/* generate schedule using the specified method, result
	 * is stored in schedule					*/
	void generateSchedule();
	void printScheme(bool mode) const;
	void printTree() const;
};
#endif
