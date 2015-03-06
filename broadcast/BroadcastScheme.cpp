#include "BroadcastScheme.h"
#include "iostream"
#include "List.h"
#include "math.h"
#include "algorithm"
#include "memory"
using namespace std;

BroadcastScheme::BroadcastScheme(int method, int stationCount, int source)
				: schedule(stationCount), set(stationCount) {
	this->method = method;
	this->stationCount = stationCount;
	this->source = source;

	stations = new Station[stationCount];
	set.setStations(stations);
	initializeStations();
}

BroadcastScheme::~BroadcastScheme() {
	delete[] stations;
}

void BroadcastScheme::merge(int low, int high, int* array, int* buffer, bool standard) const{
	// mid point is the pivot that divide the
	// space into two parts: [low, mid], [low + 1, high]
	int mid = (low + high) >> 1,
		size = high - low + 1;
	int p = low, q = mid + 1, i = 0;
	
	while(p <= mid && q <= high) {
		if(standard == SORT_BY_DIST2SOURCE) {
			long power1, power2;
			power1 = stations[array[p]].getPower(stations[source]);
			power2 = stations[array[q]].getPower(stations[source]); 
			if(power1 < power2 
				|| (power1 == power2 && array[p] < array[q]))
				buffer[i ++] = array[p ++];
			else buffer[i ++] = array[q ++];
		} else {
			int x1, x2;
			x1 = stations[array[p]].getX();
			x2 = stations[array[q]].getX();
			if(x1 < x2 || (x1 == x2 && array[p] < array[q]))
				buffer[i ++] = array[p ++];
			else buffer[i ++] = array[q ++];
		}
	}
	// concatenate the rest to buffer
	if(p > mid) {// part1 exausted 
		for(; i < size; i ++)
			buffer[i] = array[q ++];
	} else { // part2 exausted
		for(; i < size; i ++)
			buffer[i] = array[p ++];
	}
	// copy buffer to original array
	for(int i = low; i <= high; i ++)
		array[i] = buffer[i - low];
}

void BroadcastScheme::recursivelyMergeSort(int low, int high, int* array, int* buffer, bool standard) const {
	if(high == low) return;
	int mid = (low + high) >> 1;
	recursivelyMergeSort(low, mid, array, buffer, standard);
	recursivelyMergeSort(mid + 1, high, array, buffer, standard);
	this->merge(low, high, array, buffer, standard);
}

void BroadcastScheme::sort(int* array, bool standard) const{
	int* buffer = new int[stationCount];
	this->recursivelyMergeSort(0, stationCount - 1, array, buffer, standard);
}

void BroadcastScheme::initializeStations() {
	// initialize the coordinate of all stations
	int x, y;
	for(int i = 0; i < stationCount; i ++) {
		cin >> x >> y;
		stations[i].setCoordinate(x, y);
	}
}

List* BroadcastScheme::getEncompassed(int relayIndex, long power, bool add) {
	Station* target, 
			*relay = &stations[relayIndex];
	// if the power is 0, no stations will
	// be transmit to from this station
	if(power == 0) {
		relay->clearReachable();
		relay->setRelayPower(0);
		return relay->getRelayTarget();
	}
	Station* justAdded = &stations[set.get(set.getSize() - 1)];
	int toAdd = justAdded->getIndex();
	long relayPower, crtRelayPower = 0;
	// update the reachable list, do only
	// when the relay power of this station changes.
	// We use a reachable list to speed up the
	// process of obtaining all stations that 
	// lie in the power radius of the station
	// indexed by relayIndex. We notice that 
	// most of time, the relay power does not
	// change. So we dont really need to traverse
	// the whole set to calculate the power of 
	// stations. According to experiment, this
	// optimization contribute to 50% perfor-
	// mance gain
	if(power > relay->getRelayPower()) {
		// the new power is greater than the
		// relay's current relay power, need
		// to update the whole reachable list
		relay->clearReachable();
		for (int i = 1; i < this->set.getSize(); i ++) {
			// since souce must have been marked, the next if
			// branch will never be reached for the source.
			// start from the 2nd in the set
			// toAdd should be explicitly excluded
			if(this->set.get(i) == relayIndex
				|| this->set.get(i) == toAdd) continue;
			target = &stations[this->set.get(i)];
			relayPower = target->getPower(*relay);
			if (relayPower <= power) 
				relay->addReachable(target->getIndex());
		}
	} else {
		// eliminate the stations that now lie outside
		// the power radius
		if(power < relay->getRelayPower()) {
			for(ListNode* it = relay->getReachable(); 
		    	it != NULL;) {
				target = &(stations[it->getData()]);
				relayPower = target->getPower(*relay);
				if(relayPower > power)
					it = relay->removeReachable(it);
				else it = it->getNext();
			}
		}
	}
	// toAdd is the station just added to the set.
	// Since the set may be later rolled back, we 
	// treat the toAdd particularly.
	if(toAdd != relayIndex
		&& relay->getPower(*justAdded) <= power) {
		// we add toAdd to reachable
		if(add) relay->commit(toAdd);
		else {
			target = &stations[toAdd];
			relayPower = target->getPower(*relay);
			if(!target->isMarked()) {
				// we will add this to relayTarget
				relay->addRelayTarget(justAdded->getIndex());
				crtRelayPower = relayPower;
			}
		}
	} 
	// select from reachable the staion that
	// is not marked and lie in power radius
	// remove those stations that lie outside
	// power radius and already marked
	// if the 
	for(auto it = relay->getReachable(); 
	    it != NULL; it = it->getNext()) {
		target = &stations[it->getData()];
		relayPower = target->getPower(*relay);
		if(!target->isMarked()) {
			relay->addRelayTarget(it->getData());
			if(relayPower > crtRelayPower)
				crtRelayPower = relayPower;
		}
	}
	relay->setRelayPower(crtRelayPower);
	relay->sortRelayTarget();
	return relay->getRelayTarget();
}

void BroadcastScheme::printTree() const {
	this->recursivelyPrintNodes(source);
	cout << endl;
}

void BroadcastScheme::recursivelyPrintNodes(int root) const {
	Station* rootStation = &stations[root];
	List* targetList = rootStation->getRelayTarget();
	cout << rootStation->getIndex() << '(';
	for(ListNode* target = targetList->begin(); 
		target != NULL;
		target = target->getNext())
		this->recursivelyPrintNodes(target->getData());
	cout << ')';
}

long BroadcastScheme::generateTree(bool add) {
	List queue;
	List* unmarkedList;
	int target; // index of the sma1llest indexed station
	this->set.unmark(); // unmark all stations in station set
	long totalPower = 0;
	queue.push_back(source);
	stations[source].setMark(STATION_MARK);
	
	for(int relay; !queue.empty(); ) {
		relay = queue.pop_front();
		// get a list of all unmarked stations in set
		// that in the radius power of the relay
		unmarkedList = getEncompassed(relay,
						  			  schedule.getPower(relay),
									  add);
		totalPower += stations[relay].getRelayPower();
		while(!unmarkedList->empty()) {
			// get index of the smallest indexed station
			target = unmarkedList->pop_front();
			stations[target].setMark(STATION_MARK);
			// if add bit is set, the change to schedule 
			// should be commited.
			if(add) schedule.setParent(target, relay);
			queue.push_back(target);
		}
	}
	return totalPower;
}

void BroadcastScheme::addRelayStation(int relay, int target) {
	int relayPower = max((long) schedule.getPower(relay),
					stations[relay].getPower(stations[target]));
	schedule.modify(relay, relayPower);
	// need not to rollback, commit the change
	// to the schedule in the tree precedure.
	this->generateTree(COMMIT);
	this->updateSchedule();
}

void BroadcastScheme::updateSchedule() {
	int index;
	for(int i = 0; i < set.getSize(); i ++) {
		index = set.get(i);
		schedule.setPower(stations[index].getIndex(),
						  stations[index].getRelayPower());
	}
}

int BroadcastScheme::findRelayFromGroup(int target, 
	int grpStart, int grpEnd, int* group) {
	int smallestRelayIndex, candidate;
	long totalPower, relayPower, lowestPower = GREATEST_POWER;
	// use source
	relayPower = max((long) schedule.getPower(candidate),
					 stations[candidate].
					 getPower(stations[target]));
	for(int j = grpStart; j <= grpEnd; j ++) {
		candidate = group[j];
		// calculate the new power needed 
		// had candidate become the relay
		relayPower = max((long) schedule.getPower(candidate),
						 stations[candidate].
						 getPower(stations[target]));
		schedule.modify(candidate, relayPower);
		// generate Transmission tree and get
		// the total power
		totalPower = this->generateTree(ATTEMP);
		if(totalPower < lowestPower 
			|| (totalPower == lowestPower 
				&& candidate < smallestRelayIndex)) {
			// record the lowerPower and corresponding relay
			lowestPower = totalPower;
			smallestRelayIndex = candidate;
		}
	}						
	return smallestRelayIndex;
}

int BroadcastScheme::findRelayStation(int target) {
	int smallestRelayIndex, candidate;
	long totalPower, relayPower, lowestPower = GREATEST_POWER;
	for(int j = 0; j < set.getSize(); j ++) {
		candidate = set.get(j);
		if(candidate == target) continue; 
		// calculate the new power needed had candidate become the relay
		relayPower = max((long) schedule.getPower(candidate),
						 stations[candidate].getPower(stations[target]));
		schedule.modify(candidate, relayPower);
		// generate Transmission tree and get
		// the total power
		totalPower = this->generateTree(ATTEMP);
		if(totalPower < lowestPower 
			|| (totalPower == lowestPower 
				&& candidate < smallestRelayIndex)) {
			// record the lowerPower and corresponding relay
			lowestPower = totalPower;
			smallestRelayIndex = candidate;
		}
		schedule.rollback(); // rollback the modification
				     // so the schedule can be reused
	}						
	return smallestRelayIndex;
}

void BroadcastScheme::generateSchedule() {
	// all stations are unmarked upon initialization
	// no need to unmark again.
	switch(method) {
	case METHOD_ACSF:
		AddClosest2SourceFirst();
		break;
	case METHOD_ACSET:
		AddClosest2SetFirst();
		break;
	case METHOD_ACF:
		AddCheapestFirst();
		break;
	case METHOD_DC:
		DivideConquer();
		break;
	case METHOD_OWN:
		ownMethod();
		break;
	case METHOD_DP_RECOMP:
		dynamicProgRecomp();
		break;
	case METHOD_DP_NO_RECOMP:
		dynamicProgNoRecomp();
		break;
	case METHOD_DP_ITER:
		dynamicProgIter();
		break;
	}
}

void BroadcastScheme::AddClosest2SourceFirst() {
	int target; // the station to be added to set in each iteration
	int relay; // the station selected to transmit to crtStation
	int* priorityQueue = new int[stationCount];
	int queueHead = 0;
	for(int i = 0; i < stationCount; i ++)
		priorityQueue[i] = i;
	// sort the priority queue by distance to source
	this->sort(priorityQueue, SORT_BY_DIST2SOURCE);
	// eliminate the source from the queue
	queueHead ++;
	// add source station to set
	set.insert(source);
	schedule.setParent(source, STATION_SOURCE);
	for(int i = 1; i < stationCount; i ++) {
		// iterate stationCount - 1. 1st element 
		// of set is always source. start from the 2nd
		// get the next station that is closet to source
		target = priorityQueue[queueHead ++];
		//priorityQueue.pop_front();
		// add new target station to set
		set.insert(target);
		// find relay station to transmit to the target station
		relay = this->findRelayStation(target);
		// use relay to transmit to target station
		this->addRelayStation(relay, target);
	}
} 

void BroadcastScheme::AddClosest2SetFirst() {
	int target; // the station to be added to set in each iteration
	long newPower, lowestPower = GREATEST_POWER;
	int relay; // the station selected to transmit to target station
	int smallestTargetIndex = GREATEST_INDEX;
	/* REASON why not to use heap:
	 * Everytime a new station is added to the set, the power of 
	 * all stations to the closest station in the set should be
	 * recalculated, taking into consideration the newly added 
	 * station. It is highly possible that many has changed.
	 * This results in O(nlogn) to restore the heap. If we only use
	 * array, we just do a linear scan to find out the closest
	 * station to the set, while updating the power. This takes
	 * O(n)														 */
	// power is used to store the shortest distance from the set
	// to the station, once added, the distance become GREATEST_POWER
	// This means that the station is no long considered. Everytime
	// a new station is added to the set, do a scan and update each
	// station's power
	long* power = new long[stationCount];
	for(int newTarget = 0; newTarget < stationCount; newTarget ++) {
		// source is explicitly excluded
		if(newTarget == this->source) continue; 
		power[newTarget] = stations[newTarget].getPower(stations[this->source]);
		if(power[newTarget] < lowestPower
		   || (power[newTarget] == lowestPower
				&& newTarget < smallestTargetIndex)) {
			lowestPower = power[newTarget];
			smallestTargetIndex = newTarget;
		}
	}
	target = smallestTargetIndex;
	// add source station to set
	set.insert(source);
	schedule.setParent(source, STATION_SOURCE);
	power[source] = STATION_IN_SET;
	for(int i = 1; i < stationCount; i ++) {
		// iterate stationCount - 1 times,
		// source is implicitly excluded because it
		// will never be chosen as target
		set.insert(target);
		power[target] = STATION_IN_SET;
		//find relay to transmit to target
		relay = this->findRelayStation(target);
		// use relay to transmit to target
		this->addRelayStation(relay, target);
		
		// new station has been added to the set,
		// update power array and closet array
		lowestPower = GREATEST_POWER;
		smallestTargetIndex = GREATEST_INDEX;
		for(int newTarget = 0; newTarget < stationCount; newTarget ++) {
			if(power[newTarget] == STATION_IN_SET) continue;
			// calculate the power required to transmit from
			// newTarget station to the newly added station, if
			// it is less than the current power, power should
			// be updated, so is the closest.
			newPower = stations[newTarget].getPower(stations[target]);
			if(newPower < power[newTarget])
				power[newTarget] = newPower;
			if(power[newTarget] < lowestPower
				|| (power[newTarget] == lowestPower
					&& newTarget < smallestTargetIndex)) {
				lowestPower = power[newTarget];
				smallestTargetIndex = newTarget;
			} 
		}
		target = smallestTargetIndex;
	}
}

void BroadcastScheme::reorderGroup(int* group, int lstGrpStart,
	int lstGrpEnd, int newGrpStart, int newGrpEnd) {
	long lowestPower, crtPower, lowestPowerIndex;
	long* power = new long[newGrpEnd - newGrpStart + 1];
	// calculate the lowest power needed to transmit
	// from the last group to new stations in new group.
	for(int i = newGrpStart; i <= newGrpEnd; i ++) {
		lowestPower = GREATEST_POWER;
		for(int j = lstGrpStart; j <= lstGrpEnd; j ++) {
			crtPower = stations[group[j]].
						getPower(stations[group[i]]);
			if(crtPower < lowestPower)
				lowestPower = crtPower;
		}
		power[i - newGrpStart] = lowestPower;
	}
	// sort using selective sort
	for(int i = newGrpStart; i < newGrpEnd; i ++) {
		lowestPower = power[i];
		lowestPowerIndex = i;
		for(int j = i + 1; j <= newGrpEnd; j ++) {
			if(lowestPower > power[j]) {
				lowestPower = power[j];
				lowestPowerIndex = j;
			}
		}
		// swap the lowestPowerIndex and current i
		long temp = group[i];
		group[i] = group[lowestPowerIndex];
		group[lowestPowerIndex] = temp;
	}
	delete[] power;
}

void BroadcastScheme::ownMethod() {
	// the station to be added to set in each iteration
	int target; 
	// the station selected to transmit to target station
	int relay;
	int groupCount = 6,
		standardSize = (int)floor(
						((float)stationCount / groupCount)),
		groupSize = standardSize;
	int* priorityQueue = new int[stationCount];	
	for(int i = 0; i < stationCount; i ++)
		priorityQueue[i] = i;
	// sort the priority queue by distance to source
	this->sort(priorityQueue, SORT_BY_DIST2SOURCE);
	// add source station to set
	set.insert(this->source);
	
	for(int group = 0; group < groupCount; group ++) {
		// do a linear scan over inSet to find the unmarked station
		// each iteration add one station in a group
		// while loop iterates until all stations in a
		// group have been added to the set
		if(group == 0) {
			// for the stations in first group, 
			// add closest to source first
			for(int i = 1; i < groupSize; i ++) {
				target = priorityQueue[i];
				set.insert(target);
				// find relay station to transmit to the target station
				relay = this->findRelayStation(target);
				// use relay to transmit to target station
				this->addRelayStation(relay, target);
			}
			continue;
		}

		this->reorderGroup(priorityQueue, 
						   (group - 1) * standardSize,
						   group * standardSize - 1,
							group * standardSize, 
							group * standardSize + groupSize - 1);
		
		for(int i = 0; i < groupSize; i ++) {
			target = priorityQueue[group * standardSize + i];
			set.insert(target);
			// find the relay to transmit to target station
			// findRelay function makes sure the return value is 
			// the best relay station with the smallest index
			relay = this->findRelayStation(target);
			// use relay to transmit to target station
			this->addRelayStation(relay, target);
		}
		if(group == groupCount - 2)
			groupSize = stationCount - groupSize * (groupCount - 1);
	}

	schedule.setParent(source, STATION_SOURCE);
	delete[] priorityQueue;
}

void BroadcastScheme::AddCheapestFirst() {
	 // the target station with the smallest index
	int smallestTargetIndex = GREATEST_INDEX;
	long totalPower, lowestPower, relayPower;
	int candidate; // temp station to relay
	// inSet is used to recorde all stations that has not
	// been added to the set. Those not in the set are marked
	// false. notInSetCount tells how many stations that are 
	// not in the set;
	bool* inSet = new bool[stationCount];
	int* relay = new int[stationCount];
	set.insert(source);
	schedule.setParent(source, STATION_SOURCE);
	inSet[source] = true;

	while(!set.full()) {
		lowestPower = GREATEST_POWER;
		smallestTargetIndex = GREATEST_INDEX;
		// do a linear scan over inSet to find the unmarked station
		for(int target = 0; target < stationCount; target ++) {
			if(target == source || inSet[target]) continue;
			set.insert(target);
			// find the relay to transmit to target station
			// findRelay function makes sure the return value is 
			// the best relay station with the smallest index
			candidate = this->findRelayStation(target);
			relay[target] = candidate;
			// use relay to transmit to target station
			relayPower = max((long) schedule.getPower(candidate),
							stations[candidate].getPower(stations[target]));
			schedule.modify(candidate, relayPower); // need not to rollback
			totalPower = this->generateTree(ATTEMP);
			schedule.rollback();
			if(totalPower < lowestPower
				|| (totalPower == lowestPower
					&& target < smallestTargetIndex)) {
				lowestPower = totalPower;
				smallestTargetIndex = target;
			}
			set.rollback();
		}
		// smallestTargetIndex is currently the target station that
		// use the least power to transmit to.
		// Its best relay station is stored in relay array
		// add the smallestTargetIndex station to inSet
		set.insert(smallestTargetIndex);
		inSet[smallestTargetIndex] = true;
		this->addRelayStation(relay[smallestTargetIndex], smallestTargetIndex);
	}
	delete[] inSet;
	delete[] relay;
}

// Consider what if there is only one group

void BroadcastScheme::DivideConquer() {
	int groupCount = 4,
		groupSize,
		standardSize = (int)floor(((float)stationCount / groupCount));
	int addCount;
	int target, candidate, smallestTargetIndex;
	long totalPower, lowestPower, relayPower;
	bool* inSet = new bool[stationCount];
	int* relay = new int[stationCount];
	int* priorityQueue = new int[stationCount];
	for(int i = 0; i < stationCount; i ++)
		priorityQueue[i] = i;
	// sort the priority queue by distance to source
	this->sort(priorityQueue, SORT_BY_DIST2SOURCE);
	// add source station to set
	set.insert(source);
	inSet[source] = true;
	schedule.setParent(source, STATION_SOURCE);
	
	groupSize = standardSize;
	addCount = groupSize - 1;
	for(int group = 0; group < groupCount; group ++) {
		// do a linear scan over inSet to find the unmarked station
		// each iteration add one station in a group
		// while loop iterates until all stations in a
		// group have been added to the set
		while(addCount > 0 && !set.full()) {	
			lowestPower = GREATEST_POWER;
			smallestTargetIndex = GREATEST_INDEX;
			for(int i = 0; i < groupSize; i ++) {
				target = priorityQueue[group * standardSize + i];
				if(inSet[target]) continue;
				set.insert(target);
				// find the relay to transmit to target station
				// findRelay function makes sure the return value is 
				// the best relay station with the smallest index
				candidate = this->findRelayStation(target);
				relay[target] = candidate;
				// use relay to transmit to target station
				relayPower = max((long) schedule.getPower(candidate),
				   				 stations[candidate].getPower(stations[target]));
				schedule.modify(candidate, relayPower); // need not to rollback
				totalPower = this->generateTree(ATTEMP);
				schedule.rollback();
				if(totalPower < lowestPower
					|| (totalPower == lowestPower
					&& target < smallestTargetIndex)) {
					lowestPower = totalPower;
					smallestTargetIndex = target;
				}
				set.rollback();
			}
			// smallestTargetIndex is currently the target 
			// station that use the least power to transmit to.
			// Its best relay station is stored in relay array
			// add the smallestTargetIndex station to inSet
			set.insert(smallestTargetIndex);
			inSet[smallestTargetIndex] = true;
			this->addRelayStation(relay[smallestTargetIndex], smallestTargetIndex);
			addCount --;
		}
		// get ready for the last group, so it 
		// has to be groupCount - 2
		if(group == groupCount - 2) {
			addCount = stationCount - groupSize * (groupCount - 1);
			groupSize = addCount;
		}
		else addCount = groupSize;
	}
	delete[] inSet;
	delete[] relay;
}

void BroadcastScheme::dynamicProgRecomp() {
	int* transmitSequence = new int[stationCount];
	for(int i = 0; i < stationCount; i ++)
		transmitSequence[i] = i;
	this->sort(transmitSequence, SORT_BY_X_AXIS);
	this->leftMost = transmitSequence[0];
	this->constructScheduleRecomp(0, stationCount - 1,
								   transmitSequence,
								   source + 1);
	schedule.setParent(transmitSequence[0], STATION_SOURCE);
	delete[] transmitSequence;
}

long BroadcastScheme::recursivelyDPRecomp(int startPos, int endPos,
	int* sequence, int maxTransmission, int* pivot, int* firstHalf) {
	int start = sequence[startPos];
	int end = sequence[endPos];
	if(endPos - startPos <= 1 || maxTransmission == 1) {
		*pivot = startPos;
		*firstHalf = 1;
		return stations[start].getXDiff(stations[end]);	
	}

	maxTransmission =  min(endPos - startPos, maxTransmission);

	int partition = startPos,
		partitionCount = 1;
	long power1, power2, totalPower,
		// direct transmission is considered in particular
		 minPower = stations[start].getXDiff(stations[end]);
	
	for(int i = 1; i < maxTransmission; i ++) {
		// ensure 1st and 2nd half has at least
		// one transmision therefore no infinite
		// recursion
		for(int j = startPos + 1; j < endPos; j ++) {
			// transmission is always from left to right
			power1 = this->recursivelyDPRecomp(startPos, j, 
												sequence, i, pivot, firstHalf);
			power2 = this->recursivelyDPRecomp(j, endPos,
												sequence,
												maxTransmission - i, pivot, firstHalf);
			totalPower = power1 + power2;
			if(totalPower < minPower) {
				minPower = totalPower;
				partition = j;
				partitionCount = i;
			}
		}
	}
	*pivot = partition;
	*firstHalf = partitionCount;
	return minPower;
}

void BroadcastScheme::constructScheduleRecomp(int startPos,
		int endPos, int* sequence, int maxTransmission) {
	int start = sequence[startPos];
	int end = sequence[endPos];
	int pivot, firstHalf;
	this->recursivelyDPRecomp(startPos, endPos,
			 				  sequence,
							  maxTransmission,
							  &pivot, &firstHalf);

	if(pivot != startPos) {
		// if not direct transmission from
		// startPos to endPos
		this->constructScheduleRecomp(startPos, pivot,
									  sequence,
									  firstHalf);
		this->constructScheduleRecomp(pivot, endPos,
									  sequence,
									  maxTransmission - firstHalf);
	} else {
		schedule.modify(start, 
						stations[start].getXDiff(stations[end]));
		schedule.setParent(end, start);
	}
}

void BroadcastScheme::constructScheduleNoRecomp(
int startPos, int endPos, int* sequence, int m, 
long*** val, int*** keyP, int*** keyQ) {
	int start = sequence[startPos];
	int end = sequence[endPos];
	m = min(endPos - startPos, m);
	if(keyP[startPos][endPos][m] != startPos) {
		// if not direct transmission from
		// startPos to endPos
		this->constructScheduleNoRecomp(startPos, 
										keyP[startPos][endPos][m],
										sequence,
										keyQ[startPos][endPos][m],
							   			val, keyP, keyQ);
		this->constructScheduleNoRecomp(keyP[startPos][endPos][m],
										endPos,
										sequence,
										m - keyQ[startPos][endPos][m],
										val, keyP, keyQ);
	} else {
		schedule.modify(start, val[startPos][endPos][m]);
		schedule.setParent(end, start);
	}
}

void BroadcastScheme::dynamicProgNoRecomp() {
	// keyP takes down from start to end
	// which station should be relay when
	// the keyQ relays can be used. keyQ
	// is the least number of relay needed
	// to get a minimum val.
	int*** keyP = new int**[stationCount];
	int*** keyQ = new int**[stationCount];
	long*** val = new long**[stationCount];
	for(int i = 0; i < stationCount; i ++) {
		keyP[i] = new int*[stationCount];
		keyQ[i] = new int*[stationCount];
		val[i] = new long*[stationCount]; 
		for(int j = 0; j < stationCount; j ++) {
			keyP[i][j] = new int[this->source + 2];
			keyQ[i][j] = new int[this->source + 2];
			val[i][j] = new long[this->source + 2];
		}
	}
	// initialize val array with distance between stations
	for(int i = 0; i < stationCount; i ++) {
		for(int j = 0; j < stationCount; j ++) {
			for(int k = 1; k < this->source + 2; k ++) {
				keyP[i][j][k] = i;
				keyQ[i][j][k] = 1;
				val[i][j][k] = STATION_UNDETERMINED;
			}
		}
	}

	int* transmitSequence = new int[stationCount];
	for(int i = 0; i < stationCount; i ++)
		transmitSequence[i] = i;
	this->sort(transmitSequence, SORT_BY_X_AXIS);
	this->leftMost = transmitSequence[0];
	this->recursivelyDPNoRecomp(0, stationCount - 1,
								transmitSequence,
								this->source + 1,
								val, keyP, keyQ);
	this->constructScheduleNoRecomp(0, stationCount - 1,
									transmitSequence,
									this->source + 1,
									val, keyP, keyQ);
	schedule.setParent(transmitSequence[0], STATION_SOURCE);	
	for(int i = 0; i < stationCount; i ++) {
		for(int j = 0; j < stationCount; j ++) {
			delete[] keyP[i][j];
			delete[] keyQ[i][j];
			delete[] val[i][j];
		}
		delete[] keyP[i];
		delete[] keyQ[i];
		delete[] val[i];
	}
	delete[] keyP;
	delete[] keyQ;
	delete[] val;
	delete[] transmitSequence;
}

long BroadcastScheme::recursivelyDPNoRecomp(int startPos, int endPos, int* sequence, int m, long*** val, int*** keyP, int*** keyQ) {
	/*int interval = endPos - startPos;
	if(m > interval) {
		if(val[startPos][endPos][interval] == STATION_UNDETERMINED) {
			this->recursivelyDPNoRecomp(startPos, endPos, 
										sequence, interval,
										val, keyP, keyQ);
		}
		keyP[startPos][endPos][m] = keyP[startPos][endPos][interval];
		keyQ[startPos][endPos][m] = keyQ[startPos][endPos][interval];
		val[startPos][endPos][m] = val[startPos][endPos][interval]; 
		cout << startPos << ',' << endPos << ',' 
					 	 << m << " = " << val[startPos][endPos][m]
					 	 << " via " << keyP[startPos][endPos][m] 
					 	 << " firstHalf = " << keyQ[startPos][endPos][m] 
					 	 << endl;
		return val[startPos][endPos][m];
	}*/
	m = min(endPos - startPos, m);
	if(val[startPos][endPos][m] != STATION_UNDETERMINED)
		return val[startPos][endPos][m];
	int start = sequence[startPos];
	int end = sequence[endPos];
	if(m == 1) {
		keyP[startPos][endPos][1] = startPos;
		keyQ[startPos][endPos][1] = 1;
		val[startPos][endPos][1] = stations[start].getXDiff(stations[end]); 
		return val[startPos][endPos][1];
	}

	int pivot = startPos,
		firstHalf = 1;
	long power1, power2, totalPower, minPower;
		// direct transmission is considered in particular
		if(val[startPos][endPos][1] == STATION_UNDETERMINED) {
		 	minPower = stations[start].getXDiff(stations[end]);
			val[startPos][endPos][1] = minPower;
			keyP[startPos][endPos][1] = startPos;
			keyQ[startPos][endPos][1] = 1;
		} else minPower = val[startPos][endPos][1];
	for(int i = 1; i < m; i ++) {
		// ensure 1st and 2nd half has at least
		// one transmision therefore no infinite
		// recursion
		for(int j = startPos + 1; j < endPos; j ++) {
			// transmission is always from left to right
			power1 = this->recursivelyDPNoRecomp(startPos, j, 
												 sequence, i,
												 val, keyP, keyQ);
			power2 = this->recursivelyDPNoRecomp(j, endPos,
												 sequence, m - i,
												 val, keyP, keyQ);
			totalPower = power1 + power2;
			if(totalPower < minPower) {
				minPower = totalPower;
				pivot = j;
				firstHalf = i;
			}
		}
	}
	val[startPos][endPos][m] = minPower;
	keyP[startPos][endPos][m] = pivot;
	keyQ[startPos][endPos][m] = firstHalf;
/*	cout << startPos << ',' << endPos << ',' 
					 << m << " = " << val[startPos][endPos][m]
					 << " via " << keyP[startPos][endPos][m] 
					 << " firstHalf = " << keyQ[startPos][endPos][m] 
					 << endl;*/
	
	return val[startPos][endPos][m];
}

void BroadcastScheme::dynamicProgIter() {
	int*** keyP = new int**[stationCount];
	int*** keyQ = new int**[stationCount];
	long*** val = new long**[stationCount];
	for(int i = 0; i < stationCount; i ++) {
		keyP[i] = new int*[stationCount];
		keyQ[i] = new int*[stationCount];
		val[i] = new long*[stationCount];
		for(int j = 0; j < stationCount; j ++) {
			keyP[i][j] = new int[this->source + 2];
			keyQ[i][j] = new int[this->source + 2];
			val[i][j] = new long[this->source + 2];
		}
	}
	int start, end, pivot, firstHalf;
	long minPower, crtPower;
	int* transmitSequence = new int[stationCount];
	for(int i = 0; i < stationCount; i ++)
		transmitSequence[i] = i;
	this->sort(transmitSequence, SORT_BY_X_AXIS);
	this->leftMost = transmitSequence[0];
	for(int m = 1; m < this->source + 2; m ++) {
		// m means # of intermediate transmitters
		for(int j = 1; j < stationCount; j ++) {
			// j means # of intervals
			for(int i = 0; i + j < stationCount; i ++) {
				// i means start position
				start = transmitSequence[i];
				end = transmitSequence[i + j];
				if(m == 1) {
					val[i][i + j][1] = stations[start].
						getXDiff(stations[end]);
					keyP[i][i + j][1] = i;
					keyQ[i][i + j][1] = 1;
				} else if(m > j) {
					val[i][i + j][m] = val[i][i + j][j];
					keyP[i][i + j][m] = keyP[i][i + j][j];
					keyQ[i][i + j][m] = keyQ[i][i + j][j];
				} else {
					minPower = val[i][i + j][1];
					pivot = i;
					firstHalf = 1;
					for(int l = 1; l < m; l ++) {
						// l means # of intermediate transmitters
						// in first half
						for(int n = 1; n < j; n ++) {
							// n means # of intervals
							crtPower = val[i][i + n][l]
										+ val[i + n][i + j][m - l];
							if(crtPower < minPower) {
								minPower=  crtPower;
								pivot = i + n;
								firstHalf = l;
							}
						}
					}
					val[i][i + j][m] = minPower;
					keyP[i][i + j][m] = pivot;
					keyQ[i][i + j][m] = firstHalf;
				}
			/*	cout << i << ',' << i + j << ',' 
					 << m + 1 << " = " << val[i][i + j][m]
					 << " via " << keyP[i][i + j][m] 
					 << " firstHalf = " << keyQ[i][i + j][m] 
					 << endl;*/
			}
		}
	}
	this->constructScheduleNoRecomp(0, stationCount - 1,
									transmitSequence, this->source + 1,
									val, keyP, keyQ);
	schedule.setParent(transmitSequence[0], STATION_SOURCE);
	for(int i = 0; i < stationCount; i ++) {
		for(int j = 0; j < stationCount; j ++) {
			delete[] keyP[i][j];
			delete[] keyQ[i][j];
			delete[] val[i][j];
		}
		delete[] keyP[i];
		delete[] keyQ[i];
		delete[] val[i];
	}
	delete[] keyP;
	delete[] keyQ;
	delete[] val;
	delete[] transmitSequence;
}

void BroadcastScheme::printScheme(bool mode) const {
	int parent;
	cout << schedule.calculateTotalPower() << endl;
	for(int i = 0; i < stationCount; i ++) {
		parent = schedule.getParent(i);
		if(mode)
			parent = (i == this->source) ? parent : parent + 1;
		else
			parent = (i == this->leftMost) ? parent : parent + 1;
		if(schedule.getPower(i) != 0)
			cout << i + 1 << ' '
				 << schedule.getPower(i) << ' ' 
				 << parent << endl;
	}
}
