#include "BroadcastScheme.h"
#include "List.h"
#include "iostream"
using namespace std;
int main(int argc, char* args[]) {
	int method, stationCount, source;
	cin >> method >> stationCount >> source;
	BroadcastScheme broadcastScheme(method, stationCount, source - 1);
	broadcastScheme.generateSchedule();
	if(method > METHOD_OWN)
		broadcastScheme.printScheme(SOURCE_LEFTMOST);
	else broadcastScheme.printScheme(SOURCE_INPUT);
}
