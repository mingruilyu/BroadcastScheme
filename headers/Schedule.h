#ifndef SCHEDULE_H
#define SCHEDULE_H
class Schedule {
private:
	int* schedule;
	int* parent;
	int lastIndex;
	int previousValue;
	int count;
	bool modifyFlag;
public:
	Schedule(const int count);
	~Schedule();
	void modify(const int index,const int value);
	// return the power(dist ^ 2)
	void setPower(const int index, const long power);
	long getPower(const int index) const;
	int getParent(const int index) const;
	void rollback();
	void setParent(const int index, const int parent);
	long calculateTotalPower() const;
	void printSchedule() const;
};
#endif
