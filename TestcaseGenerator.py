import sys
import random
import os

stationCount = int(sys.argv[2])
distribution = 5000
testSeries1Dir = 'testcases/testseries1/'
testSeries2Dir = 'testcases/testseries2/'
def testSeries2():
	source = int(stationCount / 2)
	fileList = []
	for i in range(6, 9):
		filename = 'inp' + '_' + str(i) + '_' + str(stationCount)
		fileList.append(open(testSeries2Dir + filename, 'w'))
		header = str(i) + '\t' + str(stationCount) \
				 + '\t' + str(source) + '\n'
		fileList[i - 6].write(header)
		
	for n in range(0, stationCount):
		s = str(random.randrange(0, distribution)) + '\t' \
			+ str('0') + '\n'
		for i in range(6, 9):
			fileList[i - 6].write(s)

def testSeries1():
	source = random.randrange(0, stationCount)
	fileList = []
	for i in range(1, 6):
		filename = 'inp' + '_' + str(i) + '_' + str(stationCount)
		fileList.append(open(testSeries1Dir + filename, 'w'))
		header = str(i) + '\t' + str(stationCount) \
				 + '\t' + str(source) + '\n'
		fileList[i - 1].write(header)
		
	for n in range(0, stationCount):
		s = str(random.randrange(0, distribution)) + '\t' \
			+ str(random.randrange(0, distribution)) + '\n'
		for i in range(1, 6):
			fileList[i - 1].write(s)

if __name__ == '__main__':
	if sys.argv[1] == '1':
		testSeries1()
	elif sys.argv[1] == '2':
		testSeries2()



