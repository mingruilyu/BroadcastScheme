import os
import time
testSeries1 = 'testcases/testseries1/'
testSeries2 = 'testcases/testseries2/'
binDir = 'broadcast/'
resultDir = 'testresults/'
progName = 'prog'
def check2():
	files = os.listdir(testSeries2)
	for file in files:
		if file.startswith('inp_8_200'):
			if file[4] == '6':
				if int(file[6 : ]) > 20:
					print('Skipping file ' + file)
					continue
			startTime = time.time()
			print('Testing file ' + file)
			os.system(binDir + progName + ' < ' + testSeries2 + file + ' > ' + resultDir + 'out' + file[3 : ])
			endTime = time.time()
			timeElapse = (endTime - startTime) * 1000 # ms
			print('Testing ' + file + ' used ' + str(int(timeElapse)) + ' ms ')
			#os.system('diff output ' + testSeries1 + 'out' + file[3:])
			
def check1():
	files = os.listdir(testSeries1)
	for file in files:
		if file.startswith('inp_1'):
			startTime = time.time()
			print('Testing ' + file)
			os.system(binDir + progName + ' < ' + testSeries1 + file + ' > ' + resultDir + 'out' + file[3 : ])
			endTime = time.time()
			timeElapse = (endTime - startTime) * 1000
			print('Testing ' + file + ' used ' + str(int(timeElapse)) + ' ms ')
			#os.system('diff output ' + testSeries2 + 'out' + file[3:])
			
if __name__ == '__main__':
	check2()
	
