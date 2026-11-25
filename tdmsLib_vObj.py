from nptdms import TdmsFile
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from tser_class import TSer

def readFile(path, fName, num, grName):

	if num<10:
		numStr = "000"+str(num)
	elif num<100:
		numStr = "00"+str(num)
	elif num<1000:
		numStr = "0"+str(num)
	else:
		numStr = str(num)

	if isinstance(path, Path):
		fullPath = path / str(fName+grName+numStr+".tdms")
	elif isinstance(path, str):
		fullPath = path+fName+grName+numStr+".tdms"
	else:
		raise TypeError("ERROR: path can only be a Path obj. or str.")

	with TdmsFile.open(fullPath) as tdms_file:
		all_groups = tdms_file.groups()
		#print(all_groups)
		#print(all_groups[0])
		#all_channels = tdms_file['_counts_'].channels()
		#print(all_channels)
		chData = tdms_file[grName]['Channel 0'][:]
		dT=tdms_file[grName]['Channel 0'].properties["wf_increment"]
		t0=tdms_file[grName]['Channel 0'].properties["wf_start_offset"]

	tSerObj = TSer(chData, t0, dT, unit="sec")

	return tSerObj

def readConsecFiles(path, fName, grName, numStrt, numEnd):
	nums = nums=list(range(numStrt, numEnd+1))
	
	tSerFull = readFile(path, fName, nums[0], grName)

	for i in nums[1:]:
		tSerBlock = readFile(path, fName, i, grName)
		tSerFull.concat(tSerBlock)

	return tSerFull

def blockData_sum(blockSz, tSer):
	blockSz=int(blockSz)

	finLen = int(np.floor(len(tSer.arr)/blockSz))
	cropLen = int(finLen*blockSz)
	blockedData = np.sum(np.reshape(tSer.arr[:cropLen],(-1,blockSz)), axis=1)

	tSer_blocked = TSer(blockedData, tSer.t0+tSer.ind_to_time(blockSz/2), blockSz*tSer.dT)
	
	return tSer_blocked, len(tSer.arr)-finLen*blockSz

def blockData_ave(blockSz, tSer):
	blockSz=int(blockSz)

	finLen = int(np.floor(len(tSer.arr)/blockSz))
	cropLen = int(finLen*blockSz)
	blockedData = np.mean(np.reshape(tSer.arr[:cropLen],(-1,blockSz)), axis=1)

	tSer_blocked = TSer(blockedData, tSer.t0+tSer.ind_to_time(blockSz/2), blockSz*tSer.dT)
	
	return tSer_blocked, len(tSer.arr)-finLen*blockSz

def prepTimeTraceAx(ax, tSer, strtT=None, endT=None):
	strtInd = 0
	endInd = len(tSer.arr)

	if strtT != None:
		strtInd = tSer.time_to_ind(strtT)
	if endT != None:
		endInd = tSer.time_to_ind(endT)
	
	timeArr = np.arange(0,len(tSer.arr)*tSer.dT, tSer.dT)+tSer.t0
	ax.plot(timeArr[strtInd:endInd], tSer.arr[strtInd:endInd])

def dispTimeTrace(tSer, strtT=None, endT=None):
	fig, ax = plt.subplots()
	
	prepTimeTraceAx(ax, tSer, strtT, endT)
	
	ax.set_ylabel("Signal")
	ax.set_xlabel("Time [s]")
	plt.show()

def dispScatFluorVolt(tSer_scat=TSer(np.array([]), 0 , 1), tSer_fluor=TSer(np.array([]), 0 , 1), tSer_volt=TSer(np.array([]), 0 , 1), strtT=None, endT=None):
	#Figure out which data sets have been given.
	numAx = (len(tSer_scat.arr)>0)+(len(tSer_fluor.arr)>0)+(len(tSer_volt.arr)>0)

	fig, axs = plt.subplots(numAx, sharex=True, squeeze=False)

	curAx = int(0)
	if len(tSer_scat.arr)>0:
		prepTimeTraceAx(axs[curAx][0], tSer_scat, strtT, endT)
		axs[curAx][0].set_xlabel("Time [s]")
		axs[curAx][0].set_ylabel("Scattering [V]")
		curAx = int(curAx+1)

	if len(tSer_fluor.arr)>0:
		prepTimeTraceAx(axs[curAx][0], tSer_fluor, strtT, endT)
		axs[curAx][0].set_xlabel("Time [s]")
		axs[curAx][0].set_ylabel(f"Fluorescence \n[counts / {tSer_fluor.dT*1e3}ms]")
		curAx = int(curAx+1)

	if len(tSer_volt.arr)>0:
		prepTimeTraceAx(axs[curAx][0], tSer_volt, strtT, endT)
		axs[curAx][0].set_xlabel("Time [s]")
		axs[curAx][0].set_ylabel("Input Voltage [V]")

	plt.show()	


def compVar(scatArr , blockSize):
	
	corTime = 1e-5
	corLen = int(1e-5/(scatArr.dT)) # Complete guess of 0.1ms for correlation time.
	var = np.var( scatArr.arr , ddof=1)

	n_eff = len(scatArr.arr)/(1+2*corLen)
	var_delta =  var*( np.sqrt(2/(n_eff-1)) )

	return var, var_delta

def compAllanVar(tSer):
	arr1 = tSer.arr[0:len(tSer.arr)-1]
	arr2 = tSer.arr[1:]
	allanVar = np.sum(((arr2-arr1)/2)**2)/(len(arr1)-1)
	allanVar_delta = allanVar/np.sqrt((len(arr1))) 
	#https://liquidinstruments.com/application-notes/understanding-and-performing-allan-variance-measurements/
	#P. Lesage and C. Audoin, 1973, “Characterization of frequency stability: uncertainty due to the finite number of measurements,” IEEE Transactions on Instrumentation and Measurement, IM-22, 157–161.
	return allanVar, allanVar_delta

def findDevPoint(tSer, strtT=None, endT=None, method="scat"):
	devSigma = 3 #max allowed deviation

	if method=="scat":	
		blLen = int(100*1E2) #number of points in baseline
		gapLen = int(10*1E2) #number of points to skip
		checkLen = int(10) #number of points to check
	if method=="fluor":
		devSigma = 4 #max allowed deviation
		blLen = int(20) #number of points in baseline
		gapLen = int(5) #number of points to skip
		checkLen = int(2) #number of points to check

	if strtT==None:
		strtT = 0

	if endT==None:
		endT = len(tSer.arr)

	tSerCrop = tSer.dupl()
	tSerCrop.crop(strtT, endT, method="time")

	devInd = None

	ind=-1
	while devInd == None and ind<len(tSerCrop.arr)-(blLen+gapLen+checkLen):
		ind = ind+1
		meanVal = np.mean(tSerCrop.arr[ind:ind+blLen])
		stdDev = np.std(tSerCrop.arr[ind:ind+blLen], ddof=1)
		if np.abs(np.mean(tSerCrop.arr[ind+gapLen+blLen:ind+gapLen+blLen+checkLen])-meanVal) > devSigma*stdDev :
			devInd = ind+blLen+gapLen

	return devInd, meanVal, stdDev/np.sqrt(checkLen)

def findStep(tSer, strtT=None, endT=None, method="scat", graph=False):

	stepStrtT = None
	strtVal = None
	strtVal_delta = None

	stepEndT = None
	endVal = None
	endVal_delta = None

	#Look for start of transition:
	tSerCrop = tSer.dupl()
	tSerCrop.crop(strtT, endT, method="time")
	strtInd, strtVal, strtVal_delta = findDevPoint(tSerCrop, method=method)
	stepStrtT = tSerCrop.ind_to_time(strtInd)

	#Look for end of transition:
	tSerCropInv = TSer(tSerCrop.arr[::-1], 0, 1) #Note: t0 and dT are not set correctly here because dT can't be negative.
	endInd_inv, endVal, endVal_delta = findDevPoint(tSerCropInv, method=method)
	endInd = len(tSerCrop.arr)-endInd_inv
	stepEndT = tSerCrop.ind_to_time(endInd)

	transLen = stepEndT-stepStrtT

	if graph==True:
		fig, ax = plt.subplots()
		prepTimeTraceAx(ax, tSer, strtT=stepStrtT-3*transLen, endT=stepEndT+3*transLen)
		ax.set_xlabel("Time [s]")
		ax.set_ylabel("Signal")
		ax.hlines([strtVal,endVal],
			xmin=[stepStrtT-3*transLen, stepEndT],
			xmax=[stepStrtT, stepEndT+3*transLen],
			colors='black', linewidths=3)
		plt.show()

	return stepStrtT, strtVal, strtVal_delta, stepEndT, endVal, endVal_delta