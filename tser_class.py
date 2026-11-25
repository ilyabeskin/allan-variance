import numpy as np

class TSer:
	def __init__(self, arr, t0=0, dT=1, unit=None):
		self.arr=arr
		self.t0=t0
		self.dT=dT
		self.unit=unit

	def ind_to_time(self,i):
		#Checks for indexes less than 0 or longer than the array.
		indVal = min( max( 0,int(i) ) , len(self.arr)-1 )
		return indVal*self.dT+self.t0

	def time_to_ind(self,t):
		indVal = int(np.round((t-self.t0)/self.dT))
		#Checks for indexes less than 0 or longer than the array.
		return min( max( 0,indVal ) , len(self.arr)-1 )

	def units(self):
		print(self.unit)
		return self.units

	def concat(self, newTSer):
		#Add check that time series have same dT
		self.arr = np.concatenate((self.arr, newTSer.arr))

	def crop(self, strt, end, method="time"):
		if method=="time":
			strtInd = self.time_to_ind(strt)
			endInd = self.time_to_ind(end)
			self.t0 = self.ind_to_time(strtInd)
			self.arr = self.arr[strtInd:endInd]
			
		elif method=="ind":
			self.t0 = self.ind_to_time(strt)
			self.arr = self.arr[strt:end]

		else:
			raise Exception("Unrecognized Method. TSer obj. unchanged.")

	def dupl(self):
		tSerCopy = TSer(self.arr, self.t0, self.dT, unit=self.unit)
		return tSerCopy