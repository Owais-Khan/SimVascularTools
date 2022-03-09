import os
import argparse
import sys
from glob import glob
from utilities import *
import numpy as np
from numpy.fft import fftfreq, fft,ifft 
class ComputeSurfaceHemodynamics():
	def __init__(self,Args):
		self.Args=Args
		self.low_cut=25
		if self.Args.OutputFileName is None:
			self.Args.OutputFileName=self.Args.InputFolder+"../SurfaceHemodynamics.vtp"	
	def Main(self):
		#Get all of the file names
		FileNames=sorted(glob(self.Args.InputFolder+"/all_results*.vtp"))

		#Read a Surface File to Get Geometry Data
		Surface=ReadVTPFile(FileNames[0])
		Surface.GetPointData().RemoveArray("GlobalNodeID")
		Surface.GetPointData().RemoveArray("average_pressure")
		Surface.GetPointData().RemoveArray("average_speed")
		Surface.GetPointData().RemoveArray("pressure")
		Surface.GetPointData().RemoveArray("timeDeriv")
		Surface.GetPointData().RemoveArray("vWSS")
		Surface.GetPointData().RemoveArray("velocity")
		Surface.GetPointData().RemoveArray("vinplane_traction")
		
		#Get the WSS and Pressure Data
		print ("--- Reading the WSS Files")
		WSS1,WSS2,WSS3,WSSmag,Pressure=self.ReadData(FileNames)
	
		#Compute Average Pressure
		print ("--- Computing Average Pressure in mmHg")
		Pressure=np.average(Pressure,axis=0)/1333.22
	
		#Time averaged WSS
		print ("--- Computing Time-average Wall Shear Stress (TA-WSS)")
		TAWSS=np.average(WSSmag,axis=0)

		#Compute OSI
		print ("--- Computing Oscillatory Shear Index (OSI)")
		WSS1avg=np.average(WSS1,axis=0)	
		WSS2avg=np.average(WSS2,axis=0)	
		WSS3avg=np.average(WSS3,axis=0)
		OSI=np.zeros(len(WSS1avg))
		for i in range(len(WSS1avg)):
			OSI[i]=0.5*(1-np.sqrt(WSS1avg[i]**2+WSS2avg[i]**2+WSS3avg[i]**2)/TAWSS[i])	
		#Compute Relative Residence Time
		print ("--- Computing Relative Residence Time (RRT)")
		RRT=1/((1-2*OSI)*TAWSS)	
	
		#Compute SPI
		print ("--- Compute Spectral Power Index (SPI)")
		#Define the ferquency array
		N_ts=len(WSS1)
		time=np.linspace(0,self.Args.Period,N_ts)
		self.W = fftfreq(N_ts, d=time[1]-time[0])
		self.dt= time[1]-time[0]
		#Initialize the array
		SPI=np.zeros(self.NPoints)
		for i in range(self.NPoints):
			SPI[i]=self.FilterSPI(WSSmag[:,i])
		


		#Add all of the Arrays
		self.AddArray("TAWSS",TAWSS,Surface)
		WriteVTPFile(self.Args.OutputFileName,Surface)	
		

     
	def AddArray(self,ArrayName,Array,Surface):
		NewArray=vtk.vtkFloatArray()
		NewArray.SetNumberOfComponents(1)
		NewArray.SetName(ArrayName)
                
		N=Surface.GetNumberOfPoints()
                
		for i in range(N):
			NewArray.InsertNextValue(Array[i])
		Surface.GetPointData().AddArray(NewArray)

	def FilterSPI(self,WSS_):
		WSS_fft=fft(WSS_-np.mean(WSS_))
		#Cut off the upper frequency to 25Hz
		WSS_fft_25Hz   = WSS_fft.copy()
		WSS_fft_25Hz[(self.W<self.low_cut)]=0
		#Cut off the lower frequency to 0Hz
		WSS_fft_0Hz =WSS_fft.copy()
		WSS_fft_0Hz[(self.W<0)]=0
		#Compute the absolute value
		Power_25Hz    =np.sum ( np.power( np.absolute(WSS_fft_25Hz),2))
		Power_0Hz     =np.sum ( np.power( np.absolute(WSS_fft_0Hz) ,2))
		if Power_0Hz<1e-5: return 0
		else:              return Power_25Hz/Power_0Hz


	
	def ReadData(self,FileNames):
		print ("------ Number of Files Detected: %d"%len(FileNames))
		FileNames=FileNames[0:5]
		NFiles=len(FileNames)
		for i in range(0,len(FileNames)):
			print ("------Loading: %s"%FileNames[i])

			#Read the input files
			SurfaceFile_=ReadVTPFile(FileNames[i])
			
			#Get the number of points and initialize arrays
			if i==0: 
				self.NPoints=SurfaceFile_.GetNumberOfPoints()
				WSS1=np.zeros(shape=(NFiles,self.NPoints))
				WSS2=np.zeros(shape=(NFiles,self.NPoints))
				WSS3=np.zeros(shape=(NFiles,self.NPoints))
				WSSmag=np.zeros(shape=(NFiles,self.NPoints))
				Pressure=np.zeros(shape=(NFiles,self.NPoints))

			#Copy over the Wall Shear Stress data
			k=0
			for j in range(0,3*self.NPoints,3):
				WSS1[i,k]=SurfaceFile_.GetPointData().GetArray("vWSS").GetValue(j)	
				WSS2[i,k]=SurfaceFile_.GetPointData().GetArray("vWSS").GetValue(j+1)	
				WSS3[i,k]=SurfaceFile_.GetPointData().GetArray("vWSS").GetValue(j+2)	
				WSSmag[i,k]=np.sqrt(WSS1[i,k]**2+WSS2[i,k]**2+WSS3[i,k]**2)
				Pressure[i,k]=SurfaceFile_.GetPointData().GetArray("pressure").GetValue(k)
				k+=1

		return WSS1,WSS2,WSS3,WSSmag,Pressure
				
			
if __name__=="__main__":
        #Arguments
	parser= argparse.ArgumentParser(description="This script will compute hemodynamic indices such as WSS and OSI.")

        #Input filename for the mesh-complete folder
	parser.add_argument('-InputFolder', '--InputFolder', type=str, required=True, default="ResultsFolder", dest="InputFolder", help="The folder that contains the VTP and VTU files.")
	
	#Write the Output FileName
	parser.add_argument('-OutputFileName', '--OutputFileName', type=str, required=False, dest="OutputFileName", help="The output filename")
        
	parser.add_argument('-Period', '--Period', type=float, required=True, dest="Period", help="Period of the cardiac cycle in seconds.")

        #Put all the arguments together
	args=parser.parse_args()

        #Call your Class
	ComputeSurfaceHemodynamics(args).Main() 
