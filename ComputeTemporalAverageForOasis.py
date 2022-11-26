import os
import argparse
import sys
from glob import glob
from utilities import *
import numpy as np
from scipy.fftpack import fftfreq,fft,ifft

class ComputeTemporalAverage():
	def __init__(self,Args):
		self.Args=Args

	def Main(self):
		#Read all of the file name
		InputFiles=sorted(glob(self.Args.InputFolder+"/all_results.vtu*.vtu"))
		InputResults=sorted(glob(self.Args.ProcessedResults+"/*.vtu"))
		
		os.system("mkdir %s"%self.Args.OutputFolder)
 
		#Loop over all of the file names and store the values
		File1=ReadVTUFile(InputFiles[0])
		NPoints=File1.GetNumberOfPoints()

		#Initialize file for surface parameters
		for file in InputResults:
			File2=ReadVTUFile(file)
			
			if File2.GetNumberOfPoints()<NPoints:
				File2=ReadVTUFile(file)
				break

		VelocityX=np.zeros(shape=(len(InputFiles),NPoints))
		VelocityY=np.zeros(shape=(len(InputFiles),NPoints))
		VelocityZ=np.zeros(shape=(len(InputFiles),NPoints))
		VelocityMag=np.zeros(shape=(len(InputFiles),NPoints))
		VelocityMag_=np.zeros(NPoints)
		#Loop over processed results files
		for file in InputResults:
			File_=ReadVTUFile(file)
			
			if File_.GetPointData().GetArrayName(0) is None:
				continue

			else:   
				dataset_name=File_.GetPointData().GetArrayName(0)
				data=File_.GetPointData().GetArray(dataset_name)
				print ("--- Looping over %s"%file)
					
                        	#Add to array
				dataVTK=numpy_to_vtk(data)
				name=file.split('/')[1]
				name=name.split('0')[0]
				dataVTK.SetName(name)

				if data.GetNumberOfValues()<NPoints:
					File2.GetPointData().AddArray(dataVTK)
				else:
					File1.GetPointData().AddArray(dataVTK) 

		#Loop over SV result files
		counter=0 #filename or timestep
		for InputFile_ in InputFiles:
			VelocityFile_=ReadVTUFile(InputFile_)
			print ("--- Looping over %s"%InputFile_)
			#Loop over all of the points
			for i in range(0,NPoints):
				VelocityX[counter,i]=VelocityFile_.GetPointData().GetArray("velocity").GetValue(i*3)	
				VelocityY[counter,i]=VelocityFile_.GetPointData().GetArray("velocity").GetValue(i*3+1)	
				VelocityZ[counter,i]=VelocityFile_.GetPointData().GetArray("velocity").GetValue(i*3+2)
			
				VelocityMag[counter,i]=np.sqrt(VelocityX[counter,i]**2+VelocityY[counter,i]**2+VelocityZ[counter,i]**2)
				VelocityMag_[i]+=np.sqrt(VelocityX[counter,i]**2+VelocityY[counter,i]**2+VelocityZ[counter,i]**2)
			counter+=1
			
		
		N_ts=len(InputFiles)
		time=np.linspace(0,self.Args.Period,N_ts)
		self.W = fftfreq(N_ts, d=time[1]-time[0])
		VelocityMag_=VelocityMag_*(1./counter)

		#Compute Frequency Analysis
		SPI=np.zeros(NPoints)
		TKE=np.zeros(NPoints)

		for i in range(NPoints):
			SPI[i]=self.filter_SPI(VelocityMag[:,i])     
			U_tke_=np.mean(np.power(self.filter_TKE(VelocityX[:,i]),2))
			V_tke_=np.mean(np.power(self.filter_TKE(VelocityY[:,i]),2))
			W_tke_=np.mean(np.power(self.filter_TKE(VelocityZ[:,i]),2))
			TKE[i]=0.5*(U_tke_+V_tke_+W_tke_)

                #Add a new array to the file
		VelocityMagVTK=numpy_to_vtk(VelocityMag_)
		VelocityMagVTK.SetName("velocity_mag_average")
		File1.GetPointData().AddArray(VelocityMagVTK)

		SPI_VTK=numpy_to_vtk(SPI)
		SPI_VTK.SetName("SPI")
		File1.GetPointData().AddArray(SPI_VTK)

		TKE_VTK=numpy_to_vtk(TKE)
		TKE_VTK.SetName("TKE")
		File1.GetPointData().AddArray(TKE_VTK)

                #Loop over all of the Boundary Files 
		#File1.GetCellData().RemoveArray("GlobalElementID")
		#File1.GetPointData().RemoveArray("GlobalNodeID")
		#File1.GetPointData().RemoveArray("pressure")
		File1.GetPointData().RemoveArray("velocity")
		File1.GetPointData().RemoveArray("vinplane_traction")
		File1.GetPointData().RemoveArray("displacement")
		File1.GetPointData().RemoveArray("wallproperty")
		File1.GetPointData().RemoveArray("vWSS")
		File1.GetPointData().RemoveArray("timeDeriv")
		File1.GetPointData().RemoveArray("average_speed")
		#File1.GetPointData().RemoveArray("average_pressure")
		
		#Write the vtu file
		WriteVTUFile(self.Args.OutputFolder+"/TemporalVolumetricAveragedResults.vtu",File1)
		WriteVTUFile(self.Args.OutputFolder+"/TemporalSurfaceAveragedResults.vtu",File2)	

	def filter_TKE(self,U):
		#For FKE
		U_fft       = fft(U)
		U_cut_fft   = U_fft.copy()
		U_cut_fft[(self.W<self.Args.CutoffFrequency)]=0
		U_ifft      =ifft(U_cut_fft)
		return U_ifft

	def filter_SPI(self,U):
		#for HI
		U_fft       = fft(U-np.mean(U))
		#Cut off the upper frequency to 25Hz
		U_fft_25Hz   = U_fft.copy()
		U_fft_25Hz[(self.W<self.Args.CutoffFrequency)]=0
		#Cut off the lower frequency to 0Hz
		U_fft_0Hz =U_fft.copy()
		U_fft_0Hz[(self.W<0)]=0
		#Compute the absolute value
		Power_25Hz    =np.sum ( np.power( np.absolute(U_fft_25Hz),2))
		Power_0Hz     =np.sum ( np.power( np.absolute(U_fft_0Hz) ,2))
		if Power_0Hz<1e-5: return 0
		else:              return Power_25Hz/Power_0Hz
	
if __name__=="__main__":
        #Arguments
	parser= argparse.ArgumentParser(description="This script will average over all of the vtu file provided in the results folder")

	parser.add_argument('-InputFolder', '--InputFolder', type=str, required=True, dest="InputFolder", help="The input folder that contains all of the results file, taged ass all_results.vtu.XXXXX.vtu")	

	parser.add_argument('-ProcessedResults', '--ProcessedResults', type=str, required=True, dest="ProcessedResults", help="The input folder that contains all of the post processed volumetric and surface results file.")
        
	parser.add_argument('-OutputFolder', '--OutputFolder', type=str, required=True, dest="OutputFolder", help="The output folder to store the time-averaged file in.")
	parser.add_argument('-Period', '--Period', type=float, required=True, dest="Period",help="The duration of the cardiac cycle in seconds.")
	parser.add_argument('-CutoffFrequency', '--CutoffFrequency', type=int, required=False, default=25, dest="CutoffFrequency",help="The cut-off frequency to compute frequency-based biomarkers.")
	
	#Put all the arguments together
	args=parser.parse_args()

        #Call your Class
	ComputeTemporalAverage(args).Main()	
