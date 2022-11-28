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
		if self.Args.OutputFolder is None:
			self.Args.OutputFolder=self.Args.InputFolder+"../TemporalAvg"
			os.system("mkdir %s"%self.Args.OutputFolder)

	def Main(self):
		#Read all of the file name
		InputFiles1=sorted(glob(self.Args.InputFolder+"/all_results.vtu*.vtu")) #volumetric files 
		InputFiles2=sorted(glob(self.Args.InputFolder+"/all_results.vtp*.vtp")) #surface files
		
 
		#Loop over all of the file names and store the values
		File1=ReadVTUFile(InputFiles1[0])
		File2=ReadVTPFile(InputFiles2[0])
		
		#Get the number of points
		NPoints1=File1.GetNumberOfPoints()
		NPoints2=File2.GetNumberOfPoints()
		
		#For volumetric data
		VelocityX=np.zeros(shape=(len(InputFiles1),NPoints1))
		VelocityY=np.zeros(shape=(len(InputFiles1),NPoints1))
		VelocityZ=np.zeros(shape=(len(InputFiles1),NPoints1))
		VelocityMag=np.zeros(shape=(len(InputFiles1),NPoints1))
		VelocityMag_=np.zeros(NPoints1)

		#For surface data
		WSSX=np.zeros(shape=(len(InputFiles2),NPoints2))
		WSSY=np.zeros(shape=(len(InputFiles2),NPoints2))
		WSSZ=np.zeros(shape=(len(InputFiles2),NPoints2))
		WSSMag=np.zeros(shape=(len(InputFiles2),NPoints2))
		WSSMag_=np.zeros(NPoints2)
		

		#Loop over SV result files
		counter=0 #filename or timestep
		for j in range(len(InputFiles1)):
			VelocityFile_=ReadVTUFile(InputFiles1[j])
			WSSFile_=ReadVTPFile(InputFiles2[j])
			
			print ("------ Looping over %s"%InputFiles1[j])
			#Loop over all of the points
			for i in range(0,NPoints1):
				VelocityX[counter,i]=VelocityFile_.GetPointData().GetArray("velocity").GetValue(i*3)	
				VelocityY[counter,i]=VelocityFile_.GetPointData().GetArray("velocity").GetValue(i*3+1)	
				VelocityZ[counter,i]=VelocityFile_.GetPointData().GetArray("velocity").GetValue(i*3+2)
				VelocityMag[counter,i]=np.sqrt(VelocityX[counter,i]**2+VelocityY[counter,i]**2+VelocityZ[counter,i]**2)
				VelocityMag_[i]+=np.sqrt(VelocityX[counter,i]**2+VelocityY[counter,i]**2+VelocityZ[counter,i]**2)
			                        
			print ("------ Looping over %s"%InputFiles2[j])
			#Loop over all of the points
			for i in range(0,NPoints2):
				WSSX[counter,i]=WSSFile_.GetPointData().GetArray("vWSS").GetValue(i*3)
				WSSY[counter,i]=WSSFile_.GetPointData().GetArray("vWSS").GetValue(i*3+1)
				WSSZ[counter,i]=WSSFile_.GetPointData().GetArray("vWSS").GetValue(i*3+2)
				WSSMag[counter,i]=np.sqrt(WSSX[counter,i]**2+WSSY[counter,i]**2+WSSZ[counter,i]**2)
				WSSMag_[i]+=np.sqrt(WSSX[counter,i]**2+WSSY[counter,i]**2+WSSZ[counter,i]**2)

			counter+=1
			

	
		N_ts=len(InputFiles1)
		time=np.linspace(0,self.Args.Period,N_ts)
		self.W = fftfreq(N_ts, d=time[1]-time[0])
		VelocityMag_=VelocityMag_*(1./counter)
		WSSMag_=WSSMag_*(1./counter)

		#Compute Frequency Analysis
		SPI=np.zeros(NPoints1)
		TKE=np.zeros(NPoints1)
		SPI_surface=np.zeros(NPoints2)
		OSI=np.zeros(NPoints2)

		print ("\n")
		print ("--- Computing Volumetric SPI and TKE...")
		for i in range(NPoints1):
			SPI[i]=self.filter_SPI(VelocityMag[:,i])     
			U_tke_=np.mean(np.power(self.filter_TKE(VelocityX[:,i]),2))
			V_tke_=np.mean(np.power(self.filter_TKE(VelocityY[:,i]),2))
			W_tke_=np.mean(np.power(self.filter_TKE(VelocityZ[:,i]),2))
			TKE[i]=0.5*(U_tke_+V_tke_+W_tke_)

		print ("\n")
		print ("--- Computing Surface SPI, WSS and OSI ...")
		for i in range(NPoints2):
			SPI_surface[i]=self.filter_SPI(WSSMag[:,i])
			if WSSMag_[i]==0.0: 
				OSI[i]=0.0
			else:
				OSI[i]=0.5*(1-np.sqrt(np.mean(WSSX[:,i])**2+np.mean(WSSY[:,i])**2+np.mean(WSSZ[:,i])**2)/WSSMag_[i])	


                #Add a new array to the Volumetric File
		VelocityMagVTK=numpy_to_vtk(VelocityMag_)
		VelocityMagVTK.SetName("velocity_mag_average")
		File1.GetPointData().AddArray(VelocityMagVTK)

		SPI_VTK=numpy_to_vtk(SPI)
		SPI_VTK.SetName("SPI")
		File1.GetPointData().AddArray(SPI_VTK)

		TKE_VTK=numpy_to_vtk(TKE)
		TKE_VTK.SetName("TKE")
		File1.GetPointData().AddArray(TKE_VTK)

		#Add new array to the Surface File
		WSSMagVTK=numpy_to_vtk(WSSMag_)
		WSSMagVTK.SetName("WSS_mag_average")
		File2.GetPointData().AddArray(WSSMagVTK)

		SPI_surface_VTK=numpy_to_vtk(SPI_surface)
		SPI_surface_VTK.SetName("SPI_WSS")
		File2.GetPointData().AddArray(SPI_surface_VTK)
		
		OSI_VTK=numpy_to_vtk(OSI)
		OSI_VTK.SetName("OSI")
		File2.GetPointData().AddArray(OSI_VTK)
		

		#clean up the volumetric file
		File1.GetPointData().RemoveArray("pressure")
		File1.GetPointData().RemoveArray("velocity")
		File1.GetPointData().RemoveArray("vinplane_traction")
		File1.GetPointData().RemoveArray("displacement")
		File1.GetPointData().RemoveArray("wallproperty")
		File1.GetPointData().RemoveArray("vWSS")
		File1.GetPointData().RemoveArray("timeDeriv")
		File1.GetPointData().RemoveArray("average_speed")
		File1.GetPointData().RemoveArray("average_pressure")
	

                #Clean up the surface file
		File2.GetPointData().RemoveArray("pressure")
		File2.GetPointData().RemoveArray("velocity")
		File2.GetPointData().RemoveArray("vinplane_traction")
		File2.GetPointData().RemoveArray("displacement")
		File2.GetPointData().RemoveArray("wallproperty")
		File2.GetPointData().RemoveArray("vWSS")
		File2.GetPointData().RemoveArray("timeDeriv")
		File2.GetPointData().RemoveArray("average_speed")
		File2.GetPointData().RemoveArray("average_pressure")	
		
		#Write the vtu file
		WriteVTUFile(self.Args.OutputFolder+"/TemporalVolumetricAveragedResults.vtu",File1)
		WriteVTPFile(self.Args.OutputFolder+"/TemporalSurfaceAveragedResults.vtp",File2)

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

	parser.add_argument('-OutputFolder', '--OutputFolder', type=str, required=False, dest="OutputFolder", help="The output folder to store the time-averaged file in.")
	parser.add_argument('-Period', '--Period', type=float, required=True, dest="Period",help="The duration of the cardiac cycle in seconds.")
	parser.add_argument('-CutoffFrequency', '--CutoffFrequency', type=int, required=False, default=25, dest="CutoffFrequency",help="The cut-off frequency to compute frequency-based biomarkers.")
	
	#Put all the arguments together
	args=parser.parse_args()

        #Call your Class
	ComputeTemporalAverage(args).Main()	
