"""
This script is the same as compute Temporal Average script, but only to compute
the average value of pressure and velocity.
This is for Coronary Validation Pipeline developed in Ana's project.
"""
import os
import argparse
from glob import glob
from utilities import *
import numpy as np
from scipy.fftpack import fftfreq,fft,ifft

class ComputeTemporalAverage():
	def __init__(self,Args):
		self.Args=Args
		if self.Args.OutputFolder is None:
			self.Args.OutputFolder=self.Args.InputFolder+"/../TemporalAvg"
			os.system("mkdir %s"%self.Args.OutputFolder)
		
	def Main(self):
		#Read all of the file name
		InputFiles1=sorted(glob(self.Args.InputFolder+"/*.vtu")) #volumetric files 
        #InputFiles2=sorted(glob(self.Args.InputFolder+"/*.vtp")) #surface files
		

		#Loop over all of the file names and store the values
		File1=ReadVTUFile(InputFiles1[0])
		#if len(InputFiles2)>0: File2=ReadVTPFile(InputFiles2[0])
		
		#Get the number of points
		NPoints1=File1.GetNumberOfPoints()
		#if len(InputFiles2)>0: NPoints2=File2.GetNumberOfPoints()
		
		#For volumetric data
		VelocityX=np.zeros(shape=(len(InputFiles1),NPoints1))
		VelocityY=np.zeros(shape=(len(InputFiles1),NPoints1))
		VelocityZ=np.zeros(shape=(len(InputFiles1),NPoints1))
		VelocityMag=np.zeros(shape=(len(InputFiles1),NPoints1))
		VelocityMag_=np.zeros(NPoints1)
		Pressure_=np.zeros(NPoints1)

		

		#Loop over SV result files
		counter=0 #filename or timestep
		for j in range(len(InputFiles1)):
			VelocityFile_=ReadVTUFile(InputFiles1[j])
			#if len(InputFiles2)>0: WSSFile_=ReadVTPFile(InputFiles2[j])
			
			print ("------ Looping over %s"%InputFiles1[j])
			#Loop over all of the points
			for i in range(0,NPoints1):
				VelocityX[counter,i]=VelocityFile_.GetPointData().GetArray("velocity").GetValue(i*3)	
				VelocityY[counter,i]=VelocityFile_.GetPointData().GetArray("velocity").GetValue(i*3+1)	
				VelocityZ[counter,i]=VelocityFile_.GetPointData().GetArray("velocity").GetValue(i*3+2)
				VelocityMag[counter,i]=np.sqrt(VelocityX[counter,i]**2+VelocityY[counter,i]**2+VelocityZ[counter,i]**2)
				VelocityMag_[i]+=np.sqrt(VelocityX[counter,i]**2+VelocityY[counter,i]**2+VelocityZ[counter,i]**2)
				Pressure_[i]+=VelocityFile_.GetPointData().GetArray("pressure").GetValue(i)

			
			counter+=1
			

	
		#N_ts=len(InputFiles1)
		#time=np.linspace(0,self.Args.Period,N_ts)
		VelocityMag_=VelocityMag_*(1./counter)
		Pressure_=Pressure_*(1./counter)


        #Add a new array to the Volumetric File
		VelocityMagVTK=numpy_to_vtk(VelocityMag_)
		VelocityMagVTK.SetName("Velocity")
		File1.GetPointData().AddArray(VelocityMagVTK)

		PressureVTK=numpy_to_vtk(Pressure_)
		PressureVTK.SetName("Pressure")
		File1.GetPointData().AddArray(PressureVTK)


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
		#if len(InputFiles2)>0:
		#	File2.GetPointData().RemoveArray("pressure")
		#	File2.GetPointData().RemoveArray("velocity")
		#	File2.GetPointData().RemoveArray("vinplane_traction")
		#	File2.GetPointData().RemoveArray("displacement")
		#	File2.GetPointData().RemoveArray("wallproperty")
		#	File2.GetPointData().RemoveArray("vWSS")
		#	File2.GetPointData().RemoveArray("timeDeriv")
		#	File2.GetPointData().RemoveArray("average_speed")
		#	File2.GetPointData().RemoveArray("average_pressure")	
		

		#Write the vtu file
		WriteVTUFile(self.Args.OutputFolder+"/TemporalVolumetricAveragedResults.vtu",File1)
		#if len(InputFiles2)>0:
		#	WriteVTPFile(self.Args.OutputFolder+"/TemporalSurfaceAveragedResults.vtp",File2)


	
if __name__=="__main__":
        #Arguments
	parser= argparse.ArgumentParser(description="This script will average over all of the vtu file provided in the results folder")

	parser.add_argument('-InputFolder', '--InputFolder', type=str, required=True, dest="InputFolder", help="The input folder that contains all of the results file, taged ass all_results.vtu.XXXXX.vtu")	

	parser.add_argument('-OutputFolder', '--OutputFolder', type=str, required=False, dest="OutputFolder", help="The output folder to store the time-averaged file in.")
	
	#Put all the arguments together
	args=parser.parse_args()

    #Call your Class
	ComputeTemporalAverage(args).Main()	
