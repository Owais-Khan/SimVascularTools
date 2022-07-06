import os
import argparse
import sys
from glob import glob
from utilities import *
import numpy as np

class ComputeTemporalAverage():
	def __init__(self,Args):
		self.Args=Args

	def Main(self):
		#Read all of the file name
		InputFiles=sorted(glob(self.Args.InputFolder+"/all_results.vtu*.vtu"))

		#Loop over all of the file names and store the values
		File1=ReadVTUFile(InputFiles[0])
		NPoints=File1.GetNumberOfPoints()

		#Create a Velocity Array
		VelocityX=np.zeros(NPoints)
		VelocityY=np.zeros(NPoints)
		VelocityZ=np.zeros(NPoints)
		VelocityMag=np.zeros(NPoints)

		#Loop over all of the files
		counter=0
		for InputFile_ in InputFiles:
			VelocityFile_=ReadVTUFile(InputFile_)
			print ("--- Looping over %s"%InputFile_)
			#Loop over all of the points
			for i in range(0,NPoints):
				velx_=VelocityFile_.GetPointData().GetArray("velocity").GetValue(i*3)	
				vely_=VelocityFile_.GetPointData().GetArray("velocity").GetValue(i*3+1)	
				velz_=VelocityFile_.GetPointData().GetArray("velocity").GetValue(i*3+2)
				velmag_=np.sqrt(velx_**2+vely_**2+velz_**2)
				VelocityMag[i]+=velmag_
			counter+=1

		#Average the velocitymag
		VelocityMag=VelocityMag*(1./counter)

	
                #Add a new array to the file
		VelocityMagVTK=numpy_to_vtk(VelocityMag)
		VelocityMagVTK.SetName("velocity_mag_average")
		File1.GetPointData().AddArray(VelocityMagVTK)

                #Loop over all of the Boundary Files 
		File1.GetCellData().RemoveArray("GlobalElementID")
		File1.GetPointData().RemoveArray("GlobalNodeID")
		#File1.GetPointData().RemoveArray("pressure")
		File1.GetPointData().RemoveArray("velocity")
		File1.GetPointData().RemoveArray("vinplane_traction")
		File1.GetPointData().RemoveArray("displacement")
		File1.GetPointData().RemoveArray("wallproperty")
		File1.GetPointData().RemoveArray("vWSS")
		File1.GetPointData().RemoveArray("timeDeriv")
		File1.GetPointData().RemoveArray("average_speed")
		File1.GetPointData().RemoveArray("average_pressure")
		
		#Write the vtu file
		WriteVTUFile(self.Args.OutputFolder+"/TemporalAveragedResults.vtu",File1)	
	

	
if __name__=="__main__":
        #Arguments
	parser= argparse.ArgumentParser(description="This script will average over all of the vtu file provided in the results folder")

	parser.add_argument('-InputFolder', '--InputFolder', type=str, required=True, dest="InputFolder", help="The input folder that contains all of the results file, taged ass all_results.vtu.XXXXX.vtu")	
        
	parser.add_argument('-OutputFolder', '--OutputFolder', type=str, required=True, dest="OutputFolder", help="The output folder to store the time-averaged file in.")

        #Put all the arguments together
	args=parser.parse_args()

        #Call your Class
	ComputeTemporalAverage(args).Main()	
