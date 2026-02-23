import os
import argparse
import sys
from glob import glob
from utilities import *
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk

class ComputeSpatialAverage():
	def __init__(self,Args):
		self.Args=Args

	def Main(self):
		#Read the Files
		print ("---Reading the volumetric and surface files...")
		InfileVolume=ReadVTUFile(self.Args.InputFileName)
		self.Args.InputFileNameSurface=self.Args.InputFileName.replace(".vtu",".vtp").replace("Volumetric","Surface")
		InfileSurface=ReadVTPFile(self.Args.InputFileNameSurface)
		
		#Compute the Average
		print ("------ Looping over: %s"%self.Args.InputFileName)
		NumberOfArrays=InfileVolume.GetPointData().GetNumberOfArrays()
		for i in range(NumberOfArrays):
			ArrayName_=InfileVolume.GetPointData().GetArrayName(i)
			ArrayData_=vtk_to_numpy(InfileVolume.GetPointData().GetArray(ArrayName_))						
			Average_=np.average(ArrayData_)
			print ("Average for %s: %.08f"%(ArrayName_,Average_))	
				
                #Compute the Average
                
		print ("------ Looping over: %s"%self.Args.InputFileNameSurface)
		NumberOfArrays=InfileSurface.GetPointData().GetNumberOfArrays()
		for i in range(NumberOfArrays):
			ArrayName_=InfileSurface.GetPointData().GetArrayName(i)
			ArrayData_=vtk_to_numpy(InfileSurface.GetPointData().GetArray(ArrayName_))
			Average_=np.average(ArrayData_)
			print ("Average for %s: %.08f"%(ArrayName_,Average_))   


if __name__=="__main__":
        #Arguments
        parser= argparse.ArgumentParser(description="This script will compute spatial averages")

        parser.add_argument('-InputFileName', '--InputFileName', type=str, required=True, dest="InputFileName", help="This script will spatially average the results and report quantities.")   
                
        #Put all the arguments together
        args=parser.parse_args()
                
    #Call your Class
        ComputeSpatialAverage(args).Main()
