# This script will convert SimVascular results to Oasis compatible Results
import sys
import os
from glob import glob
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import numpy as np

class ConvertSimVascualrToOasisResults():
	def __init__(self,Args):
		self.Args

	def Main(self,Args):
		#Load 	


if __name__=="__main__":
        #Arguments
	parser= argparse.ArgumentParser(description="This script will the Simvacsular results and convert them to oasis for postprocessing")

        #Input filename for the mesh-complete folder
	parser.add_argument('-InputFolder', '--InputFolder', type=str, required=True, dest="InputFolder", help="The path to the results folder containing the simvascular generate velocity files in .vtu format")

	parser.add_argument('-OutputFolder', '--OutputFolder', type=str, required=False, dest="OutputFolder", help="The folder path to store the convert SimVascular results to Oasis results.")

        #Put all the arguments together
	args=parser.parse_args()

        #Call your Class
	SimVascularToOasisResults(args).Main()
