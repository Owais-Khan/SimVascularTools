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
                        self.Args.OutputFileName=self.Args.InputFolder+"../VolumeHemodynamics.vtp"

	def Main(self):
		


if __name__=="__main__":
        #Arguments
	parser= argparse.ArgumentParser(description="This script will compute hemodynamic indices such as average velocity, pressure, Turbulent kinetic engery, energy dissipation etc.")

        #Input filename for the mesh-complete folder
	parser.add_argument('-InputFolder', '--InputFolder', type=str, required=True, default="ResultsFolder", dest="InputFolder", help="The folder that contains the VTP and VTU files.")

        #Write the Output FileName
	parser.add_argument('-OutputFileName', '--OutputFileName', type=str, required=False, dest="OutputFileName", help="The output filename")

	parser.add_argument('-Period', '--Period', type=float, required=True, dest="Period", help="Period of the cardiac cycle in seconds.")

        #Put all the arguments together
	args=parser.parse_args()

        #Call your Class
	ComputeSurfaceHemodynamics(args).Main()
