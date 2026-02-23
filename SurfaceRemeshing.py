import os
import sys
import argparse 

class SurfaceRemeshing():
	def __init__(self,Args):
		self.Args=Args
		self.SurfaceFileName=self.Args.InputFile.split("/")[-1]
		self.Args.OutputFolder=self.Args.InputFile.replace(self.Args.InputFile.split("/")[-1],"")
	def Main(self):
		#Compute the centerline
		os.system("vmtkcenterlines -ifile %s -ofile %s"%(self.Args.InputFile,self.Args.OutputFolder+self.SurfaceFileName.replace(".vtp","_centerlines.vtp")))
		



if __name__=="__main__":
        #Arguments
        parser= argparse.ArgumentParser(description="This script will generate a surface mesh based on diameter")

        #Input filename for the mesh-complete folder
        parser.add_argument('-InputFile', '--InputFile', type=str, required=True, dest="InputFile", help="The filename of the model surface")
        parser.add_argument('-OutputFolder', '--OutputFolder', type=str, required=False, dest="OutputFolder", help="The folder to store the output surface")

        #Input Folder where all of the results have been stored
        parser.add_argument('-ElementDensity', '--ElementDensity', type=int, required=False, dest="ElementDensity", help="The element density across a diameter.")

        #Put all the arguments together
        args=parser.parse_args()

        #Call your Class
        SurfaceRemeshing(args).Main()
