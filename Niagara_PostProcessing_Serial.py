import sys
import os
import argparse

class Niagara_PostProcessing_Serial():
	def __init__(self,Args):
		self.Args=Args
		#Create an output folder if not defined by the user
		if self.Args.OutputFolder is None:
			os.system("mkdir results/")
			self.Args.OutputFolder="./results"
	
		#Batch Size
		self.BatchSize=(self.Args.End-self.Args.Start)/self.Args.Processors


	def Main(self):
		#Open a batch file that can process the results
		infile=open("Niagara_PostSolver",'w')

		#Write Niagara-specific parameters
		infile.write("#!/bin/bash\n")
		infile.write("# SLURM submission script for multiple serial jobs on Niagara\n")
		infile.write("#SBATCH --partition=compute\n")
		infile.write("#SBATCH --nodes=1\n")
		infile.write("#SBATCH --ntasks-per-node=40\n")
		infile.write("#SBATCH --time=10:00:00\n")
		infile.write("#SBATCH --job-name serialx40\n")
		infile.write("# Turn off implicit threading in Python, R\n")
		infile.write("export OMP_NUM_THREADS=1\n")

		#Write the post-solver parameters
		infile.write("# Name of the executable you want to run on the cluster\n")
		infile.write("module purge; module load cmake lsb-release intelpython3/2019u4 gcc/8.3.0 openmpi/4.0.1 vtk/9.0.1\n")

		start_=self.Args.Start
		end_=start_+self.BatchSize
		#Jump until you reach the end time step
		while end_ <= self.Args.End:
			infile.write('(/home/k/khanmu11/khanmu11/Softwares/svSolver/BuildWithMake/Bin/svpost.exe -indir %s -outdir %s -start %d -stop %d -incr %d -vtu "all_results.vtu" -vtp "all_results.vtp") &\n'%(self.Args.InputFolder, self.Args.OutputFolder, start_, end_, self.Args.Increment))

			#Update the start and end timesteps
			start_+=self.BatchSize
			end_+=self.BatchSize	
			
		infile.write("wait")
		infile.close()
	

	



if __name__=="__main__":
	 #Arguments
        
	parser= argparse.ArgumentParser(description="This script will generate a post-processing batch script for Niagara to generate velocity surface and volumetric files.")

	parser.add_argument('-InputFolder', '--InputFolder', type=str, required=True, dest="InputFolder", help="The folder that contains the VTP and VTU files.")
	parser.add_argument('-OutputFolder', '--OutputFolder', type=str, required=False, dest="OutputFolder", help="The folder to store the output files")
	
	parser.add_argument('-Start', '--Start', type=int, required=True, dest="Start", help="The start time step.")
	
	parser.add_argument('-End', '--End', type=int, required=True, dest="End", help="The end time step.")
	
	parser.add_argument('-Increment', '--Increment', type=int, required=True, dest="Increment", help="The end time step.")
	
	parser.add_argument('-Processors', '--Processors', type=int, default=20, required=False, dest="Processors", help="The number of processors to use in serial job.")



        #Put all the arguments together
	args=parser.parse_args()

        #Call your Class
	Niagara_PostProcessing_Serial (args).Main()
