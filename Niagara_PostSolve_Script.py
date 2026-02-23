import os
import argparse

class Niagara_PostSolve_Script():
	def __init__(self,Args):
		self.Args=Args

	def Main(self):
		infile=open(self.Args.InputFolder+"/../Niagara_PostSolve_Script.sh",'w')
		infile.write("#!/bin/bash\n")
		infile.write("#SBATCH --partition=compute\n")
		infile.write("# SLURM submission script for multiple serial jobs on Niagara\n")
		infile.write("#SBATCH --nodes=1\n")
		infile.write("#SBATCH --ntasks-per-node=40\n")
		infile.write("#SBATCH --time=5:00:00\n")
		infile.write("#SBATCH --job-name serialx40\n")
		infile.write("# Turn off implicit threading in Python, R\n")
		infile.write("export OMP_NUM_THREADS=1\n")
		infile.write("module purge; module load cmake lsb-release intelpython3/2019u4 gcc/8.3.0 openmpi/4.0.1 vtk/9.0.1\n")
		os.system("mkdir vtk_files")
		
		#Number of Files
		Nfiles=int((self.Args.Stop-self.Args.Start)/self.Args.Step)
		
		#The Start timestep
		TimeStep_=self.Args.Start
		Jump_=(self.Args.Stop-self.Args.Start)/10.
		procs=0
		while TimeStep_<self.Args.Stop:
			TimeStep1_=TimeStep_+Jump_
			command="(~/Softwares/svSolver/BuildWithMake/Bin/svpost.exe -start %d -stop %d -incr %d -indir %s -outdir %s -vtp all_results.vtp -vtu all_results.vtu > log_procs_%d) &\n"%(TimeStep_,TimeStep1_,self.Args.Step,self.Args.InputFolder,self.Args.InputFolder+"../vtk_files/",procs)
			infile.write(command)
			TimeStep_=TimeStep1_
			procs+=1
		infile.write("wait")

if __name__=="__main__":
        #Arguments
	parser= argparse.ArgumentParser(description="This script will compute the vtu and vtp files and store them in the 'vtk_files' folder")

        #Input filename for the mesh-complete folder
	parser.add_argument('-InputFolder', '--InputFolder', type=str, required=True, default="./", dest="InputFolder", help="The folder that contains the restart files.")
	
        #Input Folder where all of the results have been stored
	parser.add_argument('-Start', '--Start', type=int, required=True, dest="Start", help="The start time step.")
	parser.add_argument('-Stop', '--Stop', type=int, required=True, dest="Stop", help="The stop time step.")
	parser.add_argument('-Step', '--Step', type=int, required=True, dest="Step", help="The increment time step.")

        #Put all the arguments together
	args=parser.parse_args()
        
	#Call your Class
	Niagara_PostSolve_Script(args).Main()
		
