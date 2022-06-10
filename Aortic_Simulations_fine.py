import os
import argparse

class Aortic_Simulations_Fine():
	def __init__(self,Args):
		self.Args=Args

	def Main(self):
		#Copy all of the necessary files from the coarse folder
		os.system("cp %s/inflow.flow ./"%self.Args.InputFolder)	
		os.system("cp %s/Aorta.svpre ./"%self.Args.InputFolder)	
		os.system("cp %s/solver.inp ./"%self.Args.InputFolder)	
		os.system("cp %s/numstart.dat ./"%self.Args.InputFolder)	
		os.system("cp %s/rcrt.dat ./"%self.Args.InputFolder)
		self.WriteJobScript("Niagara_FlowSolver", 'svFlowsolver', self.Args.WCT, self.Args.Nodes,40, self.Args.Email)		
		
		os.system("/home/k/khanmu11/khanmu11/Softwares/svSolver/BuildWithMake/Bin/svpre.exe Aorta.svpre")

		#Load the solver.inp file and change number of timesteps and other parameters
		infile=open("solver.inp",'r')
		outfile=open("solver_temp.inp",'w')
		for LINE in infile:
			if LINE.find("Number of Timesteps:")>=0:
				line_=int(LINE.split()[-1])
				outfile.write("Number of Timesteps: %d\n"%(line_*2))
				Nsteps=line_*2
			elif LINE.find("Time Step Size:")>=0:
				line_=float(LINE.split()[-1])
				outfile.write("Time Step Size: %.08f\n"%(line_/2.))
			elif LINE.find("Number of Timesteps between Restarts")>=0:
				outfile.write("Number of Timesteps between Restarts: %d\n"%int(Nsteps/200.))
			else: outfile.write(LINE)
		outfile.close()
		infile.close()
		os.system("mv solver_temp.inp solver.inp")	
		os.system("sbatch Niagara_FlowSolver")	
	
        
	def WriteJobScript(self,scriptName,jobName,time,nodes,procs,email):
		scriptFile = open(str(scriptName), 'w')
		scriptFile.write('#!/bin/bash\n\n')
		scriptFile.write('# Name of your job\n')
		scriptFile.write('#SBATCH --job-name='+str(jobName)+'\n')
		scriptFile.write('#SBATCH --partition=compute\n\n')
		scriptFile.write('# Specify the name of the output file. The %j specifies the job ID\n')
		scriptFile.write('#SBATCH --output='+str(jobName)+'.o%j\n\n')
		scriptFile.write('# Specify the name of the error file. The %j specifies the job ID\n')
		scriptFile.write('#SBATCH --error='+str(jobName)+'.e%j\n\n')
		scriptFile.write('# The walltime you require for your job\n')
		scriptFile.write('#SBATCH --time='+str(time)+':00:00\n\n')
		scriptFile.write('# Job priority. Leave as normal for now\n')
		scriptFile.write('#SBATCH --qos=normal\n\n')
		scriptFile.write('# Number of nodes are you requesting for your job. You can have 40 processors per node\n')
		scriptFile.write('#SBATCH --nodes='+str(nodes)+'\n\n')
		scriptFile.write('# Number of processors per node\n')
		scriptFile.write('#SBATCH --ntasks-per-node='+str(procs)+'\n\n')
		scriptFile.write('# Send an email to this address when your job starts and finishes\n')
		scriptFile.write('#SBATCH --mail-user='+str(email)+'\n')
		scriptFile.write('#SBATCH --mail-type=begin\n')
		scriptFile.write('#SBATCH --mail-type=end\n\n')
		scriptFile.write('# Name of the executable you want to run on the cluster\n')
		scriptFile.write("module purge; module load cmake lsb-release intelpython3/2019u4 gcc/8.3.0 openmpi/4.0.1 vtk/9.0.1\n")

		scriptFile.write("srun /home/k/khanmu11/ana/Softwares/svSolver/BuildWithMake/Bin/svsolver-openmpi.exe\n")
		scriptFile.close()

if __name__=="__main__":
        #Arguments
	parser= argparse.ArgumentParser(description="This script will prepare and run the fine simulation based on the results obtained from the tuned coarsed results")

        #Input filename for the mesh-complete folder
	parser.add_argument('-InputFolder', '--InputFolder', type=str, required=True, dest="InputFolder", help="The folder that contains tuned files from the coarse simulations")
        
	parser.add_argument('-Nodes', '--Nodes', type=int, required=False, default=5, dest="Nodes", help="The number of compute nodes to use when running the fine simulation. Note that each node contains 40 processors on Niagara cluster")
	
	parser.add_argument('-WCT', '--WCT', type=int, required=False, default=24, dest="WCT", help="The wall clock time to use when running the fine simulation. The default is 24 hrs, which is the maximum allowed on the Niagara cluster.")
	
	parser.add_argument('-Email', '--Email', type=str, required=False, default="", dest="Email", help="The email address to use to send Niagara solver updates.")

        #Put all the arguments together
	args=parser.parse_args()

        #Call your Class
	Aortic_Simulations_Fine(args).Main()


	
