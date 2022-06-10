from time import sleep
import vtk
import os
import glob
import sys
import time
import numpy as np
import math
from scipy.optimize import minimize
pConv = 1333.34
USER_EMAIL_ADDRESS="owaiskhan@ryerson.ca"
PRE_SOLVER_PATH = '/home/k/khanmu11/ana/Softwares/svSolver/BuildWithMake/Bin/svpre.exe'
SOLVER_PATH ='/home/k/khanmu11/ana/Softwares/svSolver/BuildWithMake/Bin/svsolver-openmpi.exe'
POST_SOLVER_PATH = '/home/k/khanmu11/ana/Softwares/svSolver/BuildWithMake/Bin/svpost.exe'
RUN_COMMAND='srun'
FLOWSOLVER_NODES=2
WALL_CLOCK_TIME=2
PROCESSORS=40
CPUs=FLOWSOLVER_NODES*PROCESSORS
BATCH_COMMAND="sbatch "
class AorticSimulations():
	def __init__(self):
		#Patient specific Assigned Parameters
		self.Pao_min=24 #minimum diastolic pressure (mmHg)
		self.Pao_max=43.5 #maximum systolic pressure (mmHg)
		self.tag="Rigid"	
		#Solver paramters
		self.Ntimesteps=1000
		self.Ncycles=2
		self.Nrestart=50

		#Some constant parameters that don't change for each patient
		self.Pao_mean=0.333*self.Pao_max + 0.666*self.Pao_min
		self.E_Aorta=665890.2229 #elastic modulus of the aorta
		self.E_AorticBranches=self.E_Aorta*2.6 #Ratio from Table2 of Coogan et al. (2013 BMM)
		self.nuvw=0.5 #poisson ratio
		self.density=1.0 #tissue density
		self.kcons=0.83333 #shear constant
		self.pressure=pConv*(0.3333*self.Pao_max+0.666*self.Pao_min)
		self.InletWaveformName="inflow.flow" #Contains flow waveform
		self.FlowSplitFileName="FlowSplit.dat" #Contains resitances for each outlet cap

		#Compute total resistance
		self.C_Aorta=0.0005 #Aortic compliance
		self.C_AorticBranches=self.C_Aorta/10.
		
		#Get the number of caps
		self.CapNames=sorted(glob.glob("./mesh-complete/mesh-surfaces/cap_*.vtp"))

		#Get the period from the flow file
		infile=open("inflow.flow")
		self.FlowRate=0.
		counter=0.
		for Line in infile: 
			line=Line.split()
			if len(line)>0: 
				self.Period=float(line[0])
				self.FlowRate+=float(line[1])
				counter+=1.
		self.FlowRate/=(counter*-1)
		print ("The Flow Rate is: %.05f"%self.FlowRate)
		print ("The Period is: %.05f"%self.Period)


		#The total resistance
		self.Rtot=(self.Pao_mean*pConv)/self.FlowRate
		print ("The total resistance is: %.05f"%self.Rtot)

        
	def Main(self):
                #Read the Flow Splits for each of the outlet
		print ("Getting Flow Splits")
		FlowSplits=self.ReadFlowSplits()
		FlowSplits_assigned=self.ReadFlowSplits()
                
		#Generate the presolver file
		print ("Writing Presolver File")
		self.WritePreSolveFile()

                #Generate the Solver file
		print ("Writing Solver File")
		self.WriteSolverFile()
                
		#Create numstart file 
		infile=open("numstart.dat",'w')
		infile.write("0")
		infile.close()

                #Run presolver
		os.system(PRE_SOLVER_PATH+" Aorta.svpre")
                
		Error=100 #Error in flow splits
		Iteration=0.
                
		outfile=open("Output.log",'w')
		while Error>3:
                        #Write rcr file
			self.WriteRCRFile(FlowSplits_assigned)
                        
			#Run FSI Simulation
			if len(glob.glob("%s-procs_case"%CPUs))>0:
				os.system("mv %s-procs_case Try%d_%s-procs_case"%(CPUs,Iteration,CPUs))
                        
			self.Run3DSIM()

                        #Compute flow split and mean pressures for last cycle
			nStart=self.Ntimesteps*(self.Ncycles-1)
			nStop=self.Ntimesteps*(self.Ncycles)
			FlowRates=self.ComputeAverage("Q",nStart,nStop)
			Pressures=self.ComputeAverage("P",nStart,nStop)

                        #Computed Flow Splits
			SumQ=0
			FlowSplits_computed={}
			for CapName_ in self.CapNames: SumQ+=FlowRates[CapName_]
			for CapName_ in self.CapNames:
				FlowSplits_computed[CapName_]=FlowRates[CapName_]/SumQ

                        #Print the assigned and computed flow splits
			Error=0
			for CapName_ in self.CapNames:
				Q_split_com=100*FlowSplits_computed[CapName_]
				Q_split_desired=100*FlowSplits[CapName_]
				print ("\n")
				outfile.write("\n")
				print ("Q_desired for %s: %.03f"%(CapName_,Q_split_desired))
				print ("Q_computed for %s: %.03f"%(CapName_,Q_split_com))
                                
				outfile.write("Q_desired for %s: %.03f\n"%(CapName_,Q_split_desired))
				outfile.write("Q_computed for %s: %.03f\n"%(CapName_,Q_split_com))
				Error+=abs(Q_split_com-Q_split_desired)

			print ("-"*20)
			print ("The cumulative error in flow split is: %.03f"%Error)
			print ("Current Iterations is: %d"%(Iteration))
			outfile.write("The cumulative error in flow split is: %.03f\n"%Error)
			outfile.write("Current Iterations is: %d\n"%(Iteration))
                        
			#Update the flow splits
			for CapName_ in self.CapNames:
				FlowSplits_assigned[CapName_]*=(FlowSplits[CapName_]/FlowSplits_computed[CapName_])
				print ("Updated Flow Split for %s is: %.05f"%(CapName_,FlowSplits_assigned[CapName_]))
				outfile.write("Updated Flow Split for %s is: %.05f\n"%(CapName_,FlowSplits_assigned[CapName_]))
                        
			Iteration+=1
                        
			#Kill the loop after 10 iterations
			if Iteration==10:
				print ("Finished running 10 iterations")
				print ("Killing the infinite loop to save our planet")
				outfile.close()
				exit(1)
		outfile.close()




	def ComputeAverage(self,Parameter,nStart,nStop):
		if Parameter=="Q":
			infile=open(sorted(glob.glob("%d-procs_case/QHist*.dat.%d"%(CPUs,nStop)))[-1],'r')
		elif Parameter=="P":
			infile=open(sorted(glob.glob("%d-procs_case/PHist*.dat.%d"%(CPUs,nStop)))[-1],'r')
		else:
			print ("Neither Pressure nor Flow Rate was defined for 3D Post-processing")
			print ("Exiting...")
			exit(1)

		print(infile.readline())
		AverageParameter={}
		for CapName_ in self.CapNames: AverageParameter[CapName_]=0
		counter=nStart
		for LINE in infile:
			line=LINE.split()
			if counter>=nStart and counter<nStop: 
				for i in range(len(self.CapNames)):
					AverageParameter[self.CapNames[i]]+=float(line[i])
				counter+=1
	
		for CapName_ in self.CapNames: AverageParameter[CapName_]/=counter
		infile.close()
		return AverageParameter
	

	def Run3DSIM(self):
		self.WriteJobScript("Niagara_FlowSolver", 'svFlowsolver', WALL_CLOCK_TIME, FLOWSOLVER_NODES, PROCESSORS)
		flow_script = open("Niagara_FlowSolver", 'a')
		flow_script.write('\n')
		flow_script.write(RUN_COMMAND + ' '+SOLVER_PATH+"\n")
		flow_script.close()	
		
		#If there is no procs-folder, then run a 3D simulation
		procs_list = glob.glob('%d-procs_*'%CPUs)
		if(len(procs_list) == 0):
			print("\n** Starting 3D SimVascular Simulation\n")
			command_string = BATCH_COMMAND + " "+ "Niagara_FlowSolver"
			print(command_string)
			os.system(command_string)
		else:
			print ("\n** ERROR: There is already a *-procs_case folder present")
			print ("Delete it before a new simulation. Exiting...")
			exit(1)	

  
		# NEED FLAG HERE TO DETERMINE STOP OF 3D SIMULATION
		total_steps = self.Ntimesteps*self.Ncycles
		procs_list = glob.glob('%d-procs_case'%(int(FLOWSOLVER_NODES*PROCESSORS)))
		print('\n** Waiting for simulation job to launch...\n')

		#Check to see if the simulation has started by finding the results folder
		while(len(procs_list) == 0):
			time.sleep(180)
			procs_list = glob.glob('%d-procs_case'%(int(FLOWSOLVER_NODES*PROCESSORS)))
			print('\n** Waiting for simulation job to launch...\n')
		sim_folder = '%d-procs_case'%(int(FLOWSOLVER_NODES*PROCESSORS))

		print('\n** Simulation has launched! Waiting for simulation to finish...\n')

		#Check to see if the AllData file has been written inside the results folder
		AllData_file=glob.glob('%d-procs_case/AllData'%(int(FLOWSOLVER_NODES*PROCESSORS)))
		while(len(AllData_file) == 0):
			time.sleep(180)
			AllData_file = glob.glob('%d-procs_case/histor.dat'%(int(FLOWSOLVER_NODES*PROCESSORS)))
			print('\n** histor.dat has been written...\n')

  
		#Check to see if the Simulation is finished
		sim_finished = False
		while(not sim_finished):
			time.sleep(100)
			AllData_check = open(sim_folder + '/histor.dat', 'r')
			for line in AllData_check: counter=int(line.split()[0])

			print ("------- Finished %d of %d timesteps"%(counter,total_steps))
			AllData_check.close()
			if(counter == total_steps):
				sim_finished = True
		print('** Simulation has finished! Moving on to post processing...\n')

		
	
	def WriteJobScript(self,scriptName,jobName,time,nodes,procs):
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
		scriptFile.write('#SBATCH --mail-user='+str(USER_EMAIL_ADDRESS)+'\n')
		scriptFile.write('#SBATCH --mail-type=begin\n')
		scriptFile.write('#SBATCH --mail-type=end\n\n')
		scriptFile.write('# Name of the executable you want to run on the cluster\n')
		scriptFile.write("module purge; module load cmake lsb-release intelpython3/2019u4 gcc/8.3.0 openmpi/4.0.1 vtk/9.0.1\n")
		scriptFile.close()

	def WriteRCRFile(self,FlowSplits):
		rcr_file=open("rcrt.dat",'w')
		rcr_file.write("2\n")
		#Get the sum of flow split
		for CapName_ in self.CapNames: print (CapName_)
		QSum=sum([FlowSplits[CapName_] for CapName_ in self.CapNames])
		
		for CapName_ in self.CapNames:
			rcr_file.write("2\n")
			rcr_file.write("%.08f\n"%(QSum/FlowSplits[CapName_]*0.09*self.Rtot))
			if CapName_.find("aorta")>=0: rcr_file.write("%.08f\n"%self.C_Aorta)
			else: rcr_file.write("%.08f\n"%self.C_AorticBranches)
			rcr_file.write("%.08f\n"%(QSum/FlowSplits[CapName_]*0.91*self.Rtot))
			rcr_file.write("0.0 0\n")
			rcr_file.write("1.0 0\n")	
		rcr_file.close()
		
	def WriteSolverFile(self):
		solver_file=open("solver.inp",'w')
	
		solver_file.write("Density: 1.06\n")
		solver_file.write("Viscosity: 0.04\n")
		solver_file.write("\n")
		solver_file.write("Number of Timesteps: %d\n"%(self.Ncycles*self.Ntimesteps))
		solver_file.write("Time Step Size: %.08f\n"%(self.Period/self.Ntimesteps))
		solver_file.write("\n")
		solver_file.write("Number of Timesteps between Restarts: %d\n"%(self.Nrestart))
		solver_file.write("Number of Force Surfaces: 1\n")
		solver_file.write("Surface ID's for Force Calculation: 1\n")
		solver_file.write("Force Calculation Method: Velocity Based\n")
		solver_file.write("Print Average Solution: True\n")
		solver_file.write("Print Error Indicators: False\n")
		solver_file.write("\n")
		solver_file.write("Time Varying Boundary Conditions From File: True\n")
		solver_file.write("\n")
		solver_file.write("Step Construction: 0 1 0 1 0 1 0 1\n")
		solver_file.write("\n")
		solver_file.write("Number of RCR Surfaces: %d\n"%len(self.CapNames))
		solver_file.write("List of RCR Surfaces: ")
		for i in range(len(self.CapNames)): solver_file.write("%d "%(i+2))
		solver_file.write("\n")
		solver_file.write("RCR Values From File: True\n")
		solver_file.write("\n")
		if self.tag=="FSI":
			solver_file.write("Deformable Wall: True\n")
			solver_file.write("Variable Wall Thickness and Young Mod: True\n")
			solver_file.write("Density of Vessel Wall: 1\n")
			solver_file.write("Poisson Ratio of Vessel Wall: 0.5\n")
			solver_file.write("Shear Constant of Vessel Wall: 0.833333\n")
		
		solver_file.write("Pressure Coupling: Implicit\n")
		solver_file.write("Number of Coupled Surfaces: %d\n"%(len(self.CapNames)))
		solver_file.write("\n")
		solver_file.write("Backflow Stabilization Coefficient: 0.2\n")
		solver_file.write("Residual Control: True\n")
		solver_file.write("Residual Criteria: 1e-5\n")
		solver_file.write("Minimum Required Iterations: 10\n")
		solver_file.write("svLS Type: NS\n")
		solver_file.write("Number of Krylov Vectors per GMRES Sweep: 200\n")
		solver_file.write("Number of Solves per Left-hand-side Formation: 1\n")
		"""solver_file.write("Tolerance on Momentum Equations: 1e-5\n")
		solver_file.write("Tolerance on Continuity Equations: 1e-5\n")
		solver_file.write("Tolerance on svLS NS Solver: 1e-5\n")
		solver_file.write("Maximum Number of Iterations for svLS NS Solver: 5\n")
		solver_file.write("Maximum Number of Iterations for svLS Momentum Loop: 5\n")
		solver_file.write("Maximum Number of Iterations for svLS Continuity Loop: 1000\n")"""

	def WritePreSolveFile(self):
		svpre_file = open('Aorta.svpre','w')
		svpre_file.write('mesh_and_adjncy_vtu mesh-complete/mesh-complete.mesh.vtu\n')
                #Write the deformable wall name
		if self.tag=="FSI":
			svpre_file.write("deformable_wall_vtp mesh-complete/walls_combined.vtp\n")
			svpre_file.write("fix_free_edge_nodes_vtp mesh-complete/walls_combined.vtp\n")


		#Set Surface Id
		svpre_file.write('set_surface_id_vtp mesh-complete/mesh-complete.exterior.vtp 1\n')
		for i in range(len(self.CapNames)):
			svpre_file.write('set_surface_id_vtp %s %d\n'%(self.CapNames[i],i+2))
		svpre_file.write('set_surface_id_vtp ./mesh-complete/mesh-surfaces/inflow.vtp %d\n'%(i+3))
		
		#Write fluid paramters
		svpre_file.write("fluid_density 1.06\n")
		svpre_file.write("fluid_viscosity 0.04\n")
		svpre_file.write("initial_pressure 0\n")
		svpre_file.write("initial_velocity 0.0001 0.0001 0.0001\n")
		
		#Assign parabolic velocity to the inlet
		svpre_file.write("prescribed_velocities_vtp mesh-complete/mesh-surfaces/inflow.vtp\n")
		svpre_file.write("bct_analytical_shape parabolic\n")
		svpre_file.write("bct_period %.05f\n"%self.Period)
		svpre_file.write("bct_point_number 201\n")
		svpre_file.write("bct_fourier_mode_number 15\n")
		svpre_file.write("bct_create mesh-complete/mesh-surfaces/inflow.vtp inflow.flow\n")
		svpre_file.write("bct_write_dat bct.dat\n")
		svpre_file.write("bct_write_vtp bct.vtp\n")

		#Write zero pressure BC at the outlet
		for CapName_ in self.CapNames:
			svpre_file.write("pressure_vtp %s 0\n"%CapName_)
		
		if self.tag=="FSI":
			#Set the surface thickness
			Thickness=self.GetWallThickness()	
			for CapName_ in self.CapNames:	
				svpre_file.write("set_surface_thickness_vtp %s %.08f\n"%(CapName_,Thickness[CapName_]))

			#Write the WallThickness for the inlet
			InflowCap="mesh-complete/mesh-surfaces/inflow.vtp"
			InflowArea=self.GetSurfaceArea(InflowCap)
			InflowThickness=0.1*np.sqrt(InflowArea/np.pi)
			svpre_file.write("set_surface_thickness_vtp %s %.08f\n"%(InflowCap,InflowThickness))
		
			#Solve for Wall Thickness using Lapacian
			svpre_file.write("solve_varwall_thickness\n")		


			#Set the surface elastic modulus
			for CapName_ in self.CapNames:
				if CapName_.find("aorta")>=0:
					svpre_file.write("set_surface_E_vtp %s %.05f\n"%(CapName_,self.E_Aorta))
				else:
					svpre_file.write("set_surface_E_vtp %s %.05f\n"%(CapName_,self.E_AorticBranches))
			
			#Solve for Elastic Modulus
			svpre_file.write("solve_varwall_E\n")	
			
			
			#Write other stuff to do
			svpre_file.write("varwallprop_write_vtp varwallprop.vtp\n")
			svpre_file.write("deformable_nu 0.5\n")
			svpre_file.write("deformable_kcons 0.833333\n")
			svpre_file.write("deformable_pressure %.05f\n"%self.Pao_mean)
			svpre_file.write("deformable_solve_varwall_displacements\n")
			svpre_file.write("wall_displacements_write_vtp displacement.vtp\n")
		else:
			svpre_file.write("noslip_vtp mesh-complete/walls_combined.vtp\n")

		svpre_file.write("write_geombc geombc.dat.1\n")
		svpre_file.write("write_restart restart.0.1\n") 
  
	

	def GetWallThickness(self):
		Thickness={}
		outlets=sorted(glob.glob("./mesh-complete/mesh-surfaces/cap_*.vtp"))
		for outlet in outlets:
			SurfaceArea_=self.GetSurfaceArea(outlet)
			Thickness[outlet]=0.1*np.sqrt(SurfaceArea_/np.pi)
		return Thickness	
	
	#This function computes the surface area of the inlet or the outlets
	def GetSurfaceArea(self,Surface):
		masser = vtk.vtkMassProperties()
		masser.SetInputData(self.ReadVTPFile(Surface))
		masser.Update()
		return masser.GetSurfaceArea()

	def ReadVTPFile(self,FileName):
		reader=vtk.vtkXMLPolyDataReader()
		reader.SetFileName(FileName)
		reader.Update()
		return reader.GetOutput()

	#This function will read user defined flow splits inside the file
	def ReadFlowSplits(self):
		infile=open(self.FlowSplitFileName,'r')
		#First colum contains outletname and the second outlet contains flow split
		FlowSplits={}
		for Line in infile:
			line=Line.split()
			CaseName_=glob.glob("./mesh-complete/mesh-surfaces/*%s*"%line[0])[0]
			FlowSplits[CaseName_]=float(line[1])
		infile.close()
		return FlowSplits
if __name__=="__main__":
	AorticSimulations().Main()		
