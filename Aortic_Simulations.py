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
PRE_SOLVER_PATH = '/home/k/khanmu11/ana/Softwares/svSolver/BuildWithMake/Bin/svpre.exe'
SOLVER_PATH ='/home/k/khanmu11/ana/Softwares/svSolver/BuildWithMake/Bin/svsolver-openmpi.exe'
POST_SOLVER_PATH = '/home/k/khanmu11/ana/Softwares/svSolver/BuildWithMake/Bin/svpost.exe'

class AorticSimulations():
	def __init__(self)
		#Patient specific Assigned Parameters
		self.heart_rate=60 #beat per minute
		self.Pao_min=80 #minimum diastolic pressure (mmHg)
		self.Pao_max=120 #maximum systolic pressure (mmHg)
		self.StrokeVolume=60 #Stroke volume (mL)

		#Some constant parameters that don't change for each patient
		self.C_aorta=0.001 #Aortic compliance
		self.E_aorta=7000000 #elastic modulus of the aorta
		self.nuvw=0.5 #poisson ratio
		self.density=1.0 #tissue density
		self.kcons=0.83333 #shear constant
		self.pressure=pConv*(0.3333*self.Pao_max+0.666*self.Pao_min)
		self.InletWaveformName="inflow.flow" #Contains flow waveform
		self.FlowSplitFileName="FlowSplit.dat" #Contains resitances for each outlet cap

		#Get the number of caps
		self.CapNames=sorted(glob("./mesh-complete/mesh-surfaces/cap_*.vtp"))

		#Get the period from the flow file
		infile=open("inflow.flow")
		for Line in infile: 
			line=Line.split()
			if len(line)>0: self.Period=line[0]	
		print ("The Period is: %.05f"%self.Period)

	def Main(self):
		#Read the Flow Splits for each of the outlet
		FlowSplits=self.ReadFlowSplits()
		
		#Generate the 
		

	def WritePreSolveFile(self):
		svpre_file = open('AortaFSI.svpre','w')
		svpre_file.write('mesh_and_adjncy_vtu ' + mesh_complete + '\n')
 
		#Set Surface Id
		svpre_file.write('set_surface_id_vtp mesh-complete/mesh-complete.exterior.vtp 1\n')
		for i in range(len(self.CapNames)):
			svpre_file.write('set_surface_id_vtp %s %d\n'%(self.CapNames[i],i+2))
		
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
		
		#Write the deformable wall name
		svpre_file.write("deformable_wall_vtp mesh-complete/walls_combined.vtp\n")

		#Set the surface thickness
		Thickness=self.GetWallThickness()	
		for CapName_ in self.CapNames:	
			svpre_file.write("set_surface_thickness_vtp %s %.08f\n"%(CapName_,Thickness[CapName_]))
		#Write the WallThickness for the inlet
		InflowCap="mesh-complete/mesh-surfaces/inflow.vtp"
		InflowArea=self.GetSurfaceArea(InflowCap)
		InflowThickness=0.1*np.sqrt(InflowArea/np.pi)
		svpre_file.write("set_surface_thickness_vtp %s %.08f\n"%(InflowCap,InflowThickness))
		
		
			
			
		 
  
	

	def GetWallThickness(self):
		Thickness={}
		outlets=sorted(glob("./mesh-complete/mesh-surfaces/cap_*.vtp"))
		for outlet in outlets:
			SurfaceArea_=self.GetSurfaceArea(outlet)
			Thickness[outlet]=0.1*np.sqrt(SurfaceArea_/np.pi)
		return Thickness	
	
	#This function computes the surface area of the inlet or the outlets
	def GetSurfaceArea(Surface):
		masser = vtk.vtkMassProperties()
		masser.SetInputData(Surface)
		masser.Update()
		return masser.GetSurfaceArea()

	#This function will read user defined flow splits inside the file
	def ReadFlowSplits(self):
		infile=open(self.FlowSplitFileName,'r')
		#First colum contains outletname and the second outlet contains flow split
		FlowSplits={}
		for Line in FlowSplits:
			line=Line.split()
			FlowSplits[line[0]]=float(line[1])

		infile.close()

		
