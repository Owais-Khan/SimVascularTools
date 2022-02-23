import os
import argparse
import sys
from glob import glob
from utilities import *
import numpy as np
import vtk

class ComputeFlowDivision():
	def __init__(self,Args):
                self.Args=Args

	def Main(self):
		#Read all of the filenames from the surface file
		FileNames=sorted(glob(self.Args.ResultsFolder+"/all_results*.vtp"))
		
		#Skip best 	
		FileNames=FileNames[::self.Args.Skip]
		
		print ("Number of Files Detected: %d"%len(FileNames))
		
		#Loop over all of the files
		print ("Getting the Node Ids from all of the Caps")
		CapNodeIds={}
		CapArea={}
		CapNames=[self.Args.MeshFolder+"/mesh-surfaces/inflow.vtp"]
		CapNames+=sorted(glob(self.Args.MeshFolder+"/mesh-surfaces/cap*.vtp"))
		for CapName_ in CapNames:
			print ("--- Processing: %s"%CapName_.split("/")[-1])
			CapSurface_=ReadVTPFile(CapName_)
			CapNodeIds[CapName_]=[CapSurface_.GetPointData().GetArray("GlobalNodeID").GetValue(k) for k in range(CapSurface_.GetNumberOfPoints())]
			CapArea[CapName_]=GetSurfaceArea(CapSurface_)


		#Get the Cap Nodes on the the Mesh Surface
		print ("Reading Mesh Exterior Surface")
		MeshSurface=ReadVTPFile(FileNames[0])
		MeshSurfaceNodeIds=np.array([MeshSurface.GetPointData().GetArray("GlobalNodeID").GetValue(k) for k in range(MeshSurface.GetNumberOfPoints())])
		
		print ("Mapping Cap Nodes to Mesh Exterior Surface Nodes")
		CapToSurfaceIds={}
		for CapName_ in CapNames:
			print ("--- Looping over: %s"%CapName_)
			CapToSurfaceIds[CapName_]=np.zeros(len(CapNodeIds[CapName_]))
			counter=0
			for Value_ in CapNodeIds[CapName_]:
				CapToSurfaceIds[CapName_][counter]=np.where(MeshSurfaceNodeIds==Value_)[0]
				counter+=1
		
		#Loop over all of the Hemodynamics Files
		#Initialize the arrays
		print ("Initializing Arrays for Velocity and Pressures")
		CapVelocityMag={}
		CapPressure={}
		for CapName_ in CapNames: 
			CapVelocityMag[CapName_]=np.zeros(len(CapNodeIds[CapName_]))
			CapPressure[CapName_]=np.zeros(len(CapNodeIds[CapName_]))
		

		#Initialize an arry to store temporal data
		CapVelocityTime={}
		CapFlowRateTime={}
		CapPressureTime={}
		for CapName_ in CapNames:
			CapVelocityTime[CapName_]=[]
			CapFlowRateTime[CapName_]=[]
			CapPressureTime[CapName_]=[]

		#Loop over all of the files
		for FileName_ in FileNames:
			print ("---Looping over File: %s"%FileName_)
			SurfaceData_=ReadVTPFile(FileName_)
			for CapName_ in CapNames:
				VelMag_=[]
				Pressure_=[]
				print ("------ Looping over: %s"%CapName_)
				for i in range(0,len(CapVelocityMag[CapName_])):
					Id_=int(CapToSurfaceIds[CapName_][i])
					VelX_=SurfaceData_.GetPointData().GetArray("velocity").GetValue(3*Id_)
					VelY_=SurfaceData_.GetPointData().GetArray("velocity").GetValue(3*Id_+1)
					VelZ_=SurfaceData_.GetPointData().GetArray("velocity").GetValue(3*Id_+2)
					VelMag_.append((VelX_**2+VelY_**2+VelZ_**2)**0.5)
					Pressure_.append(SurfaceData_.GetPointData().GetArray("pressure").GetValue(Id_))
					CapVelocityMag[CapName_][i]+=(VelX_**2+VelY_**2+VelZ_**2)**0.5
					CapPressure[CapName_][i]+= SurfaceData_.GetPointData().GetArray("pressure").GetValue(Id_)
					
				CapVelocityTime[CapName_].append(np.average(VelMag_))
				CapFlowRateTime[CapName_].append(np.average(VelMag_)*CapArea[CapName_])
				CapPressureTime[CapName_].append(np.average(Pressure_)/1333.22)

		#Divide all Variables by the Number of Time Steps
		for CapName_ in CapNames:
			CapVelocityMag[CapName_]=CapVelocityMag[CapName_]/len(FileNames)
			CapPressure[CapName_]=CapPressure[CapName_]/len(FileNames)/1333.22

		#Now Compute the cross-sectional averaged velocity and flow rate
		CapFlowRate={}; CapVelocityAvg={}; CapPressureAvg={}
		for CapName_ in CapNames:
			CapVelocityAvg[CapName_]=np.average(CapVelocityMag[CapName_])	
			CapFlowRate[CapName_]   =np.average(CapVelocityMag[CapName_])*CapArea[CapName_]
			CapPressure[CapName_]=np.average(CapPressure[CapName_])	


		#Write the outflow waveforms
		print ("Writing Flow and Pressure Waveforms to FlowWaveforms Folder")
		os.system("mkdir "+self.Args.MeshFolder+"/../FlowWaveforms")
		Time=np.linspace(0,self.Args.Period,len(FileNames))
		for CapName_ in CapNames:
			outfile_=open(self.Args.MeshFolder+"/../FlowWaveforms/"+CapName_.split("/")[-1]+".dat",'w')
			outfile_.write("Time Velocity[cm/s] FlowRate[mL/s] Pressure[mmHg]\n")
			for i in range(0,len(FileNames)):
				outfile_.write("%.05f %.05f %.05f %.05f\n"%(Time[i],CapVelocityTime[CapName_][i],CapFlowRateTime[CapName_][i],CapPressureTime[CapName_][i]))
			outfile_.close()
	
		print ("\n\n")	
		for CapName_ in CapNames:
			print ("-"*20)
			print ("Flow Rate for %s is: %.02f mL/s"%(CapName_.split("/")[-1],CapFlowRate[CapName_]))	
			print ("Pressure for %s is: %.02f mmHg"%(CapName_.split("/")[-1],CapPressure[CapName_]))	
	

	


if __name__=="__main__":
        #Arguments
	parser= argparse.ArgumentParser(description="This script will extract the flow at each of the outlet.")

        #Input filename for the mesh-complete folder
	parser.add_argument('-MeshFolder', '--MeshFolder', type=str, required=True, default="./mesh-complete/", dest="MeshFolder", help="The mesh-complete folder.")
        
	#Skiping the files in the results folder
	parser.add_argument('-Skip', '--Skip', type=int, required=True, default=1, dest="Skip", help="The results file to skip.")
        
	#The folder that contains the vtp files
	parser.add_argument('-ResultsFolder', '--ResultsFolder', type=str, required=True, default="./vtp_files/", dest="ResultsFolder", help="The folder that contains the surface hemodynamics files.")
	
	parser.add_argument('-Period', '--Period', type=float, required=False, default=1, dest="Period", help="The period associated with the simulation.")

        #Put all the arguments together
	args=parser.parse_args()

        #Call your Class
	ComputeFlowDivision(args).Main()
