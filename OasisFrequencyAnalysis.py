import sys
import os
from glob import glob
from pathlib import Path
from scipy.fftpack import fftfreq,fft,ifft
import numpy as np
import argparse
from dolfin import *
from compute_flow_and_simulation_metrics import get_dataset_names
parameters["reorder_dofs_serial"] = False

class OasisFrequencyAnalysis():
	def __init__(self,Args):
		self.Args=Args
		if self.Args.OutputFolder is None:
			os.system("mkdir %s/../Results_OasisFormat"%self.Args.InputFolder)
			self.Args.OutputFolder="%s/../Results_OasisFormat"%self.Args.InputFolder


	def Main(self):
		# Read the mesh information
		file_path_u = self.Args.InputFolder+ "/u.h5"
		if len(glob(self.Args.InputFolder+"/mesh.h5"))==1:
			mesh_path = glob(self.Args.InputFolder)+ "/mesh.h5"
			mesh = Mesh()
			with HDF5File(MPI.comm_world, mesh_path, "r") as mesh_file:
				mesh_file.read(mesh, "mesh", False)
    
		elif len(glob(self.Args.InputFolder+"/mesh.xml.gz"))==1:
			mesh_path = self.Args.InputFolder+"/mesh.xml.gz"
			mesh=Mesh(mesh_path)
    
		else:
			print ("No mesh found in .h5 or .xml.gz format")
			print ("Exiting...")
			exit(1)

		f = HDF5File(MPI.comm_world, file_path_u, "r")		

    
		# Get names of data to extract
		start = 0
		if MPI.rank(MPI.comm_world) == 0:
			print("The post processing starts from", start)
    		
		dataset_names = get_dataset_names(f, start=start)

		# Function space
		print ("Creating Function Spaces")
		V = VectorFunctionSpace(mesh, "CG", self.Args.velocity_degree)
		Vv = FunctionSpace(mesh, "CG", self.Args.velocity_degree)
    
		# Functions for storing values
		u = Function(V)
		Data_fs=Function(Vv)

		#Loop over all of the velocity files
		counter=0
		for data in dataset_names:
			print ("--- Looping over Timestep: %s"%data)
        
			if MPI.rank(MPI.comm_world) == 0:
				print(data)

        
			# Time step and velocity
			f.read(u, data)
			
			#Create an array to store all velocity at all timepoints
			if counter==0:
				Npts=int(len(u.vector()[:])/3.)
				Velocity_X=np.zeros(shape=(len(dataset_names),Npts))
				Velocity_Y=np.zeros(shape=(len(dataset_names),Npts))
				Velocity_Z=np.zeros(shape=(len(dataset_names),Npts))
				Velocity_Mag=np.zeros(shape=(len(dataset_names),Npts))

				N_ts=len(dataset_names)
				time=np.linspace(0,self.Args.Period,N_ts)
				self.W = fftfreq(N_ts, d=time[1]-time[0])


			#Store the velocity in the global arrays
			Velocity_X[counter,:]=u.vector()[0:Npts]
			Velocity_Y[counter,:]=u.vector()[Npts:2*Npts]
			Velocity_Z[counter,:]=u.vector()[2*Npts:]
			Velocity_Mag[counter,:]=np.sqrt(Velocity_X[counter,:]**2+Velocity_Y[counter,:]**2+Velocity_Z[counter,:]**2)
			counter += 1
		f.close()

		#Compute my frequency based analysis
		SPI=np.zeros(Npts)
		TKE=np.zeros(Npts)	
		for i in range(Npts):
			SPI[i]=self.filter_SPI(Velocity_Mag[:,i])
			U_tke_=np.mean(np.power(self.filter_TKE(Velocity_X[:,i]),2))		
			V_tke_=np.mean(np.power(self.filter_TKE(Velocity_Y[:,i]),2))		
			W_tke_=np.mean(np.power(self.filter_TKE(Velocity_Z[:,i]),2))		
			TKE[i]=0.5*(U_tke_+V_tke_+W_tke_)	

		#Store this in the results folder
		Data_fs.vector()[:]=SPI
		File("%s/SPI.pvd"%self.Args.OutputFolder)<<Data_fs
		Data_fs.vector()[:]=TKE
		File("%s/TKE.pvd"%self.Args.OutputFolder)<<Data_fs


        
	def filter_TKE(self,U):
                #For FKE
		U_fft       = fft(U)
		U_cut_fft   = U_fft.copy()
		U_cut_fft[(self.W<self.Args.CutoffFrequency)]=0
		U_ifft      =ifft(U_cut_fft)
		return U_ifft

        
	def filter_SPI(self,U):
                #for HI
		U_fft       = fft(U-np.mean(U))
                #Cut off the upper frequency to 25Hz
		U_fft_25Hz   = U_fft.copy()
		U_fft_25Hz[(self.W<self.Args.CutoffFrequency)]=0
                #Cut off the lower frequency to 0Hz
		U_fft_0Hz =U_fft.copy()
		U_fft_0Hz[(self.W<0)]=0
                #Compute the absolute value
		Power_25Hz    =np.sum ( np.power( np.absolute(U_fft_25Hz),2))
		Power_0Hz     =np.sum ( np.power( np.absolute(U_fft_0Hz) ,2))
		if Power_0Hz<1e-5: return 0
		else:              return Power_25Hz/Power_0Hz


if __name__=="__main__":
        #Description
	parser = argparse.ArgumentParser(description="This script will take results in Oasis format and compute frequency based biomarkers, such as FKE and SPI.")
	parser.add_argument('-InputFolder', '--InputFolder', type=str, required=True, dest="InputFolder",help="The path to folder containing the Oasis Results")
	parser.add_argument('-OutputFolder', '--OutputFolder', type=str, required=False, dest="OutputFolder",help="The folder to store the results.")
	parser.add_argument('-Period', '--Period', type=float, required=True, dest="Period",help="The duration of the cardiac cycle in seconds.")
	parser.add_argument('-CutoffFrequency', '--CutoffFrequency', type=int, required=False, default=25, dest="CutoffFrequency",help="The cut-off frequency to compute frequency-based biomarkers.")
	parser.add_argument('-velocity_degree', '--velocity_degree', type=int, required=False, default=1, dest="velocity_degree",help="The velocity degree for the function space")
	args=parser.parse_args()
	OasisFrequencyAnalysis(args).Main()


