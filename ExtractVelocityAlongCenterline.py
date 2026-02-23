import os
import argparse
import sys
from glob import glob
from utilities import *
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import copy

class ExtractVelocityAlongCenterline():
	def __init__(self,Args):
		self.Args=Args
		os.system("mkdir %s"%self.Args.OutputFolder)
		outfolder_=self.Args.OutputFolder.split("/")[-1] 
		self.OutputFolder=self.Args.InputFolder+"/../VelocityAlongCenterlineRadiusRatio%.02f"%self.Args.RadiusRatio
		self.OutputFolder2=self.Args.InputFolder+"/../VelocityAlongCenterlineRadiusRatio%.02f_Masked"%self.Args.RadiusRatio
		os.system("mkdir %s"%self.OutputFolder)
		os.system("mkdir %s"%self.OutputFolder2)

	def Main(self):
		#Read the Files
		print ("\n--- Getting Velocity File Names...")
		InputFiles=sorted(glob(self.Args.InputFolder+"/*.vtu"))
		
		#Read the Surface file containing vmtkcenterlinemeshsections	
		CenterlinesFile=ReadVTPFile(self.Args.CenterlinesFile)

		#Get the Coordinates from Centerlines File
		Npts=CenterlinesFile.GetNumberOfPoints()
		PointsCL=np.zeros(shape=(Npts,3))
		MISR=np.zeros(Npts)
		for i in range(Npts):
			point_=CenterlinesFile.GetPoints().GetPoint(i)
			misr_=CenterlinesFile.GetPointData().GetArray("MaximumInscribedSphereRadius").GetValue(i)
			PointsCL[i,0]=point_[0]	
			PointsCL[i,1]=point_[1]	
			PointsCL[i,2]=point_[2]	
			MISR[i]=misr_	
	
		#Loop over all of the sections
		Coords=vtk_to_numpy(ReadVTUFile(InputFiles[0]).GetPoints().GetData())
		
		#Loop over all of the section ids and get corresponding velocity ids for the mesh
		print ("\n--- Getting the Closest Point Along the Centerline with RadiusRatio")	
		PointIds=[]
		for i in range(Npts):
			dist_ = np.sum((Coords - PointsCL[i])**2, axis=1)
			dist_ = dist_**0.5
			for j in range(len(dist_)):
				if dist_[j]<(MISR[i]*self.Args.RadiusRatio): #only points within the MISR*radiusratio
					PointIds.append(j)
		PointIds=list(set(PointIds))
		print ("------ No. of Points Selected for Radius Ratio of %.02f is: %d"%(self.Args.RadiusRatio,len(PointIds)))	
			
		#Velocity Array for each timestep
		velocity_=np.zeros(shape=(len(PointIds),3))

		print ("\n")
		for i in range(len(InputFiles)):
			print ("--- Looping over: %s"%InputFiles[i])
			#Read the VTU File
			infile_=ReadVTUFile(InputFiles[i])
			
			#Store data in a text file
			outfile_=open(self.OutputFolder+"/Velocity_%03d.dat"%i,'w')
			outfile_.write("X Y Z U V W\n")
			
			#Loop over all of the slice points
			for j in range(len(PointIds)):
				id_=PointIds[j]
				u_=infile_.GetPointData().GetArray("Velocity").GetValue(id_*3)
				v_=infile_.GetPointData().GetArray("Velocity").GetValue(id_*3+1)
				w_=infile_.GetPointData().GetArray("Velocity").GetValue(id_*3+2)
				coord_=infile_.GetPoint(id_)
				if np.linalg.norm([u_,v_,w_])!=0.0:
						outfile_.write("%.08f %.08f %.08f %.08f %.08f %.08f\n"%(coord_[0],coord_[1],coord_[2],u_,v_,w_))


			for k in range(len(Coords)):
				if k not in PointIds:
					infile_.GetPointData().GetArray("Velocity").SetValue(k*3,  0)				
					infile_.GetPointData().GetArray("Velocity").SetValue(k*3+1,0)				
					infile_.GetPointData().GetArray("Velocity").SetValue(k*3+2,0)
			infile_.Modified()
			WriteVTUFile(self.OutputFolder2+"/"+InputFiles[i].split("/")[-1],infile_)	
			outfile_.close()
					
			
						




if __name__=="__main__":
        #Arguments
        
	parser= argparse.ArgumentParser(description="This script will extract velocity at given cross sections obtained from vmtkcenterlinemeshsections file. It will not interpolate but find the closest point velocity. The script will output surface section files + ascii file with X Y Z U V W for PINNs trainning")

	parser.add_argument('-InputFolder', '--InputFolder', type=str, required=True, dest="InputFolder", help="Folder containing the velocity data over a cardiac cycle.")
	parser.add_argument('-CenterlinesFile', '--CenterlinesFile', type=str, required=True, dest="CenterlinesFile", help="File with Centerlines")
	parser.add_argument('-RadiusRatio', '--RadiusRatio', type=float, required=True, dest="RadiusRatio", help="Ratio of the maximum inscribed sphere radius to collect data (between 0 and 1).")
	parser.add_argument('-OutputFolder', '--OutputFolder', type=str, required=True, dest="OutputFolder", help="The output folder to store the velocity and coordinate data.")

	
        #Put all the arguments together
	args=parser.parse_args()

    	#Call your Class
	ExtractVelocityAlongCenterline(args).Main()

