import os
import argparse
import sys
from glob import glob
from utilities import *
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import copy
class ExtractVelocityAwayFromWall():
	def __init__(self,Args):
		self.Args=Args
		os.system("mkdir %s"%self.Args.OutputFolder)
		outfolder_=self.Args.OutputFolder.split("/")[-1] 
		self.SkipOutFolder=self.Args.OutputFolder.replace(outfolder_,"/PINNsVelocity_Skip%d"%self.Args.Skip)
		os.system("mkdir %s"%self.SkipOutFolder)
	def Main(self):
		#Read the Files
		print ("--- Getting Velocity File Names...")
		InputFiles=sorted(glob(self.Args.InputFolder+"/*.vtu"))
		
		#Read the Surface file containing vmtkcenterlinemeshsections	
		SectionsFile=ReadVTPFile(self.Args.SurfaceFile)
		
		#Read the number of sections
		MaxSectionId=int(SectionsFile.GetPointData().GetArray("SectionIds").GetRange()[1])+1

		#Master list to store slice ids
		SliceIds=[]
		SliceIdsFlat=[]
		#Loop over all of the sections
		Coords=vtk_to_numpy(ReadVTUFile(InputFiles[0]).GetPoints().GetData())
		
		#Loop over all of the section ids and get corresponding velocity ids for the mesh
		print ("--- Getting the Closest Point of Slices in the Volumetric Mesh")	
		SliceIdsFlat=[]
		SliceIds=[[] for i in range(MaxSectionId)]
		for i in range(SectionsFile.GetNumberOfPoints()):
			point_,id_,dist_=ClosestPoint(SectionsFile.GetPoint(i),Coords)
			sliceId_=int(SectionsFile.GetPointData().GetArray("SectionIds").GetValue(i))
			SliceIdsFlat.append(id_)
			SliceIds[sliceId_].append(id_)
		
		#Velocity Array for each timestep
		velocity_=np.zeros(shape=(SectionsFile.GetNumberOfPoints(),3))
		
		for i in range(len(InputFiles)):
			print ("--- Looping over: %s"%InputFiles[i])
			#Read the VTU File
			infile_=ReadVTUFile(InputFiles[i])
			#Loop over all of the slice points
			for j in range(SectionsFile.GetNumberOfPoints()):
				id_=SliceIdsFlat[j]
				velocity_[j][0]=infile_.GetPointData().GetArray("Velocity").GetValue(id_*3)
				velocity_[j][1]=infile_.GetPointData().GetArray("Velocity").GetValue(id_*3+1)
				velocity_[j][2]=infile_.GetPointData().GetArray("Velocity").GetValue(id_*3+2)
				
			#Write the VTP File with Projected Velocities
			velocityVTK_=numpy_to_vtk(velocity_,deep=True)
			velocityVTK_.SetName("Velocity")
			SectionsFile.GetPointData().AddArray(velocityVTK_)
			SectionsFile.Modified()
			#Append to Surface
			WriteVTPFile(self.Args.OutputFolder+"/Velocity_%03d.vtp"%i,SectionsFile)	
		
			#Write the Files for PINNS Trainning
			outfile_=open(self.SkipOutFolder+"/Velocity_%03d.dat"%i,'w')
			outfile_.write("X Y Z U V W\n")
			for j in range(1,MaxSectionId-1,self.Args.Skip):
				ids_=set(SliceIds[j])
				for id_ in ids_:
					u_=infile_.GetPointData().GetArray("Velocity").GetValue(id_*3)
					v_=infile_.GetPointData().GetArray("Velocity").GetValue(id_*3+1)
					w_=infile_.GetPointData().GetArray("Velocity").GetValue(id_*3+2)
					coord_=infile_.GetPoint(id_)
					if np.linalg.norm([u_,v_,w_])!=0.0:
						outfile_.write("%.08f %.08f %.08f %.08f %.08f %.08f\n"%(coord_[0],coord_[1],coord_[2],u_,v_,w_))
			outfile_.close()
					
			
						




if __name__=="__main__":
        #Arguments
        
	parser= argparse.ArgumentParser(description="This script will extract velocity at given cross sections obtained from vmtkcenterlinemeshsections file. It will not interpolate but find the closest point velocity. The script will output surface section files + ascii file with X Y Z U V W for PINNs trainning")

       
	parser.add_argument('-InputFolder', '--InputFolder', type=str, required=True, dest="InputFolder", help="Folder containing the velocity data over a cardiac cycle.")
        
	parser.add_argument('-SurfaceFile', '--SurfaceFile', type=str, required=True, dest="SurfaceFile", help="Surface file containing vmtkcenterlinemeshsections ")
	
	parser.add_argument('-Skip', '--Skip', type=int,  required=True, default=1, dest="Skip", help="The number of slices to use from vmtkcenterlinemeshsections.")
	
	parser.add_argument('-OutputFolder', '--OutputFolder', type=str, required=True, dest="OutputFolder", help="The output folder to store the velocity and coordinate data.")

	
        #Put all the arguments together
	args=parser.parse_args()

    #Call your Class
	ExtractVelocityAwayFromWall(args).Main()

