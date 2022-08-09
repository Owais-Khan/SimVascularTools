import sys
 import os
 from glob import glob
 from vtk.util.numpy_support import vtk_to_numpy as npvtk
 import numpy as np
 import vtk
 import argparse
 from utilities import *

 class PlotAlongCenterline():
 	def __init__(self,Args):
 		self.Args=Args


 	def Main(self):
 		#Compute the centerlines
 		CL_FileName=self.Args.OutputFolder+"/"+self.Args.InputSurface.replace(".vtp","_cl.vtp").split("/")[-1]

 		os.system("vmtkcenterlines -ifile %s -ofile %s -endpoints 1 -resampling 1 -resamplingstep 0.05"%(self.Args.InputSurface,CL_FileName))


 		#Compute centerline sections
 		OutputFileName=self.Args.OutputFolder+"/"+self.Args.InputVolume.split("/")[-1].replace(".vtu","_sections")
 		#print(OutputFileName)
 		os.system("vmtkcenterlinemeshsections -centerlinesfile %s -ifile %s -ofile %s.vtp"%(CL_FileName,self.Args.InputVolume,OutputFileName))



 		SurfaceSections=ReadVTPFile(OutputFileName+".vtp")


 		outfile = open(OutputFileName+".txt",'w')
 		outfile.write('SectionID,Length,FlowRate,Velocity_Max,Velocity_Average\n')

 		#Get the number of section ids
 		print(SurfaceSections)
 		Nstart,Nend=SurfaceSections.GetPointData().GetArray("SectionIds").GetRange()
 		Nstart=int(Nstart)
 		Nend=int(Nend)
 		Dist_=0
 		for i in range(Nstart,Nend):
 			#Get the Section of the Slice
 			section_ = ThresholdInBetween(SurfaceSections,"SectionIds",i,i)

 			#Extract a surface from the section_
 			surface_=vtk.vtkDataSetSurfaceFilter()
 			surface_.SetInputData(section_)
 			surface_.Update()

 			#Store the Slice 
 			WriteVTPFile(self.Args.OutputFolder+"/"+"temp.vtp",surface_.GetOutput())

 			#Triangulate the surface
 			os.system("vmtksurfacetriangle -ifile %s/temp.vtp -ofile %s/temp.vtp"%(self.Args.OutputFolder, self.Args.OutputFolder))

 			#Remesh the surface to a constant vensity
 			os.system("vmtksurfaceremeshing -ifile %s/temp.vtp -ofile %s/temp_remeshed.vtp -elementsizemode edgelength -edgelength 0.01"%(self.Args.OutputFolder,self.Args.OutputFolder))

 			#Project the velocity on the new mesh
 			os.system("vmtksurfaceprojection -ifile %s/temp_remeshed.vtp -rfile %s/temp.vtp -ofile %s/temp_proj.vtp"%(self.Args.OutputFolder,self.Args.OutputFolder,self.Args.OutputFolder))

 			#Read the Remeshed Surface File
 			sectionremeshed_=ReadVTPFile("%s/temp_proj.vtp"%self.Args.OutputFolder)

 			Vmin_,Vmax_ =sectionremeshed_.GetPointData().GetArray("velocity_mag_average").GetValueRange()

 			#Compute the area
 			SurfaceArea_=ComputeArea(sectionremeshed_)			
 			Vel_=0
 			for j in range(sectionremeshed_.GetNumberOfPoints()):
 				Vel_+=sectionremeshed_.GetPointData().GetArray("velocity_mag_average").GetValue(j)

 			Vel_=Vel_/sectionremeshed_.GetNumberOfPoints()

 			#Flow rate
 			Q_=Vel_*SurfaceArea_

 			#Get the length from the inlet
 			Centroid_=GetCentroid(sectionremeshed_)

 			if i==Nstart: Dist_=0
 			else: 
 				Dist_+= np.sqrt( (Centroid_[0]-Centroid_old_[0])**2 + (Centroid_[1]-Centroid_old_[1])**2 + (Centroid_[2]-Centroid_old_[2])**2 )

 			#print (Vel_,Vmax_,Q_,SurfaceArea_)
 			outfile.write('%d %.05f %.05f %.05f %.05f\n'%(i,Dist_,Q_,Vmax_,Vel_))


 			Centroid_old_=Centroid_


 			"""#Vmax.extend([Vmax_])
         		Vsection_ = []
         		Nstart_,Nend_ =section_.GetPointData().GetArray("velocity_mag_average").GetRange()
         		for j in range(int(Nstart_),int(Nend_)):
                 		Vsection_.extend([section_.GetPointData().GetArray("velocity_mag_average").GetValue(j)])
         		Vavg_ = sum(Vsection_)/len(Vsection_)
         		Vmax_ = max(Vsection_)
         		#Vavg.extend([sum(Vavg_)/len(Vavg_)])
         		with open(vel,'a') as writefile:
                 		writefile.write('%s, %s, %s\n'%(i,Vmax_,Vavg_))"""

 if __name__=="__main__":
         #Description
 	parser = argparse.ArgumentParser(description="This script will compute quantities along the centerline")
 	parser.add_argument('-InputVolume', '--InputVolume', type=str, required=True, dest="InputVolume",help="The vtu file that contains the data in vtu format")

 	parser.add_argument('-InputSurface', '--InputSurface', type=str, required=True, dest="InputSurface",help="The surface file along which to compute the data.")

 	#Output Filename 
 	parser.add_argument('-OutputFolder', '--OutputFolder', type=str, required=False, dest="OutputFolder",help="An output folder to store the centerlines and data")
 	args=parser.parse_args()
 	PlotAlongCenterline(args).Main()	
