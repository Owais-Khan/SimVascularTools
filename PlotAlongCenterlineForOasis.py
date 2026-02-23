import sys
import os
from glob import glob
from vtk.util.numpy_support import vtk_to_numpy as npvtk
import numpy as np
from scipy.fftpack import fftfreq,fft,ifft
import vtk
import argparse
from utilities import *

class PlotAlongCenterline():
	def __init__(self,Args):
		self.Args=Args
		

	def Main(self):

		if self.Args.OutputFolder is None:
			self.Args.OutputFolder="TemporalAvg"
			os.system("mkdir %s"%self.Args.OutputFolder)
		#Compute the centerlines
		CL_FileName=self.Args.OutputFolder+"/"+self.Args.InputSurface.replace(".vtp","_cl.vtp").split("/")[-1]
		InputFiles=sorted(glob(self.Args.InputFolder+"/*.vtu"))
	
		os.system("vmtkcenterlines -ifile %s -ofile %s -endpoints 1 -resampling 1 -resamplingstep 0.05"%(self.Args.InputSurface,CL_FileName))
		
		
		#Compute centerline sections
		OutputFileName=self.Args.OutputFolder+"/Results_Section"
		outfile = open(OutputFileName+".txt",'w')
		outfile.write('SectionID,Length,Surface_Area')

		AllData={}

		for file in InputFiles:
			File_=ReadVTUFile(file)

			os.system("vmtkcenterlinemeshsections -centerlinesfile %s -ifile %s -ofile %s.vtp"%(CL_FileName,file,OutputFileName))
	
			SurfaceSections=(ReadVTPFile(OutputFileName+".vtp"))

			#Get the number of section ids
			print(SurfaceSections)
			Nstart,Nend=SurfaceSections.GetPointData().GetArray("SectionIds").GetRange()
			Nstart=int(Nstart)
			Nend=int(Nend)
			Dist_=np.zeros(Nend-Nstart)
			Dist_sum=0
			SurfaceArea_=np.zeros(Nend-Nstart)

			for i in range(Nstart,Nend):
				#outfile.write('%d'%i)i

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

				#Compute the area
				SurfaceArea_[i]+=ComputeArea(sectionremeshed_)			
			
				#Get the length from the inlet
				Centroid_=GetCentroid(sectionremeshed_)

				if i==Nstart: Dist_sum=0
				else: 
					Dist_sum+= np.sqrt( (Centroid_[0]-Centroid_old_[0])**2 + (Centroid_[1]-Centroid_old_[1])**2 + (Centroid_[2]-Centroid_old_[2])**2 )
					Dist_[i]+=Dist_sum

				Centroid_old_=Centroid_
			
				#Loop over arrays
				for j in range (0,sectionremeshed_.GetPointData().GetNumberOfArrays()):
					#ArrayName
					ArrayName=sectionremeshed_.GetPointData().GetArrayName(j)					
					#Number of Points
					Npts=sectionremeshed_.GetNumberOfPoints()				

					#Write array names to text
					if i == Nstart:
						outfile.write(',%s_max,%s_avg'%(ArrayName,ArrayName))
						AllData[ArrayName+"_max"]=[]
						AllData[ArrayName+"_avg"]=[]

					#Save values to matrix
					min_value,max_value=sectionremeshed_.GetPointData().GetArray(ArrayName).GetValueRange()
					AllData[ArrayName+"_max"].append(max_value)
				
					#Get the Average Value
					avg_value=np.average(vtk_to_numpy(sectionremeshed_.GetPointData().GetArray(ArrayName)))	
					AllData[ArrayName+"_avg"].append(avg_value)

		outfile.close()
		outfile=open(OutputFileName+".txt",'r+')
		header=outfile.read()
		Array=header.split(',')[3:]
	
		outfile.write('\n')
	
		#Write values to text file
		for i in range(Nstart,Nend):
			outfile.write('%d %0.5f %0.5f '%(i,Dist_[i],SurfaceArea_[i]))

			for name in Array:
				outfile.write('%0.5f '%(AllData.get(name)[i]))

			outfile.write('\n')

		
if __name__=="__main__":
        #Description
	parser = argparse.ArgumentParser(description="This script will compute quantities along the centerline")
	parser.add_argument('-InputFolder', '--InputFolder', type=str, required=True, dest="InputFolder",help="The folder that contains the data in vtu format")
	
	parser.add_argument('-InputSurface', '--InputSurface', type=str, required=True, dest="InputSurface",help="The surface file along which to compute the data.")
        
	#Output Filename 
	parser.add_argument('-OutputFolder', '--OutputFolder', type=str, required=False, dest="OutputFolder",help="An output folder to store the centerlines and data")

	args=parser.parse_args()
	PlotAlongCenterline(args).Main()	
