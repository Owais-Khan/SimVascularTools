import vtk
from vmtk import vtkvmtk, vmtkscripts
import numpy as np
import sys
import os
from glob import glob

class ComputeCenterlines():
	def __init__(self,Args):
		self.Args=Args
	def Main(self):
		print ("--- Processing Centerlines")
		#Read Cap Names
		print ("------ Getting Cap Names")
		CapNamesInlet,CapNamesOutlet=self.GetCapNames()
		#Get Cap Surfaces
		print ("------ Getting Cap Surfaces")
		CapSurfacesInlet=[self.ReadVtpFile(CapName) for CapName in CapNamesInlet]
		CapSurfacesOutlet=[self.ReadVtpFile(CapName) for CapName in CapNamesOutlet]
		#Compute Cap Center
		print ("------ Computing Cap Centers")
		CapPointsInlet=np.array([self.GetCapCenter(CapSurface) for CapSurface in CapSurfacesInlet])
		CapPointsOutlet=np.array([self.GetCapCenter(CapSurface) for CapSurface in CapSurfacesOutlet])
		#Compute Centerlines
		print ("------ Computing Centerlines\n")
		self.ComputeCenterlines(list(CapPointsInlet.flatten()),list(CapPointsOutlet.flatten()))

	def ComputeCenterlines(self,CapPointsInlet,CapPointsOutlet):
		print (CapPointsOutlet)
		exit(1)
		centerlines=vmtkscripts.vmtkCenterlines()
		centerlines.Surface= self.ReadVtpFile(self.Args.InputMeshFolder+"/walls_combined.vtp")
		centerlines.SeedSelectorName='openprofiles'#'pointlist'
		centerlines.AppendEndPoints=1
		#centerlines.Resampling=0
		#centerlines.ResamplingStepLength=0.05
		centerlines.SourcePoints=CapPointsInlet
		centerlines.TargetPoints=CapPointsOutlet
		centerlines.Execute()
		centerlines=centerlines.Centerlines

		"""SmoothenCenterlines=vmtkscripts.vmtkCenterlineSmoothing()
		SmoothenCenterlines.Centerlines=centerlines
		SmoothenCenterlines.iterations=100
		SmoothenCenterlines.factor=0.1
		SmoothenCenterlines.Execute()
		centerlines=SmoothenCenterlines.Centerlines"""

		#Save the compute centerlines
		print ("------ Storing Centerlines in %s/Centerlines.vtp"%self.Args.OutFolder)
		self.WriteVtpFile(centerlines,self.Args.OutFolder+"/Centerlines.vtp")

	def GetCapCenter(self,CapSurface):
		Center=vtk.vtkCenterOfMass()
		Center.SetInputData(CapSurface)
		Center.SetUseScalarsAsWeights(False)
		Center.Update()
		return Center.GetCenter()


	def ReadVtpFile(self,FileName):
		reader=vtk.vtkXMLPolyDataReader()
		reader.SetFileName(FileName)
		reader.Update()
		return reader.GetOutput()

	def WriteVtpFile(self,Surface,FileName):
		writer=vtk.vtkXMLPolyDataWriter()
		writer.SetFileName(FileName)
		writer.SetInputData(Surface)
		writer.Update()

	def GetCapNames(self):
		#Go inside mesh-complete/mesh-surfaces/ and get outlets
		CapNames=sorted(glob(self.Args.InputMeshFolder+"/mesh-surfaces/*.vtp"))
		CapNamesInlet=[]
		CapNamesOutlet=[]
		for CapName in CapNames:
			if CapName.find("wall_")<0 and CapName.find("inflow")<0:
				CapNamesOutlet.append(CapName)
				print (CapName)
			if CapName.find("inflow")>=0:
				CapNamesInlet.append(CapName)
		return CapNamesInlet,CapNamesOutlet
