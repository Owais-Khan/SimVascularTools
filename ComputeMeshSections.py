import vtk
from vmtk import vtkvmtk, vmtkscripts
import numpy as np
import sys
import os
from glob import glob

class ComputeMeshSections():
	def __init__(self,Args):
		self.Args=Args

	def Main(self):
		#Read the Centerlines
		Centerlines=self.ReadVtpFile(self.Args.OutFolder+"/Centerlines.vtp")
		#Read the Mesh File
		Mesh=self.ReadVtuFile(self.Args.InputMeshFolder+"/mesh-complete.mesh.vtu")
		#ComputeMeshSections
		MeshSections=self.ComputeMeshSections(Centerlines,Mesh)

	def ComputeMeshSections(self,Centerlines,Mesh):
	


	def ReadVtuFile(self,FileName):
		reader=vtk.vtkXMLUnstructuredGridReader()
		reader.SetFileName(FileName)
		reader.Update()
		return reader.GetOutput()
        
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


