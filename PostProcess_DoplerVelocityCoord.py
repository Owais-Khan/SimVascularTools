""" This script can be used for post processing of the 3D volumetric 
    results of the fine CFD simulations using an input surface of the
    Doppler Velocity location for coordinates and radius uses.

    Anahita A. Seresti
    February, 2024
"""

import os
import argparse

import vtk
import numpy as np

from utilities import ReadVTUFile, ReadVTPFile, GetCentroid, WriteVTUFile

class PostProcessingDoplerCoord():
    def __init__(self,Args):
        """Takes the input arguments of the parser

        :param Args: InputSurface containing the coordinates and InputFolder of the 3D results 
        :type Args: .vtp .vtu
        """
        self.Args = Args

    def GetCoord(self, Surface):
        """gets the input surface and extract the coordinates and radius

        :return: returns a sphere clipper centered at the coordinate and with a double sized radius
        :rtype: vtkSphereClip
        """
        centeroid = GetCentroid(Surface)
        '''Radius = ((Surface.GetRange[1]-Surface.GetRange[0])^2 + ...
                (Surface.GetRange[3]-Surface.GetRange[2])^2 + ...
                (Surface.GetRange[5]-Surface.GetRange[4])^2)^0.5/2
        '''
        #Defining the Sphere
        Sphere = vtk.vtkSphere()
        Sphere.SetCenter(centeroid)
        Sphere.SetRadius(0.35)#Radius*2)

        #Implement vtkclipping filter "sphere"
        clipper = vtk.vtkClipDataSet()
        clipper.SetClipFunction(Sphere)
        clipper.InsideOutOn()
        clipper.GetOutputInformation(1)
        clipper.Update()
        
        self.clipper = clipper
    
    def ComputeHmDy(self,volume):
        """Applys the clipper on the volume results and returns hemodynamics

        :param File: the 3D volumetric file at the output of CFD simulation
        :type File: VTU file
        :return: Hemodynamic properties
        :rtype: Velocity (cm/s), Pressure (mmHg)
        """
        # Applying the clipper on the 3D volumetric file
        self.clipper.SetInputData(volume)
        self.clipper.Update()
        SphereOutput = self.clipper.GetOutput()
        
        # Calculating the hemodynamic parameters
        Vmin, Vmax = SphereOutput.GetPointData().GetArray("velocity_mag_average").GetValueRange()
        #Pmin, Pmax = SphereOutput.GetPointData().GetArray("Pressure").GetValueRange()
        Velocity_ = np.empty(SphereOutput.GetNumberOfPoints())
        Pmean = 0
        for j in np.arange(0,SphereOutput.GetNumberOfPoints()):
            Velocity_[j] = SphereOutput.GetPointData().GetArray("velocity_mag_average").GetValue(j)
            Pmean += SphereOutput.GetPointData().GetArray("pressure").GetValue(j)
        
        Velocity_ = Velocity_[Velocity_ != 0] # remove the wall effect
        Vmean = np.mean(Velocity_)
        V95th = np.percentile(Velocity_,95)
        V50 = np.percentile(Velocity_,50)
        
        Pmean = Pmean/SphereOutput.GetNumberOfPoints()

        return Vmean, Vmin, Vmax, V95th, V50, Pmean
    
    def Main(self):
        """loops over input 3D results and applys clipper and saves the results in a textfile
        """
        # Read the input surface and get its coordinate to define a sphere clip
        Surface = ReadVTPFile(self.Args.InputSurface)
        print("--- Reading the input surface")
        self.GetCoord(Surface)
        print("--- Defining the sphere clipper")

        # Read 3D volumetric files within the input folder
        filenames = os.listdir(self.Args.InputFolder)
        filenames = [filename for filename in filenames if "vtu" in filename]
        filenames = sorted(filenames)
        N = len(filenames)
        ofile = f"{self.Args.InputFolder}/results.txt"
        
        with open(ofile,"w") as writefile:
            writefile.writelines("File Name, Vmean, Vmin, Vmax, 95th percentile V, Vmedian, Pmean \n")
            for n in np.arange(0,N):
                print("-"*25)
                print(f"--- Reading File: {filenames[n]}")
                volume = ReadVTUFile(f"{self.Args.InputFolder}/{filenames[n]}")
                print("--- Computing and Storing the hemodynamic features")
                writefile.writelines(f"{filenames[n]}: {' '.join(str(i) for i in self.ComputeHmDy(volume))} \n")
                WriteVTUFile(f"{self.Args.InputFolder}/clip_location_{filenames[n]}", self.clipper.GetOutput())


if __name__=="__main__":
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-InputSurface", "--InputSurface", type=str, dest="InputSurface", required=True, help="The location of Doppler Velocity Prob")
    parser.add_argument("-InputFolder", "--InputFolder", type=str, dest="InputFolder", required=True, help="Folder containing the 3D results of the CFD simulations")

    args = parser.parse_args()

    PostProcessingDoplerCoord(args).Main()