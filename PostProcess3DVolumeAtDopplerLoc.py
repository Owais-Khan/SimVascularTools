""" This script can be used for post processing of the 3D volumetric 
    results of the fine CFD simulations using an input surface of the
    Doppler Velocity location for coordinates and radius uses.

    Anahita A. Seresti Github: @aseresti
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
        self.OutputDirName = "PostProcessingOutputs"

    def GetCoord(self, SurfaceDir):
        """gets the input surface and extract the coordinates and radius

        :return: returns a sphere clipper centered at the coordinate and with a double sized radius
        :rtype: vtkSphereClip
        """
        Surface = ReadVTPFile(SurfaceDir)
        centeroid = GetCentroid(Surface)
        '''Radius = ((Surface.GetRange[1]-Surface.GetRange[0])^2 + ...
                (Surface.GetRange[3]-Surface.GetRange[2])^2 + ...
                (Surface.GetRange[5]-Surface.GetRange[4])^2)^0.5/2
        '''
        #Defining the Sphere
        Sphere = vtk.vtkSphere()
        Sphere.SetCenter(centeroid)
        Sphere.SetRadius(0.3)#Radius*2)

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
        Vmin, Vmax = SphereOutput.GetPointData().GetArray("Velocity").GetValueRange()
        #Pmin, Pmax = SphereOutput.GetPointData().GetArray("Pressure").GetValueRange()
        Velocity = np.empty(SphereOutput.GetNumberOfPoints())
        Pmean = 0
        for j in np.arange(0,SphereOutput.GetNumberOfPoints()):
            Velocity[j] = SphereOutput.GetPointData().GetArray("Velocity").GetValue(j)
            Pmean += SphereOutput.GetPointData().GetArray("Pressure").GetValue(j)

        Vmean = np.mean(Velocity)
        V95th = np.percentile(Velocity,95)
        V50 = np.percentile(Velocity,50)
        
        Pmean = Pmean/SphereOutput.GetNumberOfPoints()*0.00075006157584566 #Converting from dynes/square cm to mmHg

        return Vmean, Vmin, Vmax, V95th, V50, Pmean
    
    def Main(self):
        """loops over input 3D results and applys clipper and saves the results in a textfile
        """
        # Read the input surfaces
        SurfaceName = os.listdir(self.Args.InputSurfaceFolder)
        SurfaceName = [filename for filename in SurfaceName if "vtp" in filename]
        
        

        # Read 3D volumetric files within the input folder
        filenames = os.listdir(self.Args.InputFolder)

        #Create the Output Directory
        print("\n","="*10)
        if self.OutputDirName in filenames:
            print(f"--- Output Directory: {self.OutputDirName} Alredy exists!")
            print(f"--- Deleting the previous files within the Output Directory: {self.OutputDirName}")
            os.system(f"rm -rf {self.Args.InputFolder}/{self.OutputDirName}/*")
        else:
            os.system(f"mkdir {self.Args.InputFolder}/{self.OutputDirName}")
            print(f"--- Output Directory: {self.OutputDirName} Created!")

        #Extract the vtu files within the InputFolder
        filenames = [filename for filename in filenames if "vtu" in filename]
        filenames = sorted(filenames)

        for surface in SurfaceName:
                
                print("\n","="*10)
                location = surface.split(".")
                print(f"--- Reading the input surface: {surface}")
                
                
                # get surface coordinate to define a sphere clip
                self.GetCoord(f"{self.Args.InputSurfaceFolder}/{surface}")
                print("--- Defining the sphere clipper")
                
                ofile = f"{self.Args.InputFolder}/{self.OutputDirName}/results_{location[0]}.txt"
                
                with open(ofile,"w") as writefile:
                    writefile.writelines("File Name, Vmean, Vmin, Vmax, 95th percentile V, Vmedian, Pmean \n")
                    
                    for file in filenames:
                        print("-"*25)
                        print(f"--- Reading File: {file}")
                        volume = ReadVTUFile(f"{self.Args.InputFolder}/{file}")
                        
                        print("--- Computing and Storing the hemodynamic features")
                        writefile.writelines(f"{file}: {' '.join(str(i) for i in self.ComputeHmDy(volume))} \n")
                        WriteVTUFile(f"{self.Args.InputFolder}/{self.OutputDirName}/{location[0]}_{file}", self.clipper.GetOutput())


if __name__=="__main__":
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-InputSurfaceFolder", "--InputSurfaceFolder", type=str, dest="InputSurfaceFolder", required=True, help="Folder containing the surface files at the location of Doppler Velocity Prob defined as a single or several vtp files")
    parser.add_argument("-InputFolder", "--InputFolder", type=str, dest="InputFolder", required=True, help="Folder containing the 3D results of the CFD simulations")

    args = parser.parse_args()

    PostProcessingDoplerCoord(args).Main()