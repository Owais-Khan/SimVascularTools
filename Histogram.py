import sys
import os
from glob import glob
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import numpy as np
import vtk
import argparse
from utilities import *
import matplotlib.pyplot as plt
import csv
import scipy.stats as st

class ImageAnalysisHistogramTecplot():
	def __init__(self,Args):
		self.Args=Args
		if self.Args.OutputFolder is None:
			self.Args.OutputFolder=self.Args.InputFileName.replace(self.Args.InputFileName.split("/")[-1],"/")
			
		self.Args.OutputFileName=self.Args.OutputFolder+"/"+self.Args.InputFileName.split("/")[-1][:-4]
		print ("Making a Directory to store results: %s"%self.Args.OutputFileName)
		os.system("mkdir %s"%self.Args.OutputFileName)

	def Main(self):
                #Read the vtu file
		print ("--- Reading %s"%self.Args.InputFileName)
		if self.Args.InputFileName[-4:]==".vtu":
			VelocityFile=ReadVTUFile(self.Args.InputFileName) #Read a VTU unstructured grid
		else:
			print ("The extension %s is not valid for volume"%self.Args.InputFileName[-4:])
			print ("Exiting...")
			exit(1)

		
		#Compute All of the Scalar quantities
		print ("--- Reading Scalar Values")
		DataArray=vtk_to_numpy(VelocityFile.GetPointData().GetArray(self.Args.ArrayName))
		DataArray=(DataArray[:,0]**2+DataArray[:,1]**2+DataArray[:,2]**2)**0.5
		DataArray = DataArray[(DataArray > 1e-3)]
		
		
		#Compute Mean and standard deviation
		print ("--- Computing Mean and Standard Deviation")
		mu, std = st.norm.fit (DataArray)
		print ("--- Computing 75th Percentile")
		Velocity75Perc=np.percentile(DataArray,75)
		print ("--- Computing the Median")
		Velocity50Perc=np.percentile(DataArray,50)
		print ("--- Performing Shapiro-Wilks test for Normality")
		ShapiroPvalue=st.shapiro(DataArray)[1]
		print ("--- Performing Kolmogorov Smirnov test for Normality")
		KolmogorovPvalue=st.kstest(DataArray,'norm')[1]
		print ("--- Computing Kurtosis and Skewness")
		Kurtosis=st.kurtosis(DataArray)
		Skewness=st.skew(DataArray)
		
		print ("--- Writing the Statistics to An OutputFile")
		outfile=open("%s_statistics.txt"%(self.Args.OutputFileName),'w')
		outfile.write("Mean Velocity:               %.05f\n"%mu)	
		outfile.write("Std Velocity:                %.05f\n"%std)	
		outfile.write("50Percentile:                %.05f\n"%Velocity50Perc)	
		outfile.write("75Percentile:                %.05f\n"%Velocity75Perc)	
		outfile.write("ShapriWilks p-value:         %.05f\n"%ShapiroPvalue)
		outfile.write("Kolmogorov Smirnov p-value:  %.05f\n"%KolmogorovPvalue)
		outfile.write("Kurtosis:                    %.05f\n"%Kurtosis)
		outfile.write("Skewness:                    %.05f\n"%Skewness)
		outfile.close()
	
		#Compute Histogram
		x_axis=np.linspace(DataArray.min(),DataArray.max(),self.Args.Bins)
		kde=st.gaussian_kde(DataArray)
		kde_pdf=kde.pdf(x_axis)

		#Now write the Tecplot Verson of the data
		print ("--- Writing Tecplot Histogram in %s_HistogramTecplot.dat"%(self.Args.OutputFileName))
		outfile=open("%s_Histogram_Tecplot.dat"%(self.Args.OutputFileName),'w')
		outfile.write('TITLE="Velocity Magnitude"\n')
		outfile.write('VARIABLES = "Velocity","PDF"\n')
		
		#First Write Normal Distribution of the Data
		outfile.write('Zone T= "NormalDistribution", I=%d, F=POINT\n'%self.Args.Bins)
		pdf_normal = st.norm.pdf(x_axis, mu, std)	
		for i in range(self.Args.Bins):
			outfile.write("%.05f %.05f\n"%(x_axis[i],pdf_normal[i]))

		#Now write the raw data
		outfile.write('Zone T= "DataDistribution", I=%d, F=POINT\n'%self.Args.Bins)
		for i in range(self.Args.Bins):
			outfile.write("%.05f %.05f\n"%(x_axis[i],kde_pdf[i]))
			
		#Now write the Mean and Standard Deviation Line
		outfile.write('Zone T= "Mean", I=2, F=POINT\n')
		outfile.write('%.05f %.05f\n'%(mu,0))
		outfile.write('%.05f %.05f\n'%(mu,1000))
		
		outfile.write('Zone T= "Std_negative", I=2, F=POINT\n')
		outfile.write('%.05f %.05f\n'%(mu-std,0))
		outfile.write('%.05f %.05f\n'%(mu-std,1000))
		
		outfile.write('Zone T= "Std_negative", I=2, F=POINT\n')
		outfile.write('%.05f %.05f\n'%(mu+std,0))
		outfile.write('%.05f %.05f\n'%(mu+std,1000))

if __name__=="__main__":
        #Description
	parser = argparse.ArgumentParser(description="This script will compute data statistics on velocity file produced by simvascular. A tecplot file will also be produced to visualize the distributions.")

	parser.add_argument('-InputFileName', '--InputFileName', type=str, required=True, dest="InputFileName",help="The vtu file that contains the Velocity Data.")
        
	parser.add_argument('-ArrayName', '--ArrayName', type=str, required=False, dest="ArrayName", default="velocity", help="The array name where the data is stored")
	
	parser.add_argument('-Bins', '--Bins', type=int, required=False, default=300, dest="Bins", help="The number of bins for the histogram.")

	#Output Filename 
	parser.add_argument('-OutputFolder', '--OutputFolder', type=str, required=False, dest="OutputFolder",help="The folder in which to write the results files.")	

	args=parser.parse_args()
	ImageAnalysisHistogramTecplot(args).Main()
