import numpy as np
import os
from scipy.interpolate import splrep, splev
from scipy.fftpack import fftfreq,fft,ifft
import argparse
class ProcessFlowWavefrom:
	def __init__(self,Args):
		self.Args=Args
		if self.Args.OutputFile is None:
			self.Args.OutputFile=self.Args.InputFile.replace(".","_filtered.")
	def Main(self):
		#Read the flow waveofrm
		print ("Reading the flow waveform from %s"%self.Args.InputFile)
		WaveformFile=open(self.Args.InputFile,'r')
		X=[]; Y=[]
		for Line in WaveformFile:
			line=Line.split()
			X.append(float(line[0]))
			Y.append(float(line[1]))
		WaveformFile.close()

		#Convert the last Y values to first y value
		Y[-1]=Y[0]

		#Apply a fourier transform to the data
		print ("Applying a Spline Fitting to the flow waveform")	
		X_new,Y_new=self.SplineFitting(X,Y)
	
		#Create three cycles of X_new and Y_new
		dt=X_new[1]-X_new[0]
		X_new2=[]+list(X_new)
		Y_new2=[]+list(Y_new)
		for i in range(0,len(X_new)):
			X_new2.append(X_new2[-1]+dt)
			Y_new2.append(Y_new2[i])
		for i in range(0,len(X_new)):
			X_new2.append(X_new2[-1]+dt)
			Y_new2.append(Y_new2[i])

		#Now applying a fourier transform to remove higher modes
		print ("Applying a fourier filtered using %d modes"%self.Args.FourierModes)
		Y_filtered=self.FourierFiltering(X_new2,Y_new2)
		
		#Output the file
		print ("Writing the filtered waveform in: %s"%self.Args.OutputFile)
		outfile=open(self.Args.OutputFile,'w')
		for i in range(len(X_new)):
			outfile.write("%.05f %.05f\n"%(X_new[i]/1000.,-1*Y_filtered[i+len(X_new)]))
		outfile.close()

		"""plt.plot(X_new2,Y_new2)#,'-k',label="Original")
		#plt.plot(X_new2,Y_filtered,'-k',label="Filtered_FourierModes%d"%self.Args.FourierModes)
		plt.xlabel("Time")
		plt.ylabel("Flow Rate")
		plt.show()"""

	def FourierFiltering(self,X,Y):
		W = fftfreq(self.Args.NoOfPoints*3, d=X[1]-X[0])
                
		#Compute the FFR and IFR
		U_fft       = fft(Y)
		U_cut_fft   = U_fft.copy()
		U_cut_fft[(W>self.Args.FourierModes)]=0
		U_ifft      =ifft(U_cut_fft)
		return U_ifft.real	
		

	def SplineFitting(self,X,Y):
		X_new=np.linspace(X[0],X[-1],self.Args.NoOfPoints)
		spl=splrep(X,Y)
		Y_new=splev(X_new, spl)
		return X_new,Y_new

				




if __name__=="__main__":
        #Description
	parser = argparse.ArgumentParser(description="This script will apply four transform to the time vs velocity/flow rate series and generate a smoothed, periodic flow waveform.")

        #Provide a path to the Magnitude Images
	parser.add_argument('-InputFile', '--InputFile', type=str, required=True, dest="InputFile",help="The input file that contains the flow waveform data.")
        
	#The number of fourier modes
	parser.add_argument('-FourierModes', '--FourierModes', type=int, required=False, default=15, dest="FourierModes",help="The number of fourier modes to use.")
	
	#The number of output points
	parser.add_argument('-NoOfPoints', '--NoOfPoints', type=int, required=False, default=100, dest="NoOfPoints",help="The number of output points.")
	
	#The output file name
	parser.add_argument('-OutputFile', '--OutputFile', type=str, required=False, dest="OutputFile",help="The output file to store the processed flow waveform.")

	args=parser.parse_args()
	ProcessFlowWavefrom(args).Main()

