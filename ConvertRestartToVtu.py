import sys
import os
import argparse

class ConvertRestartToVtu():
	def __init__(self,Args):
		self.Args=Args

	def Main(self):
		os.system("/usr/local/sv/svsolver/2021-09-30/svpost -outdir '%s' -indir '%s' -start %d -stop %d -incr %d -vtu 'all_results.vtu' -vtp 'all_results.vtp'"%(self.Args.InputFolder,self.Args.InputFolder,self.Args.Start,self.Args.Stop, self.Args.Increment))


if __name__=="__main__":
	 #Arguments
	parser= argparse.ArgumentParser(description="This script will convert the restart files to vtu files")

        #Input Folder where all of the results have been stored
	parser.add_argument('-InputFolder', '--InputFolder', type=str, required=True, dest="InputFolder", help="The path to the folder that contains all of the restart files")

        #Start and End Argumenets                
	parser.add_argument('-Start', '--Start', type=int, required=True, dest="Start", help="The starting timestep to process the SimVascular files")

        #End Argument
	parser.add_argument('-Stop', '--Stop', type=int, required=True, dest="Stop", help="The End timestep to process the SimVascular files")
	
	#Increment
	parser.add_argument('-Increment', '--Increment', type=int, required=True, dest="Increment", help="The timestep increment")

	#Put all the arguments together
	args=parser.parse_args()

        #Call your Class
	ConvertRestartToVtu(args).Main()
		
