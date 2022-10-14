#This function will extract the probe point data
#from vtu files

class ExtractProbePoints():
	def __init__(self,args):
		self.Args=args







if __name__=="__main__":
        #Arguments
        parser= argparse.ArgumentParser(description="This script will the extract probe points from SimVascular vtu files")

        #Input filename for the mesh-complete folder
        parser.add_argument('-InputFolder', '--InputFolder', type=str, required=True, dest="InputFolder", help="The path to the results folder containing the simvascular generate velocity files in .vtu format")

        parser.add_argument('-OutputFolder', '--OutputFolder', type=str, required=False, dest="OutputFolder", help="The folder path to store the convert SimVascular results to Oasis results.")

        #Put all the arguments together
        args=parser.parse_args()

        #Call your Class
        SimVascularToOasisResults(args).Main()


