#This function will extract the probe point data
#from vtu files

class ExtractProbePoints():
	def __init__(self,args):
		self.Args=args







if __name__=="__main__":
        #Arguments
        parser= argparse.ArgumentParser(description="This script will the extract probe points from SimVascular vtu files")

        #Input filename for the mesh-complete folder
        parser.add_argument('-InputFileName', '--InputFileName', type=str, required=True, dest="InputFileName", help="This file contains the the x y z coordinates of the probe points.")
        
	
	parser.add_argument('-Radius', '--Radius', type=str, required=True, dest="Radius", help="The radius around the probe points to collect the data points.")


        parser.add_argument('-OutputFolder', '--OutputFolder', type=str, required=False, dest="OutputFolder", help="The folder path to store the probe points")

        #Put all the arguments together
        args=parser.parse_args()

        #Call your Class
        SimVascularToOasisResults(args).Main()


