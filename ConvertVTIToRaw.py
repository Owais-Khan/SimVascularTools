if __name__=="__main__":
        #Description
        parser = argparse.ArgumentParser(description="This script will convert a 3D image file from a .vti (3d) format to .raw format.")

        #Provide a path to the Magnitude Images
        parser.add_argument('-InputFolder', '--InputFolder', type=str, required=True, dest="InputFolder",help="This folder contains the velocity files in vtu format.")

        parser.add_argument('-InputFile', '--InputFile', type=str, required=True, dest="InputFile",help="The input file contains a column list of x y z R (coordinates and Radius) for which data is needed.")

