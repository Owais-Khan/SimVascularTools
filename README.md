# SimVascularTools
## SimVascular Scripts for Coronary Artery Simulations
### Automated Tuning of Lumped Parameter Network (LPN) for Coronary Simulation
We can increase the accuracy of a 3D coronary simulation by adding a lumped parameter network (LPN) to specify the boundary conditions. The LPN are described through a set of ordinary differential equations that use a circuit element representation of the heart and distal circulation. The equations consist of several parameters (e.g., resistance, compliance, inductance, etc) that need to be tuned to match clinical targets of the patients, such as blood pressure, heart rate, stroke volume and others. Manually tuning all of these parameters can be tedious.

We have developed an LPN parameter tuning framework that allows us to automatically tune these parameters to match the clinical targets of the patient. The pipeline requires runing a series of 3D simulations (rigid and FSI) until the parameters from the simulation match the clinical targets. The script ```CT_coronary_automated_tuning.py``` can perform the entire tuning pipeline, and has been tested on SciNet's Niagara cluster.

### Post-processing 3D results
Once the results folder is copied to a local computer, you can run the following scripts to process the data.
1. Convert the Restart Files to VTU and VTP files and start them in the InputFolder
```console
foo@bar:~$ python ConvertRestartToVtu.py -Start 3000 -Stop 4000 -Increment 10 -InputFolder 320-procs-folder
``` 
2. Extract Cut plans where flow data needs to be measured. This can be done in Paraview by slicing a plane normal to the vessel direction. These planes must be stored in a folder named ```ClippedPlanes``` inside the current working directory.
3. Copy over the ```mesh-complete``` folder in the current working directory.
4. Process the vtu and vtp files to compute average surface and volume files, and generate flow/pressure waveforms at all the outlets and user-defined clipped planes.
```console
foo@bar:~$ python PostProcess3DResults.py
```
The above script will generate a folder ```PostProcessedData``` which will contain the following files:
1. ```SurfaceData_TimeAveragedResults.vtp```: Contains all the surface quantities averaged over the cardiac cycle.
2. ```VolumeData_TimeAveragedResults.vtu```: Contains all the volumetric parameters averaged over the cardiac cycle.
3. ```CapData.dat```: File header contains time-averaged quantities (e.g., pressure, velocity) at each of the outlet. The remain file contains time-varying waveforms that can be loaded in TecPlot360.
4. ```NonCapData.dat```: Same as above but for user-defined clipped planes
5. ```CapSurfaceData```: Folder that contains cap surfaces at each of the outlet for all timesteps. Use for visualizing in boundary conditions in Paraview.

## SimVascular Scripts for Aortic Simulations
### Tuning RCR Parameters
1. Create a new folder in the SciNet ```$SCRATCH``` directory and copy over the ```mesh-complete``` folder, ```Aortic_Simulations.py```,```inflow.flow```,and ```FlowSplit.dat``` files from your computer directory. Make sure that all surface vtp file names in the ```mesh-complete``` folder match the ```FlowSplit.dat``` file
2. Load the necessary modules using the following command line
```console
foo@bar:~$ module purge; module load cmake lsb-release intelpython3/2019u4 gcc/8.3.0 openmpi/4.0.1 vtk/9.0.1
```
3. Update the flow splits inside the ```FlowSplit.dat``` file 
4. Run the simulation and adjust the period 
```console
foo@bar:~$ python Aortic_Simulation.py -Period 0.5
```
### Postprocessing Simulation Data 
