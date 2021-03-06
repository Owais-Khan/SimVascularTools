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



---
---
---
## 2.0 SimVascular Scripts for Aortic Simulations
### 2.1 Tuning RCR Parameters on a Coarse mesh
We can utilize the pipeline below to ensure that the assigned flow rates match the flow rates predicted by the SimVascular CFD simulation. The script is hard-coded for the Niagara cluster on SciNet but can easily be modified for other clusters with minor modifications. To run the script, you will need 4 files
1. ```Aortic_Simulations.py```: The script that can be obtained from this Github repository. 
2. ```FlowSplit.dat```: File containing the desird flow splits for each outlet.
3. ```mesh-complete```: Folder generated by SimVascular that contains the volumetric mesh, and surfaces (refer to SimVascular documentations). **Note**: Please make sure this mesh is coarse (<350,000 elements).
4. ```inflow.flow```: File that contains the time vs. flow rate data (refer to SimVascular documentations).

Below is an example of what a typical ```FlowSplit.dat``` file will look like:
```console
cap_mesh-surfaceA.vtp 0.1
cap_mesh-surfaceB.vtp 0.5
cap_mesh-surfaceC.vtp 0.25
cap_mesh-surfaceD.vtp 0.15
```
---
Here are the specific instructions to run the tuning framework on Niagara cluster:
1. Create a new folder in the SciNet `$SCRATCH` directory and copy over the `mesh-complete` folder, `Aortic_Simulations.py`,`inflow.flow`,and `FlowSplit.dat` files from your computer directory. Make sure that all surface vtp file names in the `mesh-complete` folder match the `FlowSplit.dat` file
2. Load the necessary modules to run SimVascular using the following command line
```console
foo@bar:~$ module purge; module load cmake lsb-release intel/2019u4 intelmpi/2019u4 intelpython3/2019u4 gcc/8.3.0 openmpi/4.0.1 vtk/9.0.1; LoadSimVascularModules
``` 
3. Run the simulation 
```console
foo@bar:~$ python Aortic_Simulation.py
```
The log will be printed out on the screen and also written to the `Output.log` file.

--- 

### 2.2 Running the fine simulation using the tuned parameters.
Now that the tuning framework has finished, hopefully producing the flow rates that you prescribed in `FlowSplit.dat` file, we can use the files generated on the last iteration to run a fine simulation. You need to create a new folder to add the fine `mesh-complete` folder and copy the `Aortic_Simulations_fine.py` from this repository. You can run the following command to submit the simulation job:
```console
foo@bar:~$ python Aortic_Simulations_fine.py -InputFolder /path/to/the/Coarse/simulation/folder
```
This will copy all the necessary files from the Coarse simulation, and submit the job. Note that the temporal resolution will be twice that of the coarse simulation and 200 solution files will be stored for post-processing the results.

You can check the status of your job by typing the following command: `squeue -u [username]` 

You have the following options that you can modify for the script:
1. ```-Nodes```: The default number of nodes is 5. On Niagara, each node has 40 processors, so the simulation will run with 200 processors in total. You may choose to increase this if you have substantially large mesh (e.g., >5million elements).
2. ```-WCT```: This is the wall-clock-time parameter. The default is 24 hours. However, if you anticipate that your simulation will finish in less time, you can update it to a different number. 

---
### 2.3 Post-Processing the Results for visualization
Once your simulation has finished, you can run the postprocessing script to convert the results into file formats that you can load into Paraview. To do so, please type the following command inside your simulation working directory:

```console
foo@bar:~$ mkdir results
foo@bar:~$ degbugjob
foo@bar:~$ /home/k/khanmu11/khanmu11/Softwares/svSolver/BuildWithMake/Bin/svpost.exe -indir 200-procs_case/ -outdir results/ -start 12000 -stop 16000 -incr 20 -vtu "all_results.vtu" -vtp "all_results.vtp"
```
Note that `mkdir results` will create an emppty directory to store the post-processed results. `debugjob` will submit an interactive job on SciNet to run the postsolve. `svpost.exe` command will loop through 12000 to 16000 time step (i.e. 4th cardiac cycle) in increments of 80 time steps. This will generate 200 .vtu (volume) and 200.vtp (surface) files that we can use for computing various hemodynamic quantities. You can also download these files locally and visualize using Paraview.

---

