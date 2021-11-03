"""
This script will automatically initialize the "restart.0.1" and "geombc.dat.1"
needed to run an CMM-FSI simulation in SimVascular. Before running a CMM-FSI
simulation, the velocities and pressures at every node in the model need to be
initialized from a rigid simulation first under physiologic conditions.
This process is typically done manually, which commonly leads to user errors.
This script will automate this process to ease the burden on users who want
to perform CMM-FSI simulations as part of their analysis.

Justin S. Tran - December 2017
"""

import vtk
import os
import glob
import sys
import time
import numpy as np
import math

heart_rate = 49.0 #BPM
Pao_min = 57.0 # mmHg
Pao_max = 119.0 # mmHg
strokeVol = 99.0 # mL
ejectFract = 0.69
sys_cor_split=6.37

meanPressure = (0.333)*Pao_max + (0.667)*Pao_min # mmHg
Ppul_mean = 14.0 # mmHg
Qla_ratio = 'None'#1.07 # Ratio of 'early' to 'late' flows into LV (i.e. E/A wave)
mit_valve = 'None'#0.56 # Fraction of heart cycle that mitral valve is open for
aor_valve = 'None'#0.39 # Fraction of heart cycle that aortic valve is open for
pul_valve = 'None'#0.374 # Fraction of heart cycle that pulmonary valve is open for
Pra_mean = 'None'#3.0 # mmHg - IVC right atrial mean pressure
meanFlow = strokeVol*(float(heart_rate)/60.0) # mL/s
Cam_scale = 0.89
Ca_scale = 0.11

#Solver Parameters
Period =60/heart_rate
NumberOfCycles =4
NumberOfTimesteps=1000
SaveFrequency=10
ResidualCriteria=0.0001
# ******************* DEFINE USER INPUTS IN THIS BLOCK ************************

# Input the mean pressure in your patient/system in mmHg. Default value is
# 93.33 mmHg for a typical systemic blood pressure of 120/80. If you are
# simulating pulmonary anatomy, note that your mean pressures will typically
# be lower, in the range of 10 - 20 mmHg (unless you are simulating PH).
Pao_mean = meanPressure # mmHg

# Input the mean flow rate for your patient/system. For systems which include
# either the aortic valve or pulmonary valve, this will correspond to the
# cardiac output of the patient. For other anatomies, make sure to set this
# value correctly according to your clinical or literature data
Qmean = meanFlow # mL/s

# ** Mesh information
MESH_SURFACES_PATH = os.getcwd()+"/mesh-complete/mesh-surfaces"

# Remove the cap tags if present
CapFileNames=glob.glob(MESH_SURFACES_PATH+"/cap_*.vtp")
for CapFileName in CapFileNames:
        os.system("mv %s %s"%(CapFileName, CapFileName.replace("cap_","")))


INFLOW_TAG = 'inflow'

# ** Executables and file paths
SVSOLVER_PATH    = '/home/k/khanmu11/khanmu11/Softwares/svSolver/BuildWithMake/Bin/svsolver-openmpi.exe'
PRESOLVER_PATH = '/home/k/khanmu11/khanmu11/Softwares/svSolver/BuildWithMake/Bin/svpre.exe'
POSTSOLVER_PATH = '/home/k/khanmu11/khanmu11/Softwares/svSolver/BuildWithMake/Bin/svpost.exe'


# ** Cluster settings

FLOWSOLVER_NODES = 4
BATCH_COMMAND = 'sbatch'
RUN_COMMAND = 'srun' # usually defined inside your batch scripts
SOLVER_SCRIPT = 'Niagara_flowsolver_script'
USER_EMAIL_ADDRESS = 'owaiskhan@ryerson.edu'

# ************************ USER INPUTS END HERE *******************************

def readVTPFile( file_in ):
  reader = vtk.vtkXMLPolyDataReader()
  reader.SetFileName(file_in)
  reader.Update()
  poly = reader.GetOutput()
  return poly

#-------------------------------------------------------------------------------

def findVTPArea( file_in ):
  reader = vtk.vtkXMLPolyDataReader()
  reader.SetFileName(file_in)
  reader.Update()
  poly = reader.GetOutputPort()
  masser = vtk.vtkMassProperties()
  masser.SetInputConnection(poly)
  masser.Update()
  return masser.GetSurfaceArea()

#-------------------------------------------------------------------------------

def findVTPMassCenter( file_in ):
  reader = vtk.vtkXMLPolyDataReader()
  reader.SetFileName(file_in)
  reader.Update()
  poly = reader.GetOutput()
  num_points = poly.GetNumberOfPoints()
  x_list = []
  y_list = []
  z_list = []
  for i in range(num_points):
    x_list.append(float(poly.GetPoints().GetPoint(i)[0]))
    y_list.append(float(poly.GetPoints().GetPoint(i)[1]))
    z_list.append(float(poly.GetPoints().GetPoint(i)[2]))
   
  return (np.average(x_list), np.average(y_list), np.average(z_list))

#-------------------------------------------------------------------------------

def runScript_base(scriptName, jobName, time, nodes, procs):
  scriptFile = open(str(scriptName), 'w')
  scriptFile.write('#!/bin/bash\n\n')
  scriptFile.write('# Name of your job\n')
  scriptFile.write('#SBATCH --job-name='+str(jobName)+'\n')
  scriptFile.write('#SBATCH --partition=compute\n\n')
  scriptFile.write('# Specify the name of the output file. The %j specifies the job ID\n')
  scriptFile.write('#SBATCH --output='+str(jobName)+'.o%j\n\n')
  scriptFile.write('# Specify the name of the error file. The %j specifies the job ID\n')
  scriptFile.write('#SBATCH --error='+str(jobName)+'.e%j\n\n')
  scriptFile.write('# The walltime you require for your job\n')
  scriptFile.write('#SBATCH --time=2:00:00\n\n')
  scriptFile.write('# Job priority. Leave as normal for now\n')
  scriptFile.write('#SBATCH --qos=normal\n\n')
  scriptFile.write('# Number of nodes are you requesting for your job. You can have 24 processors per node\n')
  scriptFile.write('#SBATCH --nodes='+str(nodes)+'\n\n')
  scriptFile.write('# Number of processors per node\n')
  scriptFile.write('#SBATCH --ntasks-per-node='+str(procs)+'\n\n')
  scriptFile.write('# Send an email to this address when your job starts and finishes\n')
  scriptFile.write('#SBATCH --mail-user='+str(USER_EMAIL_ADDRESS)+'\n')
  scriptFile.write('#SBATCH --mail-type=begin\n')
  scriptFile.write('#SBATCH --mail-type=end\n\n')
  scriptFile.write('# Name of the executable you want to run on the cluster\n')
  scriptFile.close()
  
#-------------------------------------------------------------------------------

def makeFlowsolver_script():
  runScript_base(SOLVER_SCRIPT, 'svFlowsolver', 1, FLOWSOLVER_NODES, 40)
  flow_script = open(SOLVER_SCRIPT, 'a')
  flow_script.write('\n')
  flow_script.write(RUN_COMMAND + ' ' + SVSOLVER_PATH + '\n')
  flow_script.close()

#-------------------------------------------------------------------------------

def makePresolver_script():
  runScript_base('presolve_batch', 'svPre', 1, 1, 40)
  flow_script = open('presolve_batch', 'a')
  flow_script.write('\n')
  #flow_script.write(RUN_COMMAND + ' ' + PRESOLVER_PATH + ' coronaryLPN.svpre\n')
  flow_script.write(PRESOLVER_PATH + ' coronaryLPN.svpre\n')#MOK CHANGED
  flow_script.close()
  
#-------------------------------------------------------------------------------

def makePostsolver_script():
  runScript_base('postsolve_batch', 'svPost', 1, 1, 1)
  flow_script = open('postsolve_batch', 'a')
  flow_script.write('\n')
  flow_script.write(RUN_COMMAND + ' ' + POSTSOLVER_PATH + ' -sn 200 -ph -vtp all_results.vtp -vtu all_results.vtu\n')
  flow_script.close()
  
#-------------------------------------------------------------------------------

def rigidPresolve():

  # Get all the surfaces inside the mesh-surfaces folder. Append a * to ensure
  # we read in all the .vtp files
  input_path = MESH_SURFACES_PATH
  if(input_path[-1:] is not '/'):
    input_path = input_path + '/'
  if(input_path[-1:] is not '*'):
    input_path = input_path + '*'
    
  filelist_raw = glob.glob(input_path)
  filelist_raw.sort()

  filelist = []
  for trial in filelist_raw:
    if(trial[-4:] == ".vtp"):
      filelist.append(trial)
  
  # Compute the areas for all the outlets, and store their name
  # Note that we skip over the inlet face since we should NOT apply
  # a resistance here!    
  base_len = len(input_path)
  outlet_names = []
  outlet_areas = []
  inflow_name = ''
  
  for vtp_file in filelist:
    temp_name = vtp_file[base_len-1:]
    if(temp_name[:len(INFLOW_TAG)] == INFLOW_TAG):
      inflow_name = vtp_file
    elif(temp_name[:4] == 'wall'):
      continue
    else:
      outlet_names.append(vtp_file)
      outlet_areas.append(findVTPArea(vtp_file))
      
  # Using the specified mean pressure and mean flowrate, compute a total
  # resistance for the system
  Rtot = Pao_mean*1333.34/Qmean
  
  # Split this total resistance in parallel among the different outlets,
  # scaling the resistances inversely proportional to area (i.e. larger
  # outlets will have less resistance
 
  #Separate the total resistance to the left and right coronary tree
  #using a 70-30 flow split
  Rcor_right=Rtot/(0.3*sys_cor_split/100.)
  Rcor_left =Rtot/(0.7*sys_cor_split/100.)
  Raorta    =Rtot/(1-sys_cor_split/100.)  
 
  #Sum over all of the areas to the right coronary artery 
  A_sum_right = 0.0
  A_sum_left  = 0.0
  A_sum_aorta = 0.0
  for i in range(0, len(outlet_areas)):
    if outlet_names[i].find("rca")>=0:
      A_sum_right = A_sum_right + np.sqrt(float(outlet_areas[i]))**2.66
    elif outlet_names[i].find("lca")>=0:
      A_sum_left  = A_sum_left  + np.sqrt(float(outlet_areas[i]))**2.66
    elif outlet_names[i].find("aorta")>=0:
      A_sum_aorta = A_sum_aorta + float(outlet_areas[i])
    else:
       print ("Error: The outlet file name can't be assigned")
       exit(1) 
  
  #Split the resistances according to murray's law
  outlet_resistances = np.zeros(len(outlet_areas))

  for i in range(0, len(outlet_areas)):
    if outlet_names[i].find("rca")>=0:
      outlet_resistances[i] = Rcor_right*(A_sum_right/np.sqrt(outlet_areas[i])**2.66)
    elif outlet_names[i].find("lca")>=0:
      outlet_resistances[i] = Rcor_left*(A_sum_left/np.sqrt(outlet_areas[i])**2.66)
    elif outlet_names[i].find("aorta")>=0:
      outlet_resistances[i] = Raorta*(A_sum_aorta/outlet_areas[i])

  nFaces = len(outlet_areas)
  
  # Make the solver.inp
  solver_writer = open('solver.inp','w')

  solver_writer.write('# NOT FOR GENERAL USE, ONLY FOR INITIALIZING CMM-FSI\n\n')
  solver_writer.write('# Phasta Version 1.5 Input File\n')
  solver_writer.write('# Produced by initialize_CMM.py\n')

  solver_writer.write('# ****** SOLUTION CONTROL ****** #\n')
  solver_writer.write('Equation of State: Incompressible\n')
  solver_writer.write('Number of Timesteps: 200\n')
  solver_writer.write('Time Step Size: {0:.6f}\n\n'.format(0.004))

  solver_writer.write('# ***** FLUID PROPERTIES ***** #\n')
  solver_writer.write('Viscosity: 0.04\n')
  solver_writer.write('Density: 1.06\n\n')

  solver_writer.write('# ***** OUTPUT CONTROL ***** #\n')
  solver_writer.write('Number of Timesteps between Restarts: 50\n')
  solver_writer.write('Print ybar: True\n\n')

  solver_writer.write('# ****** CARDIOVASCULAR MODELING PARAMETERS ***** #\n')
  solver_writer.write('Time Varying Boundary Conditions From File: True\n')
  solver_writer.write('Number of Coupled Surfaces: {0}\n'.format(nFaces))
  solver_writer.write('Pressure Coupling: Implicit\n')
  solver_writer.write('Number of Resistance Surfaces: {0}\n'.format(nFaces))
  neumann_list = 'List of Resistance Surfaces: '
  for i in range(2, nFaces+2, 1):
    neumann_list = neumann_list + str(i) + ' '
  solver_writer.write(neumann_list + '\n')
  resistance_values = 'Resistance Values: '
  for i in range(0, nFaces, 1):
    resistance_values = resistance_values + str(outlet_resistances[i]) + ' '
  solver_writer.write(resistance_values + '\n')
  solver_writer.write('Find the GenBC Inside the Running Directory: False\n')
  solver_writer.write('Number of Timesteps for GenBC Initialization: 0\n')

  solver_writer.write('# ***** LINEAR SOLVER (SHOULD LEAVE THIS ALONE) ****** #\n')
  solver_writer.write('Number of Solves per Left-hand-side Formation: 1\n\n')
  solver_writer.write('Minimum Required Iterations : 3\n')
  solver_writer.write('Maximum Number of Iterations for svLS NS Solver : 10\n')
  solver_writer.write('Maximum Number of Iterations for svLS Momentum Loop: 2\n')
  solver_writer.write('Maximum Number of Iterations for svLS Continuity Loop: 400\n\n')

  solver_writer.write('# ***** DISCRETIZATION CONTROL ***** #\n')
  solver_writer.write('Time Integration Rule: Second Order\n')
  solver_writer.write('Time Integration Rho Infinity: 0.2\n')
  solver_writer.write('Flow Advection Form: Convective\n')
  solver_writer.write('Number of Force Surfaces: 1\n')
  solver_writer.write("Surface ID's for Force Calculation: 1\n\n")

  solver_writer.write('Residual Control: True\n')
  solver_writer.write('Residual Criteria: 0.001\n')
  solver_writer.write('Step Construction : 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1\n\n')

  solver_writer.write('# NOTE: Recommended to change this to GMRES for Deformable simulations #\n')
  solver_writer.write('svLS Type: NS\n')
  solver_writer.write('Backflow Stabilization Coefficient: 0.25\n\n')

  solver_writer.close()
  
  # Make a .flow file for the inlet flow
  flow_file = open('inflow.flow', 'w')
  
  flow_file.write('0.0 ' + str(Qmean) + '\n')
  flow_file.write('1.0 ' + str(Qmean) + '\n')
  
  flow_file.close()
  
  # Generate the .svpre script for the pre-solver
  svpre_file = open('coronaryLPN.svpre','w')
   
  mesh_complete = input_path[:-1] + '../' + 'mesh-complete.mesh.vtu'
  walls_combined = input_path[:-1] + '../' + 'walls_combined.vtp'
  mesh_exterior = input_path[:-1] + '../' + 'mesh-complete.exterior.vtp'
  zp_string = 'zero_pressure_vtp '
  pv_string = 'prescribed_velocities_vtp '
  id_string = 'set_surface_id_vtp '

  svpre_file.write('mesh_and_adjncy_vtu ' + mesh_complete + '\n')
  svpre_file.write('noslip_vtp ' + walls_combined + '\n\n')

  for name in outlet_names:
    temp_string = zp_string + name + '\n'
    svpre_file.write(temp_string)
  temp_string = pv_string + inflow_name + '\n\n'
  svpre_file.write(temp_string)
  
  temp_string = id_string + mesh_exterior + ' 1\n'
  svpre_file.write(temp_string)
  count = 2

  for name in outlet_names:
    temp_string = id_string + name + ' ' + str(count) + '\n'
    svpre_file.write(temp_string)
    count = count + 1
  
  svpre_file.write(id_string + inflow_name + ' ' + str(count) + '\n\n')
  
  svpre_file.write('bct_analytical_shape parabolic\n')
  svpre_file.write('bct_period 1.0\n')
  svpre_file.write('bct_point_number 201\n')
  svpre_file.write('bct_fourier_mode_number 10\n')
  svpre_file.write('bct_create ' + inflow_name + ' inflow.flow\n')
  svpre_file.write('bct_write_dat bct.dat\n')
  svpre_file.write('bct_write_vtp bct.vtp\n')

  svpre_file.write('write_geombc geombc.dat.1\n')
  svpre_file.write('initial_pressure 119990\n')
  svpre_file.write('write_restart restart.0.1')

  svpre_file.close()
  
  # Generate numstart.dat
  numstart_writer = open('numstart.dat', 'w')
  numstart_writer.write('0')
  numstart_writer.close()
  
  # Run the pre-solver
  makePresolver_script()
  command_string = BATCH_COMMAND + ' presolve_batch'
  print(command_string)
  os.system(command_string)
  
  presolve_start = False
  presolve_job_ID = '-1'
  while(not presolve_start):
  
    time.sleep(30)
    presolve_check = glob.glob('svPre.o*')
    if(len(presolve_check) > 0):
      presolve_start = True
      presolve_job_ID = presolve_check[0][7:]
      
  print('Pre-solve started! Waiting for it to finish...')
  print('jobID = ' + presolve_job_ID)
  presolve_finished = False
  while(not presolve_finished):
    
    time.sleep(30)
    command_string = 'squeue -j ' + presolve_job_ID + ' > queue_status.txt'
    os.system(command_string)
    
    status_check = open('queue_status.txt', 'r')
    line_check = 0
    for line in status_check:
      line_check = line_check+1
    
    status_check.close()
    command_string = 'rm queue_status.txt'
    os.system(command_string)
    
    if(line_check < 1):
      presolve_finished = True
  
  # Check to make sure all the solver files are here
  bctThere = False
  restartThere = False
  geombcThere = False
  solverThere = False
  numstartThere = False
  
  bct_check = glob.glob('bct.dat')
  if(len(bct_check) > 0):
    bctThere = True
  restart_check = glob.glob('restart.0.1')
  if(len(restart_check) > 0):
    restartThere = True
  geombc_check = glob.glob('geombc.dat.1')
  if(len(geombc_check) > 0):
    geombcThere = True
  solver_check = glob.glob('solver.inp')
  if(len(solver_check) > 0):
    solverThere = True
  numstart_check = glob.glob('numstart.dat')
  if(len(numstart_check) > 0):
    numstartThere = True
  assert(bctThere)
  assert(restartThere)
  assert(geombcThere)
  assert(solverThere)
  assert(numstartThere)
  
  # Move all the solver files into a new folder
  command_string = 'mkdir -p solver_files'
  print(command_string)
  os.system(command_string)
  
  command_string = 'mv bct.dat restart.0.1 geombc.dat.1 solver.inp numstart.dat solver_files/'
  print(command_string)
  os.system(command_string)
  
  command_string = 'cd solver_files/'
  print(command_string)
  os.chdir('solver_files/')
  
  # Run svSolver
  makeFlowsolver_script()
  command_string = BATCH_COMMAND + ' ' + SOLVER_SCRIPT
  print(command_string)
  os.system(command_string)
  
  solver_start = False
  solver_job_ID = '-1'
  while(not solver_start):
  
    time.sleep(30)
    solver_check = glob.glob('svFlowsolver.o*')
    if(len(solver_check) > 0):
      solver_start = True
      solver_job_ID = solver_check[0][14:]
      
  print('Solve started! Waiting for it to finish...')
  print('jobID = ' + solver_job_ID)
  solver_finished = False
  while(not solver_finished):
    
    time.sleep(60)
    command_string = 'squeue -j ' + solver_job_ID + ' > queue_status.txt'
    os.system(command_string)
    
    status_check = open('queue_status.txt', 'r')
    line_check = 0
    for line in status_check:
      line_check = line_check+1
    
    status_check.close()
    command_string = 'rm queue_status.txt'
    os.system(command_string)
    
    if(line_check < 1):
      solver_finished = True
  
  # Run post-solver inside the n-procs_case
  proc_folder = glob.glob('*procs_case')
  proc_folder = proc_folder[0]
  
  print('svSolver has finished! Running post-solver...')
  command_string = 'cd ' + proc_folder
  print(command_string)
  os.chdir(proc_folder)
  
  makePostsolver_script()
  command_string = BATCH_COMMAND + ' postsolve_batch'
  print(command_string)
  os.system(command_string)
  
  postsolve_start = False
  postsolve_job_ID = '-1'
  while(not postsolve_start):
  
    time.sleep(30)
    postsolve_check = glob.glob('svPost.o*')
    if(len(postsolve_check) > 0):
      postsolve_start = True
      postsolve_job_ID = postsolve_check[0][8:]
      
  print('Post-solve started! Waiting for it to finish...')
  postsolve_finished = False
  while(not postsolve_finished):
    
    time.sleep(30)
    command_string = 'squeue -j ' + postsolve_job_ID + ' > queue_status.txt'
    os.system(command_string)
    
    status_check = open('queue_status.txt', 'r')
    line_check = 0
    for line in status_check:
      line_check = line_check+1
    
    status_check.close()
    command_string = 'rm queue_status.txt'
    os.system(command_string)
    
    if(line_check < 1):
      postsolve_finished = True
  
  # Move the post-processed restart file into a new folder, then move it out
  command_string = 'mkdir -p fsi_presolve/'
  print(command_string)
  os.system(command_string)
  
  command_string = 'mv restart.200.0 fsi_presolve/restart.0.1'
  print(command_string)
  os.system(command_string)
  
  command_string = 'mv fsi_presolve ../../'
  print(command_string)
  os.system(command_string)
  
  command_string = 'cd ../../'
  print(command_string)
  os.chdir('../../')
  
  # TODO: Print instruction message on what to do with the restart file,
  #       how to assign material properties, etc.
  
  # TODO: Produce template .svpre file for FSI simulations
  
  # TODO: Produce README.txt with all relevant instructions
  
if __name__ == '__main__':

########################### Run Rigid Presolve #######################
  rigidPresolve()
  

########################### Run FSI Presolve ########################
  #Get the current foldername
  WorkingDirectory=os.getcwd()
  ProjectName=WorkingDirectory.split("/")[-1]
  ProjectNameCoarse=WorkingDirectory+"/../%s"%ProjectName.replace("fine","coarse")

  #Copy the coronaryLPN_FSI.svpre file to fine folder
  infile=open("%s/rigid_coronary_tuning_and_simulation/coronaryLPN_FSI.svpre"%ProjectNameCoarse,'r')
  outfile=open("./coronaryLPN_FSI.svpre",'w')
  for line in infile:
    outfile.write(line.replace("coarse/mesh-complete","fine/mesh-complete"))
  outfile.close()


  #Copy the folders from the coarse mesh
  os.system("cp %s/FSI_coronary_tuning_and_simulation/coronaryModel.txt ./"%ProjectNameCoarse)    
  os.system("cp %s/FSI_coronary_tuning_and_simulation/coronaryParams.txt ./"%ProjectNameCoarse)    
  os.system("cp %s/FSI_coronary_tuning_and_simulation/GenBC ./"%ProjectNameCoarse) 

  #Copy the presolve_batch script to fsi_presolv
  os.system("cp coronaryLPN_FSI.svpre fsi_presolve/")
  infile=open("presolve_batch",'r')
  outfile=open("fsi_presolve/presolve_batch",'w')
  for line in infile:
    outfile.write(line.replace("coronaryLPN.svpre","coronaryLPN_FSI.svpre"))
  outfile.close()

  os.chdir('./fsi_presolve')

  ########################### Run FSI Presolve ########################print(command_string)
  # Run the pre-solver
  command_string = BATCH_COMMAND + ' presolve_batch'
  os.system(command_string)

  presolve_start = False
  presolve_job_ID = '-1'
  while(not presolve_start):

    time.sleep(30)
    presolve_check = glob.glob('svPre.o*')
    if(len(presolve_check) > 0):
      presolve_start = True
      presolve_job_ID = presolve_check[0][7:]

  print('Pre-solve started! Waiting for it to finish...')
  print('jobID = ' + presolve_job_ID)
  presolve_finished = False
  while(not presolve_finished):

    time.sleep(30)
    command_string = 'squeue -j ' + presolve_job_ID + ' > queue_status.txt'
    os.system(command_string)

    status_check = open('queue_status.txt', 'r')
    line_check = 0
    for line in status_check:
      line_check = line_check+1

    status_check.close()
    command_string = 'rm queue_status.txt'
    os.system(command_string)

    if(line_check < 1):
      presolve_finished = True

  os.chdir("../")

############################ Run the Solver ##########################################

  #Read the solver.inp file from the coarse folder but now change some of the parameters
  infile=open("%s/FSI_coronary_tuning_and_simulation/solver.inp"%ProjectNameCoarse,'r')
  outfile=open("solver.inp",'w')
  for line in infile:
    #Assing the number of timesteps
    if line.find("Number of Timesteps:")>=0:
      outfile.write("Number of Timesteps: %d\n"%(NumberOfCycles*NumberOfTimesteps))
    #Assign the timestep size
    elif line.find("Time Step Size:")>=0:
      outfile.write("Time Step Size: %.05f\n"%(Period/NumberOfTimesteps))
    #Assign the Number Of Stored Points
    elif line.find("Number of Timesteps between Restarts:")>=0:
      outfile.write("Number of Timesteps between Restarts: %d\n"%SaveFrequency)
    elif line.find("Residual Criteria:")>=0:
      outfile.write("Residual Criteria: %.05f\n"%ResidualCriteria)
    else:
      outfile.write(line)
  outfile.close()


  #Copy the solver batch file to current folder 
  os.system("cp fsi_presolve/geombc.dat.1 ./")
  os.system("cp fsi_presolve/restart.0.1 ./")
  
  #Create a solver batch file
  infile=open("solver_files/Niagara_flowsolver_script",'r')
  outfile=open("Niagara_flowsolver_script",'w')
  for line in infile:
    if line.find("#SBATCH --partition=debug")>=0:
      outfile.write("#SBATCH --partition=compute\n")

    elif line.find("#SBATCH --time=")>=0:
      outfile.write("#SBATCH --time=24:00:00\n")

    elif line.find("#SBATCH --nodes=4")>=0:
      outfile.write("#SBATCH --nodes=8\n")

    else:
      outfile.write(line)
  outfile.close()
  
  #Write numberstart file
  infile=open("numstart.dat",'w')
  infile.write("0")
  infile.close()
 
  os.system(BATCH_COMMAND + " Niagara_flowsolver_script") 
