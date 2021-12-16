from time import sleep
import vtk
import os
import glob
import sys
import time
import numpy as np
import math
from scipy.optimize import minimize
pConv = 1333.34

username="ana"
# ******************* DEFINE USER INPUTS IN THIS BLOCK ************************
# ** Patient clinical information (if not known, indicate 'NONE')

heart_rate = 60 #BPM
Pao_min    = 89.8 # mmHg
Pao_max    = 151.3 # mmHg
strokeVol = 80.85 # mL
ejectFract = 0.69
sys_cor_split= 6.37 #The fraction of cardiac output to all coronary arteries

meanPressure = (0.333)*Pao_max + (0.667)*Pao_min # mmHg
Ppul_mean = 14.0 # mmHg
Qla_ratio = 'None' #1.07 # Ratio of 'early' to 'late' flows into LV (i.e. E/A wave)
mit_valve = 'None' #0.763  #0.56 # Fraction of heart cycle that mitral valve is open for
aor_valve = 'None' #0.311  #0.39 # Fraction of heart cycle that aortic valve is open for
pul_valve = 'None' #0.377  #0.374 # Fraction of heart cycle that pulmonary valve is open for
Pra_mean = 'None'  #3.0 # mmHg - IVC right atrial mean pressure
meanFlow = strokeVol*(float(heart_rate)/60.0) # mL/s
Cam_scale = 0.89
Ca_scale = 0.11
Crcr_estim = 100e-6 #compliance of the aorta (estimate from 3ewk)

# ************************* Left and Right Coronary Split *******************
#Write  file to store flow Split data
Left_Cor_split_assigned=70
Right_Cor_split_assigned=100-Left_Cor_split_assigned

f_l=(sys_cor_split/100.)*(Left_Cor_split_assigned/100.)
f_r=(sys_cor_split/100.)*(Right_Cor_split_assigned/100.)

Left_Cor_split_computed=0 #Don't touch
Aor_Cor_split_computed=0 #Don't touch

#Write the data to a file
Flow_split_file=open("LeftRightFlowSplit.dat",'w')
Flow_split_file.write("Iteration AortaCoronaryAssigned AortaCoronaryMeasured LeftSplitAssigned LeftSplitMeasured LeftResistance RightResistance\n")
Flow_split_file.close()

# *********************** Coronary Flow Split ************************
splitting_scheme="owais"
if splitting_scheme=="justin":sys_cor_split = 4.0
if splitting_scheme=="owais": sys_cor_split = (f_l+f_r)*100

# *************************** SOLVER PARAMETERS *************************
N_timesteps=200
N_cycles=3

# ** FSI parameters

UNIFORM_WALL_PROPERTIES = False
# For uniform wall properties only
deformable_Evw=10000000
deformable_thickness=0.05
# For variable wall properties only
deformable_aortic_E = 7000000
deformable_coronary_E = 11500000
deformable_graft_E = 50000000
deformable_lima_E = 7000000
# For all deformable simulations
deformable_nuvw=0.5
deformable_density=1.0
deformable_kcons=0.833333
deformable_pressure=119990

# ** Mesh information
WORKING_DIRECTORY  = os.getcwd()
MESH_SURFACES_PATH = os.getcwd()+"/mesh-complete/mesh-surfaces"
aorBranchTag = 'aorta'
lcaBranchTag = 'lca'
lcxBranchTag = 'lcx'
rcaBranchTag = 'rca'
inflowTag = 'inflow'

# Remove the cap tags if present
CapFileNames=glob.glob(MESH_SURFACES_PATH+"/cap_*.vtp")
for CapFileName in CapFileNames:
	os.system("mv %s %s"%(CapFileName, CapFileName.replace("cap_","")))

# ** Executables and file paths
GCODE_BINARY_DIR = '/home/k/khanmu11/khanmu11/Softwares/0DSolver/lpnbin/'
PRESOLVER_PATH = '/home/k/khanmu11/%s/Softwares/svSolver/BuildWithMake/Bin/svpre.exe'%username
POST_SOLVER_PATH = '/home/k/khanmu11/%s/Softwares/svSolver/BuildWithMake/Bin/svpost.exe'%username


# ** Cluster settings
FLOWSOLVER_NODES = 4
PROCESSORS=40
BATCH_COMMAND = "sbatch"
RUN_COMMAND = 'srun' # usually defined inside your batch scripts
OPT_SCRIPT = "Niagara_OPT_script"
DREAM_SCRIPT = "Niagara_DREAM_script"
SOLVER_SCRIPT = "Niagara_flowsolver_script"
RES_OPT_SCRIPT = "Niagara_optimize_surr_res_script"
COM_OPT_SCRIPT = "Niagara_optimize_surr_com_script"

# ** Other tuning settings

USER_EMAIL_ADDRESS = 'aseresti@ryerson.ca'
USE_OPTIMIZATION = True
LOG_SCALE = 0.0001
NUM_RESTARTS = 50
MAX_ITER = 200

# ************************ USER INPUTS END HERE *******************************

if(GCODE_BINARY_DIR[-1:] is not '/'):
  GCODE_BINARY_DIR = GCODE_BINARY_DIR + '/'

SURROGATE_COMMAND = GCODE_BINARY_DIR + 'bin/runCoronaryModel'
SURROGATE_OUTPUTS_PATH = 'outputTargets.out'
SOLVER_INPUT_PATH = 'solver.inp'
CORONARY_MODEL_PATH = 'coronaryModel.txt'
THREE_D_RESULTS_PATH = 'AllData'

#-------------------------------------------------------------------------------

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
  scriptFile.write('#SBATCH --partition=debug\n\n')
  scriptFile.write('# Specify the name of the output file. The %j specifies the job ID\n')
  scriptFile.write('#SBATCH --output='+str(jobName)+'.o%j\n\n')
  scriptFile.write('# Specify the name of the error file. The %j specifies the job ID\n')
  scriptFile.write('#SBATCH --error='+str(jobName)+'.e%j\n\n')
  scriptFile.write('# The walltime you require for your job\n')
  scriptFile.write('#SBATCH --time='+str(time)+':00:00\n\n')
  scriptFile.write('# Job priority. Leave as normal for now\n')
  scriptFile.write('#SBATCH --qos=normal\n\n')
  scriptFile.write('# Number of nodes are you requesting for your job. You can have 40 processors per node\n')
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

def makeNM_script():
  runScript_base(OPT_SCRIPT, 'NM_cor', 1, 1, 1)
  NM_script = open(OPT_SCRIPT, 'a')
  NM_script.write('\n')
  NM_script.write(RUN_COMMAND + ' ' + GCODE_BINARY_DIR + 'bin/tuneModel_NM 0 coronary.csv 2 coronaryParams.txt\n')
  NM_script.close()

#-------------------------------------------------------------------------------

def makeFlowsolver_script():
  runScript_base(SOLVER_SCRIPT, 'svFlowsolver', 1, FLOWSOLVER_NODES, PROCESSORS)
  flow_script = open(SOLVER_SCRIPT, 'a')
  flow_script.write('\n')
  flow_script.write(RUN_COMMAND + ' /home/k/khanmu11/%s/Softwares/svSolver/BuildWithMake/Bin/svsolver-openmpi.exe\n'%username)
  flow_script.close()

#-------------------------------------------------------------------------------

def makeSurrResOpt_script():
  runScript_base(RES_OPT_SCRIPT, 'SUR_OPT', 1, 1, 1)
  opt_script = open(RES_OPT_SCRIPT, 'a')
  opt_script.write('\n')
  opt_script.write(RUN_COMMAND + ' python optimizeSurrogateRes.py\n')
  opt_script.close()
  
#-------------------------------------------------------------------------------

def makeSurrComOpt_script():
  runScript_base(COM_OPT_SCRIPT, 'SUR_OPT', 1, 1, 1)
  opt_script = open(COM_OPT_SCRIPT, 'a')
  opt_script.write('\n')
  opt_script.write(RUN_COMMAND + ' python optimizeSurrogateCom.py\n')
  opt_script.close()

#-------------------------------------------------------------------------------

def runNM_OPT():
  print("\n** Running Nelder-Mead Optimization\n")
  #command_string = BATCH_COMMAND + ' ' + OPT_SCRIPT ############## MOK Changed for Niagara
  command_string = GCODE_BINARY_DIR + 'bin/tuneModel_NM 0 coronary.csv 2 coronaryParams.txt'
  print(command_string)
  os.system(command_string)

#-------------------------------------------------------------------------------

def runDREAM_MCMC():
  print("\n** Running DREAM Parallel MCMC\n")
  command_string = BATCH_COMMAND + ' ' + DREAM_SCRIPT
  print(command_string)
  os.system(command_string)

#-------------------------------------------------------------------------------

def run3D_SIM():
  print("\n** Starting 3D SimVascular Simulation\n")
  command_string = BATCH_COMMAND + ' ' + SOLVER_SCRIPT
  print(command_string)
  os.system(command_string)

#-------------------------------------------------------------------------------

def runSurrogateSim():
  print("** Running surrogate simulation")
  command_string = SURROGATE_COMMAND
  print(command_string)
  os.system(command_string)
  output_file = open(SURROGATE_OUTPUTS_PATH, 'r')
  surr_results = []
  for line in output_file:
    elements = line.split()
    if(len(elements) == 0):
      continue;
    res = elements[0]
    if(res == 'PaoMin' or res == 'PaoMax' or res == 'PaoMean' 
       or res == 'AorCorSplit' or res == 'AbsQin' or res == 'LCorMaxRatio'
       or res == 'LCorTotRatio' or res == 'LThirdFF' or res == 'LHalfFF'
       or res == 'RCorMaxRatio' or res == 'RCorTotRatio' 
       or res == 'RThirdFF' or res == 'RHalfFF'):
       surr_results.append(float(elements[2]))
  output_file.close()

  return surr_results

#-------------------------------------------------------------------------------

def runSurrogateResOptimization():
  print('\n** Optimizing surrogate resistances using Nelder-Mead\n')
  #command_string = BATCH_COMMAND + ' ' + RES_OPT_SCRIPT
  command_string = "python optimizeSurrogateRes.py"
  print(command_string)
  os.system(command_string)
  
#-------------------------------------------------------------------------------

def runSurrogateComplianceOptimization():
  print('\n** Optimizing surrogate compliance using Nelder-Mead\n')
  ##command_string = BATCH_COMMAND + ' ' + COM_OPT_SCRIPT
  command_string = "python optimizeSurrogateCom.py"
  print(command_string)
  os.system(command_string)

#-------------------------------------------------------------------------------

def getLL_params():
  standard_devs = []
  weights = []
  standard_devs.append(5.0)    # PaoMin
  standard_devs.append(5.0)   # PaoMax
  standard_devs.append(5.0)    # PaoMean
  standard_devs.append(0.05)    # AorCorSplit
  standard_devs.append(5.0)    # AbsQin
  standard_devs.append(0.8)    # LCorMaxRatio
  standard_devs.append(2.5337) # LCorTotRatio
  standard_devs.append(0.02)   # LThirdFF
  standard_devs.append(0.03)   # LHalfFF
  standard_devs.append(0.3)    # RCorMaxRatio
  standard_devs.append(1.0816) # RCorTotRatio
  standard_devs.append(0.07)   # RThirdFF
  standard_devs.append(0.07)   # RHalfFF

  weights.append(3.0) # PaoMin
  weights.append(3.0) # PaoMax
  weights.append(3.0) # PaoMean
  weights.append(5.0) # AorCorSplit MOK: Previously 1.0
  weights.append(5.0) # AbsQin
  weights.append(2.5) # LCorMaxRatio
  weights.append(2.5) # LCorTotRatio
  weights.append(2.5) # LThirdFF
  weights.append(2.5) # LHalfFF
  weights.append(2.5) # RCorMaxRatio
  weights.append(2.5) # RCorTotRatio
  weights.append(2.5) # RThirdFF
  weights.append(2.5) # RHalfFF





  """standard_devs.append(8.0)    # PaoMin
  standard_devs.append(12.0)   # PaoMax
  standard_devs.append(9.6)    # PaoMean
  standard_devs.append(0.05)    # AorCorSplit
  standard_devs.append(9.0)    # AbsQin
  standard_devs.append(0.8)    # LCorMaxRatio
  standard_devs.append(2.5337) # LCorTotRatio
  standard_devs.append(0.02)   # LThirdFF
  standard_devs.append(0.03)   # LHalfFF
  standard_devs.append(0.3)    # RCorMaxRatio
  standard_devs.append(1.0816) # RCorTotRatio
  standard_devs.append(0.07)   # RThirdFF
  standard_devs.append(0.07)   # RHalfFF

  weights.append(0.5) # PaoMin
  weights.append(0.25) # PaoMax
  weights.append(1.0) # PaoMean
  weights.append(1.0) # AorCorSplit MOK: Previously 1.0
  weights.append(0.25) # AbsQin
  weights.append(5.0) # LCorMaxRatio
  weights.append(5.0) # LCorTotRatio
  weights.append(5.0) # LThirdFF
  weights.append(5.0) # LHalfFF
  weights.append(5.0) # RCorMaxRatio
  weights.append(5.0) # RCorTotRatio
  weights.append(5.0) # RThirdFF
  weights.append(5.0) # RHalfFF"""

  return standard_devs, weights

#-------------------------------------------------------------------------------

def makeCoronaryBaseModel():
  filelist = glob.glob('*')
  if(not 'coronaryModel_base.txt' in filelist):
    os.system('cp coronaryModel.txt coronaryModel_base.txt')

#-------------------------------------------------------------------------------

def getSurrogateResistances():
  model_name = 'coronaryModel.txt'
  model_file = open(model_name, 'r')
  resistances = []
  for line in model_file:
    elements = line.split(',')
    if(elements[0] == 'SurrogateRes'):
      for i in range(1, len(elements)):
        if(not elements[i] == '\n'):
          resistances.append(float(elements[i]))
  model_file.close()
  return resistances

#-------------------------------------------------------------------------------

def getSurrogateCompliance():
  model_name = 'coronaryModel.txt'
  model_file = open(model_name, 'r')
  resistances = []
  for line in model_file:
    elements = line.split(',')
    if(elements[0] == 'Csurr'):
      for i in range(1, len(elements)):
        if(not elements[i] == '\n'):
          resistances.append(float(elements[i]))
  model_file.close()
  return resistances

#-------------------------------------------------------------------------------

def updateSurrogateCompliance(new_compliances):
  base_name = 'coronaryModel_base.txt'
  model_name = 'coronaryModel.txt'
  base_file = open(base_name, 'r')
  model_file = open(model_name, 'w')
  for line in base_file:
    elements = line.split(',')
    if(elements[0] == 'Csurr'):
      model_file.write('Csurr,')
      for i in range(len(new_compliances)):
        model_file.write('{0:.6f},'.format(new_compliances[i]))
      model_file.write('\n')
    else:
      model_file.write(line)
  base_file.close()
  model_file.close()

#-------------------------------------------------------------------------------

def updateSurrogateResistances(new_resistances):
  base_name = 'coronaryModel_base.txt'
  model_name = 'coronaryModel.txt'
  base_file = open(base_name, 'r')
  model_file = open(model_name, 'w')
  for line in base_file:
    elements = line.split(',')
    if(elements[0] == 'SurrogateRes'):
       model_file.write('SurrogateRes,')
       for i in range(len(new_resistances)):
          model_file.write('{0:.6f},'.format(new_resistances[i]))
       model_file.write('\n')
    else:
       model_file.write(line)
  base_file.close()
  model_file.close()

#-------------------------------------------------------------------------------

def post_process_3D_results(all_data_path):
  #Define global variables
  global Left_Cor_split_computed
  global Aor_Cor_split_computed

  total_steps = -1
  single_cycle = -1
  step_size_3D = -1

  solver_file = open(SOLVER_INPUT_PATH, 'r')
  for line in solver_file:
    line_split = line.split()
    if(len(line_split) > 2):
      if(line_split[2] == 'Timesteps:'):
        total_steps = int(line_split[3])
    if(len(line_split) > 1):
      if(line_split[1] == 'Step'):
        step_size_3D = float(line_split[3])
  solver_file.close()

  if(total_steps == -1 or step_size_3D == -1):
    print('Error: Could not process solver.inp. Please ensure that it is formatted correctly')
    exit(1)

  model_file = open(CORONARY_MODEL_PATH, 'r')
  line = model_file.readline()
  elements = line.split(',')
  nCOR_l = int(elements[0])
  nCOR_r = int(elements[1])
  nRCR = int(elements[2])
  Tc = float(elements[3])
  model_file.close()

  nUnknowns = 2*(nCOR_l + nCOR_r) + nRCR + 10
  nFaces = nCOR_l + nCOR_r + nRCR
  single_cycle = int(round(float(Tc/step_size_3D)/10.0)*10.0)

  all_data_file = open(all_data_path, 'r')
  line = all_data_file.readline()
  elements = line.split()
  AllData = np.zeros((total_steps, len(elements)))
  count = 0
  for el in elements:
    AllData[0, count] = float(el)
    count = count + 1

  outer_count = 1
  for line in all_data_file:
    elements = line.split()
    count = 0
    for el in elements:
      AllData[outer_count, count] = float(el)
      count = count+1
    outer_count = outer_count + 1
  all_data_file.close()
   
  auxStart = nUnknowns
  rcr_st = auxStart + 4
  l_cor_st = rcr_st + nRCR
  r_cor_st = l_cor_st + nCOR_l
  ao_valve = auxStart + nFaces + 12
  mit_valve = auxStart + nFaces + 11
  tri_valve = auxStart + nFaces + 13
  pulm_valve = auxStart + nFaces + 14

# Begin post-processing of 3D results here
  # SUM RCR FLUX
  temp = 0.0
  Q_rcr = 0.0
  for loopA in range(nRCR):
    temp = np.trapz(AllData[total_steps - single_cycle : total_steps, loopA+rcr_st], x=AllData[total_steps - single_cycle : total_steps, auxStart])
    Q_rcr = Q_rcr + temp
   
  # SUM LEFT CORONARY FLUX
  Q_lcor = 0.0
  for loopA in range(nCOR_l):
    temp = np.trapz(AllData[total_steps - single_cycle : total_steps, loopA+l_cor_st], x=AllData[total_steps - single_cycle : total_steps, auxStart])
    Q_lcor = Q_lcor + temp

  # INTEGRATE LEFT MAIN FLOW
  lmain_flow = np.trapz(AllData[total_steps - single_cycle : total_steps, l_cor_st], x=AllData[total_steps - single_cycle : total_steps, auxStart])

  # SUM RIGHT CORONARY FLUX
  Q_rcor = 0.0
  for loopA in range(nCOR_r):
    temp = np.trapz(AllData[total_steps - single_cycle : total_steps, loopA+r_cor_st], x=AllData[total_steps - single_cycle : total_steps, auxStart])
    Q_rcor = Q_rcor + temp

  # INTEGRATE RIGHT MAIN FLOW
  rmain_flow = np.trapz(AllData[total_steps - single_cycle : total_steps, r_cor_st], x=AllData[total_steps - single_cycle : total_steps, auxStart])

  # FIND THE END OF SYSTOLE
  systole_end = int((total_steps - single_cycle)/2) - 1
  for i in range(total_steps - single_cycle - 1, total_steps-1):
    if(AllData[i, ao_valve] > 0.0 and AllData[i+1, ao_valve] == 0.0):
      systole_end = i
      break

  # FIND THE START OF SYSTOLE
  systole_start = total_steps - single_cycle - 1
  for i in range(total_steps - single_cycle - 1, total_steps - 1):
    if(AllData[i, ao_valve] == 0.0 and AllData[i+1, ao_valve] > 0.0):
      systole_start = i
      break
  ao_open = systole_start

  # FIND WHEN THE MITRAL VALVE OPENS
  mit_open = total_steps - single_cycle - 1
  for i in range(total_steps - single_cycle - 1, total_steps - 1):
    if(AllData[i, mit_valve] == 0.0 and AllData[i+1, mit_valve] > 0.0):
      mit_open = i
      break

  mit_half = int( round((mit_open + total_steps)/2.0) )
  aor_half = int( round((ao_open + systole_end)/2.0) )

  # CALCULATE MAX AND TOTAL CORONARY FLOW DURING SYSTOLE
  l_cor_qmax_s = np.amax( AllData[systole_start : systole_end, l_cor_st] )
  l_cor_qtot_s = np.trapz(AllData[systole_start : systole_end, l_cor_st], x=AllData[systole_start : systole_end, auxStart])
  r_cor_qmax_s = np.amax( AllData[systole_start : systole_end, r_cor_st] )
  r_cor_qtot_s = np.trapz(AllData[systole_start : systole_end, r_cor_st], x=AllData[systole_start : systole_end, auxStart])

  # CALCULATE MAX AND TOTAL CORONARY FLOW DURING DIASTOLE
  l_cor_qmax_d = max( np.amax( AllData[systole_end : total_steps, l_cor_st] ), np.amax( AllData[total_steps - single_cycle - 1 : systole_start, l_cor_st] ) )
  l_cor_qtot_d = np.trapz(AllData[systole_end : total_steps, l_cor_st], x=AllData[systole_end : total_steps, auxStart]) + np.trapz(AllData[total_steps - single_cycle - 1 : systole_start, l_cor_st], x=AllData[total_steps - single_cycle - 1 : systole_start, auxStart])
  r_cor_qmax_d = max( np.amax( AllData[systole_end : total_steps, r_cor_st] ), np.amax( AllData[total_steps - single_cycle - 1 : systole_start, r_cor_st] ) )
  r_cor_qtot_d = np.trapz(AllData[systole_end : total_steps, r_cor_st], x=AllData[systole_end : total_steps, auxStart]) + np.trapz(AllData[total_steps - single_cycle - 1 : systole_start, r_cor_st], x=AllData[total_steps - single_cycle - 1 : systole_start, auxStart])

  # CALCULATE RATIOS (DIASTOLE TO SYSTOLE)
  l_cor_max_ratio = l_cor_qmax_d/l_cor_qmax_s
  l_cor_tot_ratio = l_cor_qtot_d/l_cor_qtot_s
  r_cor_max_ratio = r_cor_qmax_d/r_cor_qmax_s
  r_cor_tot_ratio = r_cor_qtot_d/r_cor_qtot_s

  # CALCULATE THE 1/3 FF AND 1/2 FF
  thirdCyc = int( round(single_cycle/3) )
  halfCyc = int( round(single_cycle/2) )
  if(systole_end+thirdCyc-1 < total_steps):
    r_third_FF = np.trapz( AllData[systole_end-1 : systole_end+thirdCyc, r_cor_st], x=AllData[systole_end-1 : systole_end+thirdCyc, auxStart] )/rmain_flow
    l_third_FF = np.trapz( AllData[systole_end-1 : systole_end+thirdCyc, l_cor_st], x=AllData[systole_end-1 : systole_end+thirdCyc, auxStart] )/lmain_flow
  else:
    r_third_FF = 0.0
    l_third_FF = 0.0

  if(systole_end+halfCyc-1 < total_steps):
    r_half_FF = np.trapz( AllData[systole_end-1 : systole_end+halfCyc, r_cor_st], x=AllData[systole_end-1 : systole_end+halfCyc, auxStart] )/rmain_flow
    l_half_FF = np.trapz( AllData[systole_end-1 : systole_end+halfCyc, l_cor_st], x=AllData[systole_end-1 : systole_end+halfCyc, auxStart] )/lmain_flow
  else:
    r_half_FF = 0.0
    l_half_FF = 0.0

  # COMPUTE OUTPUT QUANTITIES
  Qinlet = np.trapz( AllData[total_steps - single_cycle - 1 : total_steps, nUnknowns + nFaces + 4], x=AllData[total_steps - single_cycle - 1 : total_steps, auxStart] )
  Aor_Cor_split = ( (Q_lcor + Q_rcor) / (Q_lcor + Q_rcor + Q_rcr) ) * 100.0
  Aor_Cor_split_computed=Aor_Cor_split
  #Aor_Cor_split = ( (Q_lcor) / (Q_lcor + Q_rcor + Q_rcr) ) * 100.0
  Left_Cor_split = ( (Q_lcor) / (Q_lcor + Q_rcor))*100.0
  Left_Cor_split_computed=Left_Cor_split
  Pao_max = np.amax( AllData[total_steps - single_cycle - 1 : total_steps, 9] )
  Pao_min = np.amin( AllData[total_steps - single_cycle - 1 : total_steps, 9] )
  Pao_mean = np.mean( AllData[total_steps - single_cycle - 1 : total_steps, 9] )
  Ppul_mean = np.mean( AllData[total_steps - single_cycle - 1 : total_steps, 4] )
  EF_LV = ( np.amax( AllData[total_steps - single_cycle - 1 : total_steps, 7] ) - np.amin( AllData[total_steps - single_cycle - 1 : total_steps, 7] ) ) / np.amax( AllData[total_steps - single_cycle - 1 : total_steps, 7] )
  Prv_Pra = np.amax( AllData[total_steps - single_cycle - 1 : total_steps, auxStart + 7 + nFaces] ) - np.amax( AllData[total_steps - single_cycle - 1 : total_steps, auxStart + 9 + nFaces] )
  Ppul_Prv = np.amax( AllData[total_steps - single_cycle - 1 : total_steps, auxStart + 7 + nFaces] ) - np.amax( AllData[total_steps - single_cycle - 1 : total_steps, 4] )
  mit_valve_time = float(np.sum( AllData[total_steps - single_cycle - 1 : total_steps, mit_valve] )) / float(single_cycle)
  aor_valve_time = float(np.sum( AllData[total_steps - single_cycle - 1 : total_steps, ao_valve] )) / float(single_cycle)
  pul_valve_time = float(np.sum( AllData[total_steps - single_cycle - 1 : total_steps, pulm_valve] )) / float(single_cycle)
  Qla_ratio = np.amax( AllData[mit_open-1 : mit_half, 6] ) / np.amax( AllData[ mit_half-1 : total_steps, 6] )
  Pra_mean = np.mean( AllData[total_steps - single_cycle - 1: total_steps, auxStart + 9 + nFaces] )
  if(r_cor_max_ratio < 0 or l_cor_max_ratio < 0 or r_cor_tot_ratio < 0 or l_cor_tot_ratio < 0):
    r_cor_max_ratio = 9001.0
    l_cor_max_ratio = 9001.0
    r_cor_tot_ratio = 9001.0
    l_cor_tot_ratio = 9001.0

  # COMPUTE CONVERGENCE QUANTITIES

  # PRINT RESULT
  print("** PaoMin: \t{0:.2f}".format(Pao_min))
  print("** PaoMax: \t{0:.2f}".format(Pao_max))
  print("** PaoMean: \t{0:.2f}".format(Pao_mean))
  print("** Aor_Cor_Split: \t{0:.2f}".format(Aor_Cor_split))
  print("** Left_Cor_Split: \t{0:.2f}".format(Left_Cor_split))
  print("** Qinlet: \t{0:.2f}".format(abs(Qinlet)))
  print("** EF_LV: \t\t{0:.2f}".format(EF_LV))
  print("** Qla_ratio: \t{0:.2f}".format(Qla_ratio))
  print("** mit_valve: \t{0:.2f}".format(mit_valve_time))
  print("** aor_valve: \t{0:.2f}".format(aor_valve_time))
  print("** pul_valve: \t{0:.2f}".format(pul_valve_time))
  print("** Pra-mean: \t{0:.2f}".format(Pra_mean))
  print("** lcor_max: \t{0:.3f}".format(l_cor_max_ratio))
  print("** lcor_tot: \t{0:.3f}".format(l_cor_tot_ratio))
  print("** lcor_third: \t{0:.3f}".format(l_third_FF))
  print("** lcor_half: \t{0:.3f}".format(l_half_FF))
  print("** rcor_max: \t{0:.3f}".format(r_cor_max_ratio))
  print("** rcor_tot: \t{0:.3f}".format(r_cor_tot_ratio))
  print("** rcor_third: \t{0:.3f}".format(r_third_FF))
  print("** rcor_half: \t{0:.3f}".format(r_half_FF))

  results = [Pao_min, Pao_max, Pao_mean, Aor_Cor_split, abs(Qinlet), 
    l_cor_max_ratio, l_cor_tot_ratio, l_third_FF, l_half_FF, r_cor_max_ratio,
    r_cor_tot_ratio, r_third_FF, r_half_FF]

  return results

#-------------------------------------------------------------------------------

def evalCoronaryLL_resistance_only(resistances):
  LL = 0.0
  for res in resistances:
    LL = LL - LOG_SCALE * np.log(res)

  updateSurrogateResistances(resistances)
  surr_results = runSurrogateSim()
  true_results = post_process_3D_results('AllData')
  standard_devs, weights = getLL_params()

  for i in range(len(surr_results)):
    LL = LL + 0.5 * (true_results[i] - 
      surr_results[i])**2 / (standard_devs[i]**2 * weights[i])

  print('LL: {0:.6f}'.format(LL))
  return LL

#-------------------------------------------------------------------------------

def evalCoronaryLL_compliance_only(compliance):
  LL = 0.0
  for com in compliance:
    LL = LL - LOG_SCALE * np.log(com)
  updateSurrogateCompliance(compliance)
  # Need to get data for inflow waveform too!
  surr_results = runSurrogateSim()
  true_results = post_process_3D_results('AllData')
  standard_devs, weights = getLL_params()

  for i in range(len(surr_results)):
    LL = LL + 0.5 * (true_results[i] - surr_results[i])**2 / (standard_devs[i]**2 * weights[i])

  print('LL: {0:.6f}'.format(LL))
  return LL

#-------------------------------------------------------------------------------

def optimizeSurrogateResistances():
  resistances = np.array(getSurrogateResistances())
  for i in range(NUM_RESTARTS):
    print('\n **** STARTING NM ITERATION {0}/{1} ***** \n'.format(i, NUM_RESTARTS))
    optimized_resistances = minimize(evalCoronaryLL_resistance_only, resistances, method='nelder-mead',
      options={'xtol': 1e-8, 'disp': True, 'maxiter': MAX_ITER})
    resistances = optimized_resistances.x
  print(optimized_resistances.x)
  updateSurrogateResistances(optimized_resistances.x)
  return optimized_resistances.x

#-------------------------------------------------------------------------------

def optimizeSurrogateCompliance():
  surrogate_compliance = np.array(getSurrogateCompliance())
  for i in range(NUM_RESTARTS):
    print('\n **** STARTING NM ITERATION {0}/{1} ***** \n'.format(i, NUM_RESTARTS))
    optimized_compliance = minimize(evalCoronaryLL_compliance_only,  surrogate_compliance, method='nelder-mead', options={'xtol':1e-8, 'disp': True, 'maxiter': MAX_ITER})
    surrogate_compliance = optimized_compliance.x
  print(optimized_compliance.x)
  updateSurrogateCompliance(optimized_compliance.x)
  return optimized_compliance.x 

#-------------------------------------------------------------------------------

def coronaryTuningIteration(count, rigidFlag):

  if(rigidFlag and count == 3):
    return
  if(not rigidFlag and count == 3):
    return

  if(USE_OPTIMIZATION):
    end_name = "NM_fin.txt"
    param_name = "optParams.txt"
    if(not os.path.isfile(end_name)):
      runNM_OPT()
    '''# Check for weird Sherlock error where memory writing fails
    time.sleep(180)
    error_file = glob.glob('*.e*')
    sim_start = False
    if(len(error_file) != 0 and os.stat(error_file[0]).st_size == 0.0):
      sim_start = True
    restart_counter = 0
    while(not sim_start):
      output_file = glob.glob('*.o*')
      error_file = glob.glob('*.e*')
    
      if(len(error_file) == 0):
        print('** Tuning still in the queue. Waiting...')
        time.sleep(600)
      elif(os.stat(error_file[0]).st_size != 0.0):
        if(restart_counter > 4):
          print('Too many retries submitting the job. Aborting.')
          sys.exit(0)
        print('** Sherlock error writing files. Submitting job again...\n')
        command_string = "rm -r " + output_file[0] + " " + error_file[0] + " " + sim_folder
        print(command_string)
        os.system(command_string)
        runNM_OPT()
        time.sleep(300)
        restart_counter = restart_counter + 1
      else:
        sim_start = True'''
  else: # Using DREAM MCMC
    runDREAM_MCMC()
    end_name = "bestParams.txt"
    param_name = end_name
    # Check for weird Sherlock error where memory writing fails
    time.sleep(180)
    error_file = glob.glob('*.e*')
    sim_start = False
    if(len(error_file) != 0 and os.stat(error_file[0]).st_size == 0.0):
      sim_start = True
    restart_counter = 0
    while(not sim_start):
      output_file = glob.glob('*.o*')
      error_file = glob.glob('*.e*')
    
      if(len(error_file) == 0):
        print('** Tuning still in the queue. Waiting...')
        time.sleep(600)
      elif(os.stat(error_file[0]).st_size != 0.0):
        if(restart_counter > 4):
          print('Too many retries submitting the job. Aborting.')
          sys.exit(0)
        print('** Sherlock error writing files. Submitting job again...\n')
        command_string = "rm -r " + output_file[0] + " " + error_file[0] + " " + sim_folder
        print(command_string)
        os.system(command_string)
        runDREAM_MCMC()
        time.sleep(300)
        restart_counter = restart_counter + 1
      else:
        sim_start = True

  # Wait for optimization or DREAM to finish
  #print('\n** Waiting for tuning round ' + str(count) + ' to finish...\n')
  #while(not os.path.isfile(end_name)):
  #  time.sleep(60);
  print('** Tuning complete!\n')

  command_string = "rm " + end_name
  print(command_string)
  os.system(command_string)
  print('\n** Transforming parameter file to coronaryParams.txt\n')
  command_string = "cp " + param_name + " coronaryParams.txt"
  print(command_string)
  os.system(command_string)

  if(rigidFlag):
    folder_name = "Rigid_tuning_" + str(count)
  else:
    folder_name = 'FSI_tuning_' + str(count)

  print("\n** Cleaning up tuning files and saving them in: " + folder_name + '\n')
  if(not os.path.isdir(folder_name)):
    command_string = "mkdir " + folder_name
    print(command_string)
    os.system(command_string)
    command_string = "mv optParams.txt optValue.txt outputTargets.out " + folder_name
    print(command_string)
    os.system(command_string)
    '''output_file = glob.glob('*.o*')
    error_file = glob.glob('*.e*')
    command_string = "mv " + output_file[0] + " " + error_file[0] + " " + folder_name
    print(command_string)
    os.system(command_string)''' ########### MOK  Changed for Niagara


  # Now, all files should be ready to run a 3D simulation, so just call the flowsolver script!
  procs_list = glob.glob('*-procs_*')
  if(len(procs_list) == 0):
    run3D_SIM()

  # NEED FLAG HERE TO DETERMINE STOP OF 3D SIMULATION
  total_steps = N_timesteps*N_cycles
  procs_list = glob.glob('%d-procs_case'%(int(FLOWSOLVER_NODES*PROCESSORS)))
  print('\n** Waiting for simulation job to launch...\n')

  #Check to see if the simulation has started by finding the results folder
  while(len(procs_list) == 0):
    time.sleep(180)
    procs_list = glob.glob('%d-procs_case'%(int(FLOWSOLVER_NODES*PROCESSORS)))
    print('\n** Waiting for simulation job to launch...\n')
  sim_folder = '%d-procs_case'%(int(FLOWSOLVER_NODES*PROCESSORS))
  
  print('\n** Simulation has launched! Waiting for simulation to finish...\n')

  #Check to see if the AllData file has been written inside the results folder
  AllData_file=glob.glob('%d-procs_case/AllData'%(int(FLOWSOLVER_NODES*PROCESSORS)))
  while(len(AllData_file) == 0):
    time.sleep(180)
    AllData_file = glob.glob('%d-procs_case/AllData'%(int(FLOWSOLVER_NODES*PROCESSORS)))
    print('\n** AllData file has been written...\n')

  #Check to see if the Simulation is finished
  sim_finished = False
  while(not sim_finished):
    time.sleep(100)
    AllData_check = open(sim_folder + '/AllData', 'r')
    counter = 0
    for line in AllData_check:
      counter = counter + 1
    print ("------- Finished %d of %d timesteps"%(counter,total_steps))
    AllData_check.close()
    if(counter == total_steps):
      sim_finished = True
  print('** Simulation has finished! Moving on to post processing...\n')

  # First, copy the 3D results AllData to the current direction
  AllData_location = sim_folder + '/AllData'
  command_string = 'cp ' + AllData_location + ' .'
  print(command_string)
  os.system(command_string)

  # Make a base copy of coronaryModel
  makeCoronaryBaseModel()

  # After running 3D sim, check to see if the results match the surrogate
  # results. If they do, then break out of this function
  if(rigidFlag):
    surrogate_resistances = np.array(getSurrogateResistances())
    LL = evalCoronaryLL_resistance_only(surrogate_resistances)
    if(abs(LL) < 2.0):
      print('\n** Rigid tuning has converged! Cleaning up files...\n')
      output_file = glob.glob('*.o*')
      error_file = glob.glob('*.e*')
      command_string = "mv AllData " + output_file[0] + " " + error_file[0] + " " + folder_name
      print(command_string)
      os.system(command_string)
      return
  else:
    surrogate_compliance  = getSurrogateCompliance()
    LL = evalCoronaryLL_compliance_only(surrogate_compliance)
    if(abs(LL) < 2.0):
      print('\n** FSI tuning has converged! Cleaning up files...\n')
      output_file = glob.glob('*.o*')
      error_file = glob.glob('*.e*')
      command_string = "mv AllData " + output_file[0] + " " + error_file[0] + " " + folder_name
      print(command_string)
      os.system(command_string)
      return
  
  # Move 3D simulation files out of the folder before optimizing surrogate  
  output_file = glob.glob('*.o*')
  error_file = glob.glob('*.e*')
  command_string = "mv " + output_file[0] + " " + error_file[0] + " " + folder_name
  print(command_string)
  os.system(command_string)
  
  # Now optimize coronary surrogate resistances
  if(rigidFlag):
    runSurrogateResOptimization()
    end_name = "surr_res_fin.txt"
    ### Check for weird Sherlock error where memory writing fails
    ##time.sleep(180)
    ##error_file = glob.glob('*.e*')
    ##sim_start = False
    ##if(len(error_file) != 0 and os.stat(error_file[0]).st_size == 0.0):
    ##  sim_start = True
    ##restart_counter = 0
    ##while(not sim_start):
    
    ##  output_file = glob.glob('*.o*')
    ##  error_file = glob.glob('*.e*')
    
    ##  if(len(error_file) == 0):
    ##    print('** Simulation still in the queue. Waiting...')
    ##    time.sleep(600)
    ##  elif(os.stat(error_file[0]).st_size != 0.0):
    ##    if(restart_counter > 4):
    ##      print('Too many retries submitting the job. Aborting.')
    ##      sys.exit(0)
    ##    print('** Sherlock error writing files. Submitting job again...\n')
    ##    command_string = "rm -r " + output_file[0] + " " + error_file[0] + " " + sim_folder
    ##    print(command_string)
    ##   os.system(command_string)
    ##    runSurrogateResOptimization()
    ##    time.sleep(300)
    ##    restart_counter = restart_counter + 1
    ##  else:
    ##    sim_start = True
  else:
    runSurrogateComplianceOptimization()
    end_name = "surr_com_fin.txt"
    ### Check for weird Sherlock error where memory writing fails
    ##time.sleep(180)
    ##error_file = glob.glob('*.e*')
    ##sim_start = False
    ##if(len(error_file) != 0 and os.stat(error_file[0]).st_size == 0.0):
    ##  sim_start = True
    ##restart_counter = 0
    ##while(not sim_start):
    
    ##  output_file = glob.glob('*.o*')
    ##  error_file = glob.glob('*.e*')
    
    ##  if(len(error_file) == 0):
    ##    print('** Simulation still in the queue. Waiting...')
    ##    time.sleep(600)
    ##  elif(os.stat(error_file[0]).st_size != 0.0):
    ##    if(restart_counter > 4):
    ##      print('Too many retries submitting the job. Aborting.')
    ##      sys.exit(0)
    ##    print('** Sherlock error writing files. Submitting job again...\n')
    ##    command_string = "rm -r " + output_file[0] + " " + error_file[0] + " " + sim_folder
    ##    print(command_string)
    ##    os.system(command_string)
    ##    runSurrogateComplianceOptimization()
    ##    time.sleep(300)
    ##    restart_counter = restart_counter + 1
    ##  else:
    ##    sim_start = True

  ### Wait for optimization or DREAM to finish
  ##print('\n** Waiting for surrogate optimization round ' + str(count) + ' to finish...\n')
  ##while(not os.path.isfile(end_name)):
  ##  time.sleep(60);
  ##print('** Surrogate optimization complete!\n')
  ##command_string = "rm " + end_name
  ##print(command_string)
  ##os.system(command_string)

  # Clean up files and get ready for the next iteration
  if(count < 2):
    print('\n** Moving simulation files to ' + folder_name + '\n')
    command_string = 'mv ' + sim_folder + ' ' + folder_name
    print(command_string)
    os.system(command_string)
    ##output_file = glob.glob('*.o*')
    ##error_file = glob.glob('*.e*')
    ##command_string = "mv AllData " + output_file[0] + " " + error_file[0] + " " + folder_name
    command_string = "mv AllData " +  folder_name
    print(command_string)
    os.system(command_string)

  # Recursively call this function again until rigid tuning satisfied
  coronaryTuningIteration(count+1, rigidFlag)
  
#-------------------------------------------------------------------------------

def rigidPresolve():
  input_path = MESH_SURFACES_PATH
  if(input_path[-1:] is not '/'):
    input_path = input_path + '/'
  if(input_path[-1:] is not '*'):
    input_path = input_path + '*'
  Tc = 60.0/heart_rate

  filelist_raw = glob.glob(input_path)
  filelist_raw.sort()

  filelist = []
  for trial in filelist_raw:
    if(trial[-4:] == ".vtp"):
      filelist.append(trial)

  base_len = len(input_path)
  A_br = []
  A_RCA = []
  A_LCA = []
  A_LCX = []
  chosen_positions = []
  chosen_names = []
  inlet_position = -1
  inflow_name = ''

  # Sort through the input file list, and save the different face areas into
  # the appropriate container depending on the face name
  for vtp_file in filelist:
    temp_name = vtp_file[base_len-1:]
    if(temp_name[:len(aorBranchTag)] == aorBranchTag):
      A_br.append(findVTPArea(vtp_file))
      chosen_positions.append(findVTPMassCenter(vtp_file))
      chosen_names.append(vtp_file)
    elif(temp_name[:len(lcaBranchTag)] == lcaBranchTag):
      A_LCA.append(findVTPArea(vtp_file))
      chosen_positions.append(findVTPMassCenter(vtp_file))
      chosen_names.append(vtp_file)
    elif(temp_name[:len(lcxBranchTag)] == lcxBranchTag):
      A_LCX.append(findVTPArea(vtp_file))
      chosen_positions.append(findVTPMassCenter(vtp_file))
      chosen_names.append(vtp_file)
    elif(temp_name[:len(rcaBranchTag)] == rcaBranchTag):
      A_RCA.append(findVTPArea(vtp_file))
      chosen_positions.append(findVTPMassCenter(vtp_file))
      chosen_names.append(vtp_file)
    elif(temp_name[:len(inflowTag)] == inflowTag):
      inlet_position = findVTPMassCenter(vtp_file)
      inflow_name = vtp_file
  # Append all lcx faces to the end of the lca area list
  for lcx_area in A_LCX:
    A_LCA.append(lcx_area)

  # Calculate and store the distances between the inflow and the outflows
  distance_list = []
  for pos in chosen_positions:
    x_diff = inlet_position[0] - pos[0]
    y_diff = inlet_position[1] - pos[1]
    z_diff = inlet_position[2] - pos[2]
    dist = (x_diff**2 + y_diff**2 + z_diff**2)**(0.5)
    distance_list.append(dist)

  # Calculate and store surrogate resistances/inductances based on geometry
  model3dres = []
  L_surr = []
   
  rikka = 0
  for br_area in A_br:
    temp_res = 8.0*(distance_list[rikka])*0.04 / (math.pi * (br_area/math.pi)**2 )
    temp_res = min(temp_res, 25000)
    model3dres.append(temp_res)
    L_surr.append(5.25)
    rikka = rikka + 1

  for lca_area in A_LCA:
    temp_res = 8.0*distance_list[rikka]*0.04 / (math.pi * (lca_area/math.pi)**2)
    temp_res = min(temp_res, 25000)
    model3dres.append(temp_res)
    L_surr.append(10)
    rikka = rikka + 1

  for rca_area in A_RCA:
    temp_res = 8.0*distance_list[rikka]*0.04 / (math.pi * (rca_area/math.pi)**2)
    temp_res = min(temp_res, 25000)
    model3dres.append(temp_res)
    L_surr.append(10)
    rikka = rikka + 1

  # Assign the number of outlets, number of unknowns, types of outlets, and faceToStateMapping
  nCOR_l = len(A_LCA)
  nCOR_r = len(A_RCA)
  nCOR = nCOR_l + nCOR_r
  nRCR = len(A_br)
  nFaces = nCOR + nRCR + 1
  nUnknowns = 10 + nRCR + 2*nCOR
  face_counter = 10
  faceToStateMapping = []
  for i in range(nRCR):
    faceToStateMapping.append(face_counter)
    face_counter = face_counter + 1

  for i in range(nCOR):
    faceToStateMapping.append(face_counter)
    face_counter = face_counter + 2

  faceToStateMapping.append(9)

  # Now to calculate and assign outlet resistances. First, calculate Murray's Law scaling factors
  r_br_murray = []
  r_RCA_murray = []
  r_LCA_murray = []
  RA_rcr = 0.0
  RA_cor = 0.0

  for br_area in A_br:
    temp = (math.sqrt(br_area/math.pi))**(2.0)
    RA_rcr = RA_rcr + temp
    r_br_murray.append(temp)

  for rca_area in A_RCA:
    temp = (math.sqrt(rca_area/math.pi))**(2.66)
    RA_cor = RA_cor + temp
    r_RCA_murray.append(temp)

  for lca_area in A_LCA:
    temp = (math.sqrt(lca_area/math.pi))**(2.66)
    RA_cor = RA_cor + temp
    r_LCA_murray.append(temp)

  if splitting_scheme=="justin":
    R_tot = meanPressure * pConv / meanFlow
    beta = sys_cor_split*0.01 / (1.0 - sys_cor_split*0.01)
    R_cor = R_tot * (1.0 + beta)/beta
    R_br = R_tot * (1.0 + beta)
  elif splitting_scheme=="owais":
    R_tot  = (meanPressure * pConv)/meanFlow
    R_cor_l= R_tot/float(f_l) #flow split to the left coronary tree
    R_cor_r= R_tot/float(f_r) #flow split to the right coronary tree
    R_aorta= R_tot/(1.-f_l-f_r) #flow split to the aorta

  Rrcr_base = []
  Crcr_base = []

  # Calculate and store the systemic outlet resistances and capacitances
  for i in range(nRCR):
    if splitting_scheme=="justin": temp =     R_br*RA_rcr/r_br_murray[i]
    if splitting_scheme=="owais":  temp =  R_aorta*RA_rcr/r_br_murray[i]
    Rrcr_base.append(temp / pConv)
    Crcr_base.append((A_br[i] * Crcr_estim) / sum(A_br))

  Ram_l_base = []
  Rv_l_base = []
  Cam_l_base = []
  Ca_l_base = []
  C_lca = 0.5e-6 #MOK original was 1e-5

  # Calculate and store the left coronary outlet resistances and capacitances
  for i in range(nCOR_l):
    if splitting_scheme=="justin":temp = RA_cor*R_cor/r_LCA_murray[i]
    if splitting_scheme=="owais": temp = RA_cor*R_cor_l/r_LCA_murray[i]
    Ram_l_base.append(0.89 * temp / pConv)
    Rv_l_base.append(0.11 * temp / pConv)
    CA = C_lca / sum(A_LCA)
    Cam_l_base.append(Cam_scale * A_LCA[i] * CA * pConv)
    Ca_l_base.append(Ca_scale * A_LCA[i] * CA * pConv)

  Ram_r_base = []
  Rv_r_base = []
  Cam_r_base = []
  Ca_r_base = []
  C_rca = 0.25e-6 #MOK: original was 1e-5
   
  # Calculate and store the right coronary outlet resistances and capacitances
  for i in range(nCOR_r):
    if splitting_scheme=="justin": temp = RA_cor*R_cor/r_RCA_murray[i]
    if splitting_scheme=="owais":  temp = RA_cor*R_cor_r/r_RCA_murray[i]
    Ram_r_base.append(0.89 * temp / pConv)
    Rv_r_base.append(0.11 * temp / pConv)
    CA = C_rca / sum(A_RCA)
    Cam_r_base.append(Cam_scale * A_RCA[i] * CA * pConv)
    Ca_r_base.append(Ca_scale * A_RCA[i] * CA * pConv)

  # Need to scale the rcr surrogate resistances so surrogate simulation
  # will converge. Find smallest combinatin of Crcr * R_surrogate
  min_combo = 9999999.0
  for i in range(len(Crcr_base)):
    min_combo = min(min_combo, Crcr_base[i]*model3dres[i])

  # Scale this resistance such that a simulation with 10000 time steps per
  # cycle, and a heart rate of 60 BPM will converge
  scale = Tc * pConv / (20000 * min_combo)
  if(scale >= 1.0):
    for i in range(len(Crcr_base)):
      temp_res = scale * model3dres[i]
      model3dres[i] = min(temp_res, 25000)

  # Write out to model file for use with gcode framework
  coronaryModel = open('coronaryModel.txt','w')
  coronaryModel.write('{0},{1},{2},{3:.4f}\n'.format(nCOR_l,nCOR_r,nRCR,Tc))

  temp_string = 'FaceToStateMapping'
  for mapper in faceToStateMapping:
    temp_string = temp_string + ',{0}'.format(mapper)
  temp_string = temp_string + '\n'
  coronaryModel.write(temp_string)

  temp_string = 'Rrcr'
  for r_rcr in Rrcr_base:
    temp_string = temp_string + ',{0:.8f}'.format(r_rcr)
  temp_string = temp_string + '\n'
  coronaryModel.write(temp_string)

  temp_string = 'Crcr'
  for c_rcr in Crcr_base:
    temp_string = temp_string + ',{0:.8f}'.format(c_rcr)
  temp_string = temp_string + '\n'
  coronaryModel.write(temp_string)

  temp_string = 'Ram_l'
  for r_am_l in Ram_l_base:
    temp_string = temp_string + ',{0:.8f}'.format(r_am_l)
  temp_string = temp_string + '\n'
  coronaryModel.write(temp_string)

  temp_string = 'Rv_l'
  for r_v_l in Rv_l_base:
    temp_string = temp_string + ',{0:.8f}'.format(r_v_l)
  temp_string = temp_string + '\n'
  coronaryModel.write(temp_string)

  temp_string = 'Cam_l'
  for c_am_l in Cam_l_base:
    temp_string = temp_string + ',{0:.8f}'.format(c_am_l)
  temp_string = temp_string + '\n'
  coronaryModel.write(temp_string)

  temp_string = 'Ca_l'
  for c_a_l in Ca_l_base:
    temp_string = temp_string + ',{0:.8f}'.format(c_a_l)
  temp_string = temp_string + '\n'
  coronaryModel.write(temp_string)

  temp_string = 'Ram_r'
  for r_am_r in Ram_r_base:
    temp_string = temp_string + ',{0:.8f}'.format(r_am_r)
  temp_string = temp_string + '\n'
  coronaryModel.write(temp_string)

  temp_string = 'Rv_r'
  for r_v_r in Rv_r_base:
    temp_string = temp_string + ',{0:.8f}'.format(r_v_r)
  temp_string = temp_string + '\n'
  coronaryModel.write(temp_string)

  temp_string = 'Cam_r'
  for c_am_r in Cam_r_base:
    temp_string = temp_string + ',{0:.8f}'.format(c_am_r)
  temp_string = temp_string + '\n'
  coronaryModel.write(temp_string)

  temp_string = 'Ca_r'
  for c_a_r in Ca_r_base:
    temp_string = temp_string + ',{0:.8f}'.format(c_a_r)
  temp_string = temp_string + '\n'
  coronaryModel.write(temp_string)

  temp_string = 'SurrogateRes'
  for surr_res in model3dres:
    temp_string = temp_string + ',{0:.8f}'.format(surr_res)
  temp_string = temp_string + '\n'
  coronaryModel.write(temp_string)

  temp_string = 'L_surr'
  for surr_induct in L_surr:
    temp_string = temp_string + ',{0:.8f}'.format(surr_induct)
  temp_string = temp_string + '\n'
  coronaryModel.write(temp_string)

  coronaryModel.close()

  # Generate an .svpre script that is consistent with the ordering used above
  svpre_file = open('coronaryLPN.svpre','w')
   
  mesh_complete = input_path[:-1] + '../' + 'mesh-complete.mesh.vtu'
  walls_combined = input_path[:-1] + '../' + 'walls_combined.vtp'
  mesh_exterior = input_path[:-1] + '../' + 'mesh-complete.exterior.vtp'
  zp_string = 'zero_pressure_vtp '
  id_string = 'set_surface_id_vtp '

  svpre_file.write('mesh_and_adjncy_vtu ' + mesh_complete + '\n')
  svpre_file.write('noslip_vtp ' + walls_combined + '\n\n')

  for name in chosen_names:
    temp_string = zp_string + name + '\n'
    svpre_file.write(temp_string)

  svpre_file.write(zp_string + inflow_name + '\n\n')
  temp_string = id_string + mesh_exterior + ' 1\n'
  svpre_file.write(temp_string)
  count = 2

  for name in chosen_names:
    temp_string = id_string + name + ' ' + str(count) + '\n'
    svpre_file.write(temp_string)
    count = count + 1

  svpre_file.write(id_string + inflow_name + ' ' + str(count) + '\n\n')

  svpre_file.write('write_geombc geombc.dat.1\n')
  svpre_file.write('initial_pressure 119990\n')
  svpre_file.write('write_restart restart.0.1')

  svpre_file.close()


  # Generate a solver.inp that can be used for closed-loop simulations
  solver_writer = open('solver.inp','w')

  solver_writer.write('# Phasta Version 1.5 Input File\n')
  solver_writer.write('# Produced by generateCoronaryLPN\n')
  solver_writer.write('# Specialized for closed-loop coronary simulations\n\n')

  solver_writer.write('# ****** SOLUTION CONTROL ****** #\n')
  solver_writer.write('Equation of State: Incompressible\n')
  solver_writer.write('Number of Timesteps: %d\n'%(N_timesteps*N_cycles))
  solver_writer.write('Time Step Size: {0:.6f}\n\n'.format(Tc/float(N_timesteps)))

  solver_writer.write('# ***** FLUID PROPERTIES ***** #\n')
  solver_writer.write('Viscosity: 0.04\n')
  solver_writer.write('Density: 1.06\n\n')

  solver_writer.write('# ***** OUTPUT CONTROL ***** #\n')
  solver_writer.write('Number of Timesteps between Restarts: %d\n'%N_timesteps)
  solver_writer.write('Print ybar: True\n\n')

  solver_writer.write('# ****** CARDIOVASCULAR MODELING PARAMETERS ***** #\n')
  solver_writer.write('Time Varying Boundary Conditions From File: False\n')
  solver_writer.write('Number of Coupled Surfaces: {0}\n'.format(nFaces))
  solver_writer.write('Pressure Coupling: Implicit\n')
  solver_writer.write('Number of Resistance Surfaces: 0\n')
  solver_writer.write('Find the GenBC Inside the Running Directory: True\n')
  solver_writer.write('Number of Timesteps for GenBC Initialization: 0\n')
  solver_writer.write('Number of Neumann Surfaces: {0}\n'.format(nFaces))
  neumann_list = 'List of Neumann Surfaces: '
  for i in range(2, nFaces+2, 1):
    neumann_list = neumann_list + str(i) + ' '
  solver_writer.write(neumann_list + '\n\n')

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

  # Make a coronaryParams.txt file with default parameters

  # Generate numstart.dat
  numstart_writer = open('numstart.dat', 'w')
  numstart_writer.write('0')
  numstart_writer.close()

  # Run pre-solver to generate geombc.dat.1 and restart.0.1
  print('\nRunning SimVascular Pre-solver...\n')
  command_string = PRESOLVER_PATH + ' coronaryLPN.svpre'
  print(command_string)
  os.system(command_string)

  # Make a starting coronaryParams.txt with default values to initialize NM optimization
  paramsFile = open('coronaryParams.txt','w')
  paramsFile.write('0.4\n') # Tsa
  paramsFile.write('8.5\n') # tpwave
  paramsFile.write('0.5\n') # Erv
  paramsFile.write('1.2\n') # Elv
  paramsFile.write('5.5\n') # Elvp
  paramsFile.write('0.55\n') # Lrv_a
  paramsFile.write('50.0\n') # Rrv_a
  paramsFile.write('1.0\n') # Lra_v
  paramsFile.write('5.0\n') # Rra_v
  paramsFile.write('0.002\n') # Lla_v
  paramsFile.write('5.0\n') # Rla_v
  paramsFile.write('100.0\n') # Rlv_ao
  paramsFile.write('1.0\n') # Llv_a
  paramsFile.write('0.0\n') # Vrv_u
  paramsFile.write('0.0\n') # Vlv_u
  paramsFile.write('150.0\n') # Rpd
  paramsFile.write('3.0\n') # Cp
  paramsFile.write('0.9\n') # Cpa
  paramsFile.write('4.0\n') # Kxp_ra
  paramsFile.write('0.005\n') # Kxv_ra
  paramsFile.write('0.2\n') # Emax_ra
  paramsFile.write('0.0\n') # Vaso_ra
  paramsFile.write('8.0\n') # Kxp_la
  paramsFile.write('0.008\n') # Kxv_la
  paramsFile.write('0.3\n') # Emax_la
  paramsFile.write('0.0\n') # Vaso_la
  paramsFile.write('1.0\n') # Ram_cor
  paramsFile.write('1.0\n') # Rv_cor
  paramsFile.write('1.01\n') # Cam_l
  paramsFile.write('0.1\n') # Ca_l
  paramsFile.write('1.01\n') # Cam_r
  paramsFile.write('0.25\n') # Ca_r
  paramsFile.write('1.0\n') # Rrcr
  paramsFile.write('0.75\n') # Crcr
  paramsFile.write('0.75\n') # dPdr_r
  paramsFile.close()

  # Generate a coronary.csv file that has the target information
  targetsFile = open('coronary.csv','w')
  targetsFile.write('PaoMin,Minimum Aortic Pressure,'+str(Pao_min)+'\n')
  targetsFile.write('PaoMin_conv,Minimum Aortic Pressure Convergence,0\n')
  targetsFile.write('PaoMax, Maximum Aortic Pressure,'+str(Pao_max)+'\n')
  targetsFile.write('PaoMax_conv,Maximum Aortic Pressure Convergence,0\n')
  targetsFile.write('PaoMean,Mean Aortic Pressure,'+str(meanPressure)+'\n')
  targetsFile.write('PaoMean_conv,Mean Aortic Pressure Convergence,0\n')
  targetsFile.write('AorCorSplit,Coronary Flow Split,'+str(sys_cor_split)+'\n')
  targetsFile.write('AbsQin,Inlet flow rate,'+str(meanFlow)+'\n')
  targetsFile.write('AbsQin_conv,Inlet flow rate convergence,0\n')
  targetsFile.write('Qsystole_perc,Percentage of cardiac output in systole,0.9\n')
   
  out_string = 'PpulMean,Mean Pulmonary Pressure,'
  if(Ppul_mean == 'NONE'):
    out_string = out_string + 'none\n'
  else:
    out_string = out_string + str(Ppul_mean) + '\n'
  targetsFile.write(out_string)

  out_string = 'EFLV,Left Ventricular Ejection Fraction,'
  if(ejectFract == 'NONE'):
    out_string = out_string + 'none\n'
  else:
    out_string = out_string +str(ejectFract)+'\n'
  targetsFile.write(out_string)

  out_string = 'QlaRatio,Left atrioventricular early/late flow ratio,'
  if(Qla_ratio == 'NONE'):
    out_string = out_string + 'none\n'
  else:
    out_string = out_string +str(Qla_ratio)+'\n'
  targetsFile.write(out_string)

  out_string = 'mitValveTime,Mitral Valve Opening Time,'
  if(mit_valve == 'NONE'):
    out_string = out_string + 'none\n'
  else:
    out_string = out_string +str(mit_valve)+'\n'
  targetsFile.write(out_string)

  out_string = 'aorValveTime,Aortic Valvee Opening Time,'
  if(aor_valve == 'NONE'):
    out_string = out_string + 'none\n'
  else:
    out_string = out_string +str(aor_valve)+'\n'
  targetsFile.write(out_string)

  out_string = 'pulValveTime,Pulmonary Valve Opening Time,'
  if(pul_valve == 'NONE'):
    out_string = out_string + 'none\n'
  else:
    out_string = out_string +str(pul_valve)+'\n'
  targetsFile.write(out_string)

  out_string = 'PraMean,Mean Right Atrial Pressure,'
  if(Pra_mean == 'NONE'):
    out_string = out_string + 'none\n'
  else:
    out_string = out_string +str(Pra_mean)+'\n'
  targetsFile.write(out_string)
  targetsFile.write('LCorMaxRatio,Diastolic/Systolic maximum flow ratio in left coronary,2.4\n')
  targetsFile.write('LCorTotRatio,Diastolic/Systolic total flow ratio in left coronary,4.667\n')
  targetsFile.write('LThirdFF,Left coronary 1/3 flow fraction,0.46\n')
  targetsFile.write('LHalfFF,Left coronary 1/2 flow fraction,0.67\n')
  targetsFile.write('LGradOK,Checking amount of peaks/valleys in left coronary wave,5.0\n')
  targetsFile.write('RCorMaxRatio,Diastolic/Systolic maximum flow ratio in right coronary,1.4\n')
  targetsFile.write('RCorTotRatio,Diastolic/Systolic total flow ratio in right coronary,1.2857\n')
  targetsFile.write('RThirdFF,Right coronary 1/3 flow fraction,0.39\n')
  targetsFile.write('RHalfFF,Right coronary 1/2 flow fraction,0.57\n')
  targetsFile.write('RGradOK,Checking amount of peaks/valleys in right coronary wave,5.0\n')
  targetsFile.close()

  # Move the simulation files into a seperate directory
  rigid_sim_folder = 'rigid_coronary_tuning_and_simulation'
  print('\nMoving simulation files into: ' + rigid_sim_folder + '\n')
  command_string = 'mkdir ' + rigid_sim_folder
  print(command_string)
  os.system(command_string)
  command_string = 'mv restart.0.1 geombc.dat.1 numstart.dat solver.inp coronaryModel.txt coronaryParams.txt coronary.csv ' + rigid_sim_folder
  print(command_string)
  os.system(command_string)
  
#-------------------------------------------------------------------------------

def deformablePresolve():
  input_path = MESH_SURFACES_PATH
  if(input_path[-1:] is not '/'):
    input_path = input_path + '/'
  if(input_path[-1:] is not '*'):
    input_path = input_path + '*'
  Tc = 60.0/heart_rate

  print("This script will post process a SimVascular closed-loop coronary simulation")
  print("\nWARNING: You must be in the simulation directory for this script to work, ")
  print("and you must have properly set the paths to the SimVascular pre-solver and")
  print("post-solver before continuing")

  # First, need to get the name of the simulation folder
  procs_list = glob.glob('*-procs_*')
  if(len(procs_list) == 0):
    print("ERROR: Could not find n-procs_case simulation folder. Pleasre ensure you")
    print(" are inside your simulation folder")
    sys.exit()
  sim_folder = procs_list[0]

  # Read through the solver.inp file to get the number of timesteps in this
  # simulation
  solver_file = open('solver.inp', 'r')
  for line in solver_file:
    line_split = line.split()
    if(len(line_split) > 2):
      if(line_split[2] == 'Timesteps:'):
        total_steps = int(line_split[3])
    if(len(line_split) > 4):
      if(line_split[4] == 'Restarts:'):
        save_every = int(line_split[5])
  solver_file.close()

  # Now go into the n-procs_case folder and run the post-solver
  step_to_post = total_steps - save_every
  command_string = "cd " + sim_folder
  print(command_string)
  os.chdir(sim_folder)
  print("\n ** Post-processing rigid simulation results...\n")
  command_string = POST_SOLVER_PATH + ' -ph -sn ' + str(step_to_post)
  print(command_string)
  os.system(command_string)

  restart_file = 'restart.' + str(step_to_post) + '.0'
  command_string = 'cp ' + restart_file + ' ..'
  print(command_string)
  os.system(command_string)
  command_string = 'cd ..'
  print(command_string)
  os.chdir('..')

  # Move rigid results into a new directory
  rigid_results_dir = 'tuned_rigid_results'
  print('\nMoving tuned rigid results into ' + rigid_results_dir)
  command_string = 'mkdir ' + rigid_results_dir
  print(command_string)
  os.system(command_string)
  command_string = 'mv restart.0.1 geombc.dat.1 solver.inp ' + sim_folder + ' ' + rigid_results_dir
  print(command_string)
  os.system(command_string)
  command_string = 'cp GenBC numstart.dat ' + rigid_results_dir
  print(command_string)
  os.system(command_string)

  # Rename the post-processed restart file to prepare for FSI pre-solver
  command_string = 'mv ' + restart_file + ' restart.0.1'
  print(command_string)
  os.system(command_string)

  # Read the mesh-surfaces folder to get all the face names and areas
  # Areas are needed if applying variable wall properties
  filelist_raw = glob.glob(input_path)
  filelist_raw.sort()

  filelist = []
  for trial in filelist_raw:
    if(trial[-4:] == ".vtp"):
      filelist.append(trial)
  # Check to make sure we account for all surfaces
  num_surfaces = len(filelist)

  base_len = len(input_path)
  A_br = []
  A_RCA = []
  A_LCA = []
  A_LCX = []
  chosen_positions = []
  chosen_names = []
  inlet_position = -1
  inflow_name = ''

  # Need to save wall surfaces too, for variable wall pre-solver
  aorWallNames = []
  corWallNames = []
  graftWallNames = []
  limaWallNames = []
  inflow_area = 0.0
  aorWallTag = 'wall_'
  lcaWallTag = 'wall_'
  lcxWallTag = 'wall_'
  rcaWallTag = 'wall_'
  graftWallTag = 'wall_svg'
  limaWallTag = 'wall_lima'
  blendWallTag = 'wall_blend'
  if(aorBranchTag[:4] == 'cap_'):
    aorWallTag = aorWallTag + aorBranchTag[4:]
  else:
    aorWallTag = aorWallTag + aorBranchTag
  if(lcaBranchTag[:4] == 'cap_'):
    lcaWallTag = lcaWallTag + lcaBranchTag[4:]
  else:
    lcaWallTag = lcaWallTag + lcaBranchTag
  if(lcxBranchTag[:4] == 'cap_'):
    lcxWallTag = lcxWallTag + lcxBranchTag[4:]
  else:
    lcxWallTag = lcxWallTag + lcxBranchTag
  if(rcaBranchTag[:4] == 'cap_'):
    rcaWallTag = rcaWallTag + rcaBranchTag[4:]
  else:
    rcaWallTag = rcaWallTag + rcaBranchTag

  # Sort through the input file list, and save the different face areas into
  # the appropriate container depending on the face name
  surf_count = 0
  for vtp_file in filelist:
    temp_name = vtp_file[base_len-1:]
    if(temp_name[:len(aorBranchTag)] == aorBranchTag):
      A_br.append(findVTPArea(vtp_file))
      chosen_positions.append(findVTPMassCenter(vtp_file))
      chosen_names.append(vtp_file)
      surf_count = surf_count + 1
    elif(temp_name[:len(lcaBranchTag)] == lcaBranchTag):
      A_LCA.append(findVTPArea(vtp_file))
      chosen_positions.append(findVTPMassCenter(vtp_file))
      chosen_names.append(vtp_file)
      surf_count = surf_count + 1
    elif(temp_name[:len(lcxBranchTag)] == lcxBranchTag):
      A_LCX.append(findVTPArea(vtp_file))
      chosen_positions.append(findVTPMassCenter(vtp_file))
      chosen_names.append(vtp_file)
      surf_count = surf_count + 1
    elif(temp_name[:len(rcaBranchTag)] == rcaBranchTag):
      A_RCA.append(findVTPArea(vtp_file))
      chosen_positions.append(findVTPMassCenter(vtp_file))
      chosen_names.append(vtp_file)
      surf_count = surf_count + 1
    elif(temp_name[:len(inflowTag)] == inflowTag):
      inflow_area = findVTPArea(vtp_file)
      inlet_position = findVTPMassCenter(vtp_file)
      inflow_name = vtp_file
      surf_count = surf_count + 1
    elif(temp_name[:len(aorWallTag)] == aorWallTag):
      aorWallNames.append(vtp_file)
      surf_count = surf_count + 1
    elif(temp_name[:len(lcaWallTag)] == lcaWallTag):
      corWallNames.append(vtp_file)
      surf_count = surf_count + 1
    elif(temp_name[:len(lcxWallTag)] == lcxWallTag):
      corWallNames.append(vtp_file)
      surf_count = surf_count + 1
    elif(temp_name[:len(rcaWallTag)] == rcaWallTag):
      corWallNames.append(vtp_file)
      surf_count = surf_count + 1
    elif(temp_name[:len(graftWallTag)] == graftWallTag):
      graftWallNames.append(vtp_file)
      surf_count = surf_count + 1
    elif(temp_name[:len(limaWallTag)] == limaWallTag):
      limaWallNames.append(vtp_file)
      surf_count = surf_count + 1
    elif(temp_name[:len(blendWallTag)] == blendWallTag):
      surf_count = surf_count + 1

  # If we do not account for all surfaces, exit and notify the user
  if(not(surf_count == num_surfaces)):
    print('ERROR: Not all surfaces were found by C2. Please ensure your surface names adhere to the naming convention')
#    exit()
      

  # Append all lcx faces to the end of the lca area list
  for lcx_area in A_LCX:
    A_LCA.append(lcx_area)

  # Assign the number of outlets, number of unknowns, types of outlets, and faceToStateMapping
  nCOR_l = len(A_LCA)
  nCOR_r = len(A_RCA)
  nCOR = nCOR_l + nCOR_r
  nRCR = len(A_br)
  nFaces = nCOR + nRCR + 1

  # Now need to make the pre-solver file
  svpre_file = open('coronaryLPN_FSI.svpre','w')
  mesh_complete = input_path[:-1] + '../' + 'mesh-complete.mesh.vtu'
  walls_combined = input_path[:-1] + '../' + 'walls_combined.vtp'
  mesh_exterior = input_path[:-1] + '../' + 'mesh-complete.exterior.vtp'
  zp_string = 'zero_pressure_vtp '
  id_string = 'set_surface_id_vtp '

  svpre_file.write('mesh_and_adjncy_vtu ' + mesh_complete + '\n')
  svpre_file.write('deformable_wall_vtp ' + walls_combined +'\n')
  svpre_file.write('fix_free_edge_nodes_vtp ' + walls_combined + '\n\n')

  for name in chosen_names:
    temp_string = zp_string + name + '\n'
    svpre_file.write(temp_string)

  svpre_file.write(zp_string + inflow_name + '\n\n')
  temp_string = id_string + mesh_exterior + ' 1\n'
  svpre_file.write(temp_string)
  count = 2

  for name in chosen_names:
    temp_string = id_string + name + ' ' + str(count) + '\n'
    svpre_file.write(temp_string)
    count = count + 1

  svpre_file.write(id_string + inflow_name + ' ' + str(count) + '\n\n')

  if(UNIFORM_WALL_PROPERTIES):
    svpre_file.write('deformable_E ' + str(deformable_Evw) + '\n')
    svpre_file.write('deformable_nu ' + str(deformable_nuvw) + '\n')
    svpre_file.write('deformable_thickness ' + str(deformable_thickness) + '\n')
    svpre_file.write('deformable_kcons ' + str(deformable_kcons) + '\n') 
    svpre_file.write('write_geombc geombc.dat.1 \n')
    svpre_file.write('deformable_pressure ' + str(deformable_pressure) + '\n')
    svpre_file.write('deformable_solve_displacements\n') 
    svpre_file.write('append_displacements restart.0.1\n')
  else: # TODO: Variable wall properties
    var_thick = 'set_surface_thickness_vtp '
    var_E = 'set_surface_E_vtp '

    # First, specify the thicknesses
    rikka = 0
    thick_temp = 0.0
    for aorOut in A_br:
      thick_temp = 0.1*np.sqrt(aorOut/math.pi)
      write_string = var_thick + chosen_names[rikka] + ' ' + str(thick_temp) + '\n'
      svpre_file.write(write_string)
      rikka = rikka + 1

    param_A = 0.000387
    param_B = 0.63
    for lcaOut in A_LCA:
      thick_temp = param_A*((np.sqrt(lcaOut/math.pi)*10000.0)**(param_B))
      write_string = var_thick + chosen_names[rikka] + ' ' + str(thick_temp) + '\n'
      svpre_file.write(write_string)
      rikka = rikka + 1
    for rcaOut in A_RCA:
      thick_temp = param_A*((np.sqrt(rcaOut/math.pi)*10000.0)**(param_B))
      write_string = var_thick + chosen_names[rikka] + ' ' + str(thick_temp) + '\n'
      svpre_file.write(write_string)
      rikka = rikka + 1
    thick_temp = 0.1*np.sqrt(inflow_area/math.pi)
    write_string = var_thick + inflow_name + ' ' + str(thick_temp) + '\n'
    svpre_file.write(write_string)
    svpre_file.write('solve_varwall_thickness\n\n')

    # Now, set commands to specify elastic modulus on the walls
    for aorWall in aorWallNames:
      write_string = var_E + aorWall + ' ' + str(deformable_aortic_E) + '\n'
      svpre_file.write(write_string)
    for corWall in corWallNames:
      write_string = var_E + corWall + ' ' + str(deformable_coronary_E) + '\n'
      svpre_file.write(write_string)
    for graftWall in graftWallNames:
      write_string = var_E + graftWall + ' ' + str(deformable_graft_E) + '\n'
      svpre_file.write(write_string)
    for limaWall in limaWallNames:
      write_string = var_E + limaWall + ' ' + str(deformable_lima_E) + '\n'
      svpre_file.write(write_string)
    svpre_file.write('solve_varwall_E\n\n')

    svpre_file.write('varwallprop_write_vtp\n')
    svpre_file.write('deformable_nu ' + str(deformable_nuvw) + '\n')
    svpre_file.write('deformable_kcons ' + str(deformable_kcons) + '\n') 
    svpre_file.write('deformable_pressure ' + str(deformable_pressure) + '\n')
    svpre_file.write('deformable_solve_varwall_displacements\n')
    svpre_file.write('wall_displacements_write_vtp\n\n')
    svpre_file.write('write_geombc geombc.dat.1 \n\n')
    svpre_file.write('append_displacements restart.0.1\n')

  svpre_file.close()

  # Generate a solver.inp that can be used for closed-loop simulations
  solver_writer = open('solver.inp','w')

  solver_writer.write('# Phasta Version 1.5 Input File\n')
  solver_writer.write('# Produced by generateCoronaryLPN\n')
  solver_writer.write('# Specialized for closed-loop coronary simulations\n\n')

  solver_writer.write('# ****** SOLUTION CONTROL ****** #\n')
  solver_writer.write('Equation of State: Incompressible\n')
  solver_writer.write('Number of Timesteps: %d\n'%(N_timesteps*N_cycles))
  solver_writer.write('Time Step Size: {0:.8f}\n\n'.format(Tc/float(N_timesteps)))

  solver_writer.write('# ***** FLUID PROPERTIES ***** #\n')
  solver_writer.write('Viscosity: 0.04\n')
  solver_writer.write('Density: 1.06\n\n')

  solver_writer.write('# ***** OUTPUT CONTROL ***** #\n')
  solver_writer.write('Number of Timesteps between Restarts: %d\n'%N_timesteps)
  solver_writer.write('Print ybar: True\n\n')

  solver_writer.write('# ****** CARDIOVASCULAR MODELING PARAMETERS ***** #\n')
  solver_writer.write('Time Varying Boundary Conditions From File: False\n')
  solver_writer.write('Number of Coupled Surfaces: {0}\n'.format(nFaces))
  solver_writer.write('Pressure Coupling: Implicit\n')
  solver_writer.write('Number of Resistance Surfaces: 0\n')
  solver_writer.write('Find the GenBC Inside the Running Directory: True\n')
  solver_writer.write('Number of Timesteps for GenBC Initialization: 0\n')
  solver_writer.write('Number of Neumann Surfaces: {0}\n'.format(nFaces))
  neumann_list = 'List of Neumann Surfaces: '
  for i in range(2, nFaces+2, 1):
    neumann_list = neumann_list + str(i) + ' '
  solver_writer.write(neumann_list + '\n\n')

  solver_writer.write('# ***** LINEAR SOLVER (SHOULD LEAVE THIS ALONE) ****** #\n')
  solver_writer.write('Number of Solves per Left-hand-side Formation: 1\n\n')
  solver_writer.write('# NOTE: Recommended to change this to GMRES for Deformable simulations #\n')
  solver_writer.write('svLS Type: GMRES\n')
  solver_writer.write('Number of Krylov Vectors per GMRES Sweep: 200\n')
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
  solver_writer.write('Minimum Required Iterations : 3\n')
  solver_writer.write('Step Construction : 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1\n\n')

  solver_writer.write('Backflow Stabilization Coefficient: 0.25\n\n')

  solver_writer.write('# ***** DEFORMABLE SIMULATION SETTINGS ***** #\n')
  solver_writer.write('Deformable Wall: True\n')
  if(not UNIFORM_WALL_PROPERTIES):
    solver_writer.write('Variable Wall Thickness and Young Mod: True\n')
  else:
    solver_writer.write('Variable Wall Thickness and Young Mod: False\n')
    solver_writer.write('Thickness of Vessel Wall: ' + str(deformable_thickness) + '\n')
    solver_writer.write('Young Mod of Vessel Wall: ' + str(deformable_Evw) + '\n')
  solver_writer.write('Density of Vessel Wall: ' + str(deformable_density) + '\n')
  solver_writer.write('Number of Wall Properties per Node: 10\n')
  solver_writer.write('Poisson Ratio of Vessel Wall: {0:.6f}\n'.format(deformable_nuvw))
  solver_writer.write('Shear Constant of Vessel Wall: {0:.6f}\n'.format(deformable_kcons))
  solver_writer.write('Wall Mass Matrix for LHS: True\n')

  solver_writer.close()

  # Run the pre-solver to generate FSI simulation files
  command_string = PRESOLVER_PATH + ' coronaryLPN_FSI.svpre'
  print(command_string)
  os.system(command_string)

  # Update the coronaryModel.txt to include an inlet capacitor for FSI tuning
  modelFile = open('coronaryModel.txt','a')
  modelFile.write('Lin,4.0\n') 
  modelFile.write('Csurr,0.0006\n')
  modelFile.close()

  # Update coronaryParams.txt to lower the inlet capacitance Cpa
  print('\n** Adjusting LPN inlet compliance to account for FSI effects...\n')
  command_string = 'mv coronaryParams.txt coronaryParams_temp.txt'
  print(command_string)
  os.system(command_string)
  base_file = open('coronaryParams_temp.txt', 'r')
  model_file = open('coronaryParams.txt', 'w')
  count = 1
  for line in base_file:
    if(count == 18):
      new_Cpa = 0.25 * float(line)
      model_file.write(str(new_Cpa) + '\n')
    else:
      model_file.write(line)
    count = count + 1
  base_file.close()
  model_file.close()
  command_string = 'rm coronaryParams_temp.txt'
  print(command_string)
  os.system(command_string)

  # Move all relevant simulation files into a seperate directory
  FSI_sim_folder = 'FSI_coronary_tuning_and_simulation'
  print('\nMoving simulation files into: ' + FSI_sim_folder + '\n')
  command_string = 'mkdir ' + FSI_sim_folder
  print(command_string)
  os.system(command_string)
  command_string = 'cp restart.0.1 geombc.dat.1 numstart.dat solver.inp coronaryModel.txt coronaryParams.txt coronary.csv ' + FSI_sim_folder
  print(command_string)
  os.system(command_string)

#-------------------------------------------------------------------------------

if __name__ == "__main__":
  for Counter in range(8): #Loop to make sure left-right coronary split is maintained
    print ("-"*200)
    print('This script is for automatically tuning the boundary conditions for closed loop coronary simulations')
    print('\n First, we use the SimVascular pre-solver to process the mesh')
    rigidPresolve()
  
    # At this point, all the files we need for rigid tuning are located in a folder
    # called 'rigid_coronary_tuning_and_simulation' which is located in the current folder
    # we use system tools to move there and start the rigid tuning
    command_string = 'cd rigid_coronary_tuning_and_simulation'
    print(command_string)
    os.chdir('rigid_coronary_tuning_and_simulation')

    # Copy over necesary scripts and executables from gcode
    command_string = 'cp ' + GCODE_BINARY_DIR + 'bin/GenBC .'
    print(command_string)
    os.system(command_string)
    command_string = 'cp ' + GCODE_BINARY_DIR + 'py/optimizeSurrogateRes.py .'
    print(command_string)
    os.system(command_string)

    print('Making batch scripts for NM_optimization, SimVascular flowsolver, and surrogate resistance optimization...')
    makeNM_script()
    makeFlowsolver_script()
    makeSurrResOpt_script()
    coronaryTuningIteration(1, True)
  
    # Once rigid tuning is done, we launch into the deformable presolver
    deformablePresolve()
  
    # Now that deformable pre-solving is done, the files required for FSI_tuning
    # are located in a folder called 'FSI_coronary_tuning_and_simulation'
    # Move this folder one up then run the FSI tuning from there
    command_string = 'mv FSI_coronary_tuning_and_simulation ..'
    print(command_string)
    os.system(command_string)
  
    # Now we move into that directory
    command_string = 'cd ../FSI_coronary_tuning_and_simulation' #add ../ at the front
    print(command_string)
    os.chdir('../FSI_coronary_tuning_and_simulation')
  
    # Copy over necesary scripts and executables from gcode
    command_string = 'cp ' + GCODE_BINARY_DIR + 'bin/GenBC .'
    print(command_string)
    os.system(command_string)
    command_string = 'cp ' + GCODE_BINARY_DIR + 'py/optimizeSurrogateCom.py .'
    print(command_string)
    os.system(command_string)

    print('Making batch scripts for NM_optimization, SimVascular flowsolver, and surrogate resistance optimization...')
    makeNM_script()
    makeFlowsolver_script()
    makeSurrComOpt_script()
    coronaryTuningIteration(1, False)

    #Change the working directory to original one
    os.chdir(WORKING_DIRECTORY)

    #Save the Flow Splits to the file
    with open("LeftRightFlowSplit.dat",'a') as Flow_split_file:
      Flow_split_file.write("%d %.04f %.04f %.04f %.04f %.04f %.04f\n"%(Counter, sys_cor_split, Aor_Cor_split_computed, Left_Cor_split_assigned, Left_Cor_split_computed, f_l, f_r))  

   #Check to ensure left and right coronary split is maintained
    Diff=(Left_Cor_split_computed-Left_Cor_split_assigned)
    print ("The Left Coronary Split is: %.02f "%Left_Cor_split_computed)
    if Diff<2 and Diff>-2:
      break
    elif Diff>2:
      f_l=f_l*(Left_Cor_split_assigned/Left_Cor_split_computed)
      f_r=(sys_cor_split/100.)-f_l
      print ("The Updated Left Split is: %.03f"%f_l)
      print ("The Updated Right Split is: %.03f"%f_r)
      os.system("mkdir Coronary_Flow_Split_Iteration%d"%Counter)
      os.system("mv FSI_coronary_tuning_and_simulation Coronary_Flow_Split_Iteration%d/"%Counter)
      os.system("mv rigid_coronary_tuning_and_simulation Coronary_Flow_Split_Iteration%d/"%Counter)
      os.system("mv coronaryLPN.svpre Coronary_Flow_Split_Iteration%d/"%Counter)

    elif Diff<-2:
      f_l=f_l*(Left_Cor_split_computed/Left_Cor_split_assigned)
      f_r=(sys_cor_split/100.)-f_l
      print ("The Updated Left Split is: %.03f"%f_l)
      print ("The Updated Right Split is: %.03f"%f_r)
      os.system("mkdir Coronary_Flow_Split_Iteration%d"%Counter)
      os.system("mv FSI_coronary_tuning_and_simulation Coronary_Flow_Split_Iteration%d/"%Counter)
      os.system("mv rigid_coronary_tuning_and_simulation Coronary_Flow_Split_Iteration%d/"%Counter)
      os.system("mv coronaryLPN.svpre Coronary_Flow_Split_Iteration%d/"%Counter)

    else:
      print ("There is an error. Check the Flow Split Conditions")


    
#-------------------------------------------------------------------------------





