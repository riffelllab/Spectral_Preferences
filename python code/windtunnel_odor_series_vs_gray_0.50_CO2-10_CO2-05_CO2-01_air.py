"""
- WARNING THIS SCRIPT USES PYTHON 3 

Script to control both an Alicat mass flow controller and an arduino controled LED stimuli via python.

------ Mass Controller Info --------
Alicat have a package (https://github.com/numat/alicat), you will need to install
it using from source or using pip (pip install alicat)

If the flow controller is communicating with the python scrip, it will produce a dictionary as this:
{
  'setpoint': 0.0,         # Setpoint, either mass flow rate or pressure (values between  0 to 100%)
  'control_point': 'flow', # Either 'flow' or 'pressure'
  'gas': 'Air',            # Can be any option in `flow_controller.gases`
  'mass_flow': 0.0,        # Mass flow (in units specified at time of purchase)
  'pressure': 25.46,       # Pressure (normally in psia)
  'temperature': 23.62,    # Temperature (normally in C)
  'total_flow': 0.0,       # Optional. If totalizer function purchased, will be included
  'volumetric_flow': 0.0   # Volumetric flow (in units specified at time of purchase)
}

From the Alicat mass flow controller documentation. When you use FlowController.get() 
  the values obtained from the MAss Flow Controller are:
    * Pressure (normally in psia)
    * Temperature (normally in C)
    * Volumetric flow (in units specified at time of order)
    * Mass flow (in units specified at time of order)
    * Flow setpoint (in units of control point, values between  0 to 100% )
    * Flow control point (either 'flow' or 'pressure')
    * Total flow (only on models with the optional totalizer function)
    * Currently selected gas


WARNING: If user cannot operate with the Mass Flow Controller (Permission denied: '/dev/ttyUSB0'), 
        add the user to the usergroup dialout ($ sudo usermod -a -G dialout riffell), logout and log in again.


Before use, you will need to check the following in the Mass Flow Controller:
    * The volumetric flow units are set to SCCM (Standard Cubic Centimeter per Minute)
    * The communication is set to Serial
 (https://documents.alicat.com/manuals/Gas_Flow_Controller_Manual.pdf)


------ Arduino LED stimuli Info --------
simple program to compile and upload Arduino sketches via the command line

"""

#Loads libraries to start and stop recording for the camera
import os

#Mass controller libraries
from alicat import FlowController
from serial import Serial
import time

#Arduino libraries
import sys
import subprocess


#Libraries to output timestamps to file
from datetime import datetime
import csv
import glob

#Import functions for Google Drive Upload
from pydrive.auth import GoogleAuth
from pydrive.drive import GoogleDrive
#Specify the folder id for the google drive folder where you'd like
#the output folders uploaded. The ID can be found at the end of the google 
#drive folder URL, IE https://drive.google.com/drive/folders/1hq8xCobkucQIn_ibj3ACVcQjjwsg_2Mh
Upload_Folder_ID = '19nSXlc_PLVaWDYfPZkFh2Un49ltYZ1Kj'


#------ Mass Controller Settings -----

# Connect to the flow controllers
#***these ports identified by their usb plug because the devices are identical and have no serial number
#   if the usb connections are unplugged and not returned to the same possition the improper device
#   may be referenced below
serialPort1='/dev/serial/by-path/pci-0000:00:14.0-usb-0:2:1.0-port0'      # used to control air bypassing odor
serialPort2='/dev/serial/by-path/pci-0000:00:14.0-usb-0:4:1.0-port0'      # used to control CO2 to odor
# ttyUSB3 used to control air to odor, identifed by id and uaffected by plug position
serialPort3='/dev/serial/by-id/usb-Prolific_Technology_Inc._USB-Serial_Controller_FDCKs11CN11-if00-port0'


#------ Arduino Settings -----
arduinoProg = "/home/wtunnel/Downloads/arduino-1.8.10/arduino" #location of arduino program
projectFileDir = "/home/wtunnel/Desktop/Arduino Code/color, vs gray_0.50/" #folder location for arduino sketches
boardLine = "arduino:avr:uno" #type of board (see arduino IDE)
portLine = "/dev/ttyACM0" #port used for board (see arduino IDE)

#Sets the arduino program to use for the neutral and off states
arduino_neutral = 'a_gray_1.00_b_gray_1.00'
arduino_off = 'a_all_0.00_b_all_0.00'


#---Defines the series of arduino stimuli to upload
   
gray_VGR_black = [
   'a_gray_1.00_b_gray_1.00',
   'a_gray_1.00_b_gray_1.00',
   'a_430_1.00_b_gray_0.50',
   'a_gray_0.50_b_430_1.00',
   'a_525_1.00_b_gray_0.50',
   'a_gray_0.50_b_525_1.00',
   'a_625_1.00_b_gray_0.50',
   'a_gray_0.50_b_625_1.00',
   'a_all_0.00_b_gray_0.50',
   'a_gray_0.50_b_all_0.00'
]

black_RGV_gray = [
   'a_all_0.00_b_gray_0.50',
   'a_gray_0.50_b_all_0.00',
   'a_625_1.00_b_gray_0.50',
   'a_gray_0.50_b_625_1.00',
   'a_525_1.00_b_gray_0.50',
   'a_gray_0.50_b_525_1.00',
   'a_430_1.00_b_gray_0.50',
   'a_gray_0.50_b_430_1.00',
   'a_gray_1.00_b_gray_1.00',
   'a_gray_1.00_b_gray_1.00'
]


#------ Time Settings -----
uploadTime = 7		# Estimate of the time in s required to upload an arduino program
waitTime = 900		# Wait 30 min (1800 s) before starting air
measureBuffer = 5	# Small wait inbetween mass flow controller setting upload and measurement
odorBuffer = 300	# 5 minute buffer period to allow previous odor to clear
preCO2Time = 600  	# Record behavioral for 10 min (600 s) with no CO2
StimTime = 300        	# length of each LED color stimulus 5 min (300 s)
interStimTime = 30	# length of time between each set of stimuli, 30
postCO2Time = 600	# Record behavioral for 10 min (600 s) with no CO2

# Total duration (in seconds) of the experiments
totalExpDuration = \
  waitTime + \
  3 * measureBuffer + \
  3 * odorBuffer + \
  preCO2Time + \
  4 * len(black_RGV_gray) * (StimTime + interStimTime) 


#------ Time Stamp Output Settings -----
ts_output_folder='/home/wtunnel/FLYDRA/'
header = ['event','date time', 'timecode']
time_stamps = [header]


#------ Mass Controller Functions -----

#Function to change the mass flow control setpoint in SCCM 
def change_setpoint(flowCtrller, value):
  flowCtrller._set_setpoint(value)

#Function to set the gas type
def set_gas(flowCtrller, gas):
  flowCtrller.set_gas(gas)


#------ Arduino Functions -----

def arduino_upload(project):
	#Concatenates the full upload command for below
	arduinoCommand = arduinoProg + " --" + "upload" + " --board " + boardLine + " --port " + portLine + " " \
		+ "\'" + projectFileDir + project + "/" + project + ".ino" + "\'"
	#Print message to terminal indicating the start of the upload
	print ("\n-- Starting Upload of %s --" %(project))
	#Upload command
	presult = subprocess.call(arduinoCommand, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
	#Print a sucess or failure message to terminal
	if presult != 0:
		print ("-- Upload Failed --")
	else:
		print ("-- Success --")
  

#----- time stamp output functions -----

def add_timestamp(event):
    """
    Function to add a timestamps to a running list
    """
    #Gets timestamp when function called
    timestamp = time.time()
    #Convert timestamp to date and time
    dt = datetime.fromtimestamp(timestamp)
    #Add the "event" from the function call, and the time variables to a row list
    row = [event, dt.strftime("%m/%d/%Y, %H:%M:%S"), timestamp]
    #Add the "row" to the time stamp list
    time_stamps.append(row)

def save_timestamps():  
    """
    Function to save timestampsas a csv file to the FLYDRA folder with H5 file
    """
    try:
      #Gets the name of the most recent h5 file
      h5_path = max(glob.iglob( ts_output_folder + "*.mainbrain.h5"), key=os.path.getctime)
      #replaces the file extension with ".csv"
      csv_path = h5_path.replace(".mainbrain.h5", ".csv")
      #writes the csv from the saved timestamps
      with open(csv_path, 'w') as csv_output:
        writer = csv.writer(csv_output)
        writer.writerows(time_stamps)
      print('\nTimestamps sucsessfully saved in %s'%(csv_path))
    except Exception as e:
      print('\n ========= ')
      print('\n* ERROR: Data not saved')
      print("* Exception msg: %s"%e)
    try:
      #set google drive authorization
      gauth = GoogleAuth()           
      drive = GoogleDrive(gauth)
      #Create upload folder object for use below
      Upload_Folder = drive.CreateFile({'id': Upload_Folder_ID})
      #creates a file on google drive to receive the upload
      csv_file = drive.CreateFile({
        'title': os.path.basename(csv_path), 
        #Sets the parent directory on google drive
        'parents': [{'id': Upload_Folder_ID}]
        })
      #Sets local file path for upload
      csv_file.SetContentFile(csv_path)
      #Uploads the file
      csv_file.Upload()
      #creates a file on google drive to receive the upload
      h5_file = drive.CreateFile({
        'title': os.path.basename(h5_path), 
        #Sets the parent directory on google drive
        'parents': [{'id': Upload_Folder_ID}]
        })
      #Sets local file path for upload
      h5_file.SetContentFile(h5_path)
      #Uploads the file
      h5_file.Upload()
      #Display folder title and create link to the folder
      print('Uploaded to google drive in the folder \'%s\': \n\
      https://drive.google.com/drive/folders/%s' % 
      (Upload_Folder['title'], Upload_Folder['id']))
    except Exception as e:
      print('\n ========= ')
      print('\n* ERROR: Data not uploaded to Google Drive')
      print("* Exception msg: %s"%e)



try:
  CleanAirCtr = FlowController(port= serialPort1)
  CO2_Ctr = 	FlowController(port= serialPort2)
  OdorAirCtr = 	FlowController(port= serialPort3)

  #Set the gas for each flow controller
  set_gas(CleanAirCtr, 'Air')
  set_gas(CO2_Ctr, 'CO2')
  set_gas(OdorAirCtr, 'Air')

  #Ensure that the flow rate of both mass controllers is 0
  change_setpoint(CleanAirCtr, 0.0)
  change_setpoint(CO2_Ctr, 0.0)
  change_setpoint(OdorAirCtr, 0.0)

  #Wait 5 seconds to allow mass flow controllers to bring the flow to the set point
  time.sleep(measureBuffer) 

  #Get current setting for output to console
  Settings_CleanAirCtr = CleanAirCtr.get()
  Settings_CO2_Ctr = 	 CO2_Ctr.get()
  Settings_OdorAirCtr =  OdorAirCtr.get()

  #Add timestamp for start of experiment
  add_timestamp('start_WaitTime_' + arduino_neutral)

  #Starts recording from cameras
  os.system("python ~/ros/flydra-kinetic/src/ros_flydra/scripts/simple-console.py --start-saving-hdf5")

  print('\n ================== ')
  print('\nStarting experiment at: %s' %time.ctime(time_stamps[-1][2]+waitTime))
  print('Experiment will end at: %s' %time.ctime(time_stamps[-1][2]+totalExpDuration))
  print('\n ================== ')

  print('\nStarting initial wait period, airflow off')
  print('\nCurrSett MassFlowController 1 %s'%Settings_CleanAirCtr)
  print('\nCurrSett MassFlowController 2 %s'%Settings_CO2_Ctr)
  print('\nCurrSett MassFlowController 3 %s'%Settings_OdorAirCtr)

  #LEDs controller set to neutral state
  print('\nPutting Arduino LED controller in neutral state (all LEDs off)')
  arduino_upload(arduino_neutral)

  #Wait until the insect get used to the test section
  time.sleep(waitTime-uploadTime)

  #--------------------------------------------------------------------------------

  #Start the flow of air, CO2, and odor
  change_setpoint(CleanAirCtr, 190.0)
  change_setpoint(CO2_Ctr, 0)
  change_setpoint(OdorAirCtr, 0)

  #Wait 5 seconds to allow mass flow controllers to bring the flow to the set point
  time.sleep(measureBuffer) 

  #Get current setting for output to console
  Settings_CleanAirCtr = CleanAirCtr.get()
  Settings_CO2_Ctr = 	 CO2_Ctr.get()
  Settings_OdorAirCtr =  OdorAirCtr.get()

  #Add timestamp for start of experiment
  add_timestamp('start_PreCO2Time_' + arduino_neutral)
  print('\n ================== ')
  print('\nStarting pre-CO2 period, clean air released at %s'%time_stamps[-1][1])
  print('\nCurrSett clean air controller %s'%Settings_CleanAirCtr)
  print('\nCurrSett CO2 controller %s'%Settings_CO2_Ctr)
  print('\nCurrSett odor air controller %s'%Settings_OdorAirCtr)

  #LEDs controller set to neutral state
  print('\nPutting Arduino LED controller in neutral state')
  arduino_upload(arduino_neutral)

  #Wait pre-C02 period
  time.sleep(preCO2Time-uploadTime)

  #--------------------------------------------------------------------------------

  #Start CO2 and sets the flow of CO2 to 10% of total flow
  change_setpoint(CleanAirCtr, 190.0)
  change_setpoint(CO2_Ctr, 21.5)
  change_setpoint(OdorAirCtr, 0.0)

  #Wait 5 minutes to clear out the previous odor
  time.sleep(odorBuffer) 

  #Get current setting for output to console
  Settings_CleanAirCtr = CleanAirCtr.get()
  Settings_CO2_Ctr = 	 CO2_Ctr.get()
  Settings_OdorAirCtr =  OdorAirCtr.get()

  print('\n ================== ')
  print('\nStarting 10%% CO2 period at %s'%time_stamps[-1][1])
  print('\nCurrSett clean air controller %s'%Settings_CleanAirCtr)
  print('\nCurrSett CO2 controller %s'%Settings_CO2_Ctr)
  print('\nCurrSett odor air controller %s'%Settings_OdorAirCtr)

  print('\n ========= ')
  print('\n LED stimuli starting')
  print('\n ========= ')

  #Displays a range of test stimuli
  for i in range(len(gray_VGR_black)):
    #Add timestamp for stimulus
    add_timestamp("start_" + "CO2 10_" + gray_VGR_black[i] )
    #Upload stimulus
    arduino_upload(gray_VGR_black[i])
    #Leave LED stimulus on for the StimTime
    time.sleep(StimTime-uploadTime)
    #Add timestamp for stimulus
    add_timestamp("end_" + "CO2 10_" + gray_VGR_black[i])
    #outputs time stamp to console
    print ("-- Stimulus period complete at %s --"%time_stamps[-1][1])
    #LEDs controller set to neutral state
    arduino_upload(arduino_neutral)
    #Leave LEDs in neutral state for the interStimTime
    time.sleep(interStimTime-uploadTime)

#--------------------------------------------------------------------------------

  #Set the flow of CO2 to 5% of total flow
  change_setpoint(CleanAirCtr, 211.0)
  change_setpoint(CO2_Ctr, 10.5)
  change_setpoint(OdorAirCtr, 0.0)

  #Wait 5 minutes to clear out the previous odor
  time.sleep(odorBuffer) 

  #Get current setting for output to console
  Settings_CleanAirCtr = CleanAirCtr.get()
  Settings_CO2_Ctr = 	 CO2_Ctr.get()
  Settings_OdorAirCtr =  OdorAirCtr.get()

  print('\n ================== ')
  print('\nStarting 5%% CO2 period at %s'%time_stamps[-1][1])
  print('\nCurrSett clean air controller %s'%Settings_CleanAirCtr)
  print('\nCurrSett CO2 controller %s'%Settings_CO2_Ctr)
  print('\nCurrSett odor air controller %s'%Settings_OdorAirCtr)

  print('\n ========= ')
  print('\n LED stimuli starting')
  print('\n ========= ')

  #Displays a range of test stimuli
  for i in range(len(gray_VGR_black)):
    #Add timestamp for stimulus
    add_timestamp("start_" + "CO2 05_" + gray_VGR_black[i] )
    #Upload stimulus
    arduino_upload(gray_VGR_black[i])
    #Leave LED stimulus on for the StimTime
    time.sleep(StimTime-uploadTime)
    #Add timestamp for stimulus
    add_timestamp("end_" + "CO2 05_" + gray_VGR_black[i])
    #outputs time stamp to console
    print ("-- Stimulus period complete at %s --"%time_stamps[-1][1])
    #LEDs controller set to neutral state
    arduino_upload(arduino_neutral)
    #Leave LEDs in neutral state for the interStimTime
    time.sleep(interStimTime-uploadTime)
    
  #--------------------------------------------------------------------------------

  #Set the flow of CO2 to 1% of total flow
  change_setpoint(CleanAirCtr, 219.3)
  change_setpoint(CO2_Ctr, 2.2)
  change_setpoint(OdorAirCtr, 0.0)

  #Wait 5 minutes to clear out the previous odor
  time.sleep(odorBuffer) 

  #Get current setting for output to console
  Settings_CleanAirCtr = CleanAirCtr.get()
  Settings_CO2_Ctr = 	 CO2_Ctr.get()
  Settings_OdorAirCtr =  OdorAirCtr.get()

  print('\n ================== ')
  print('\nStarting 1%% CO2 period at %s'%time_stamps[-1][1])
  print('\nCurrSett clean air controller %s'%Settings_CleanAirCtr)
  print('\nCurrSett CO2 controller %s'%Settings_CO2_Ctr)
  print('\nCurrSett odor air controller %s'%Settings_OdorAirCtr)

  print('\n ========= ')
  print('\n LED stimuli starting')
  print('\n ========= ')

  #Displays a range of test stimuli
  for i in range(len(gray_VGR_black)):
    #Add timestamp for stimulus
    add_timestamp("start_" + "CO2 01_" + gray_VGR_black[i] )
    #Upload stimulus
    arduino_upload(gray_VGR_black[i])
    #Leave LED stimulus on for the StimTime
    time.sleep(StimTime-uploadTime)
    #Add timestamp for stimulus
    add_timestamp("end_" + "CO2 01_" + gray_VGR_black[i])
    #outputs time stamp to console
    print ("-- Stimulus period complete at %s --"%time_stamps[-1][1])
    #LEDs controller set to neutral state
    arduino_upload(arduino_neutral)
    #Leave LEDs in neutral state for the interStimTime
    time.sleep(interStimTime-uploadTime)

  #--------------------------------------------------------------------------------

  #Stops the flow of CO2
  change_setpoint(CleanAirCtr, 190.0)
  change_setpoint(CO2_Ctr, 0.0)
  change_setpoint(OdorAirCtr, 0.0)

  #Wait 5 seconds to allow mass flow controllers to bring the flow to the set point
  time.sleep(measureBuffer) 

  #Get current setting for output to console
  Settings_CleanAirCtr = CleanAirCtr.get()
  Settings_CO2_Ctr = 	 CO2_Ctr.get()
  Settings_OdorAirCtr =  OdorAirCtr.get()

  print('\n ================== ')
  print('\nStarting clean air period at %s'%time_stamps[-1][1])
  print('\nCurrSett clean air controller %s'%Settings_CleanAirCtr)
  print('\nCurrSett CO2 controller %s'%Settings_CO2_Ctr)
  print('\nCurrSett odor air controller %s'%Settings_OdorAirCtr)

  print('\n ========= ')
  print('\n LED stimuli starting')
  print('\n ========= ')

  #Displays a range of test stimuli
  for i in range(len(gray_VGR_black)):
    #Add timestamp for stimulus
    add_timestamp("start_" + "clean air_" + gray_VGR_black[i] )
    #Upload stimulus
    arduino_upload(gray_VGR_black[i])
    #Leave LED stimulus on for the StimTime
    time.sleep(StimTime-uploadTime)
    #Add timestamp for stimulus
    add_timestamp("end_" + "clean air_" + gray_VGR_black[i])
    #outputs time stamp to console
    print ("-- Stimulus period complete at %s --"%time_stamps[-1][1])
    #LEDs controller set to neutral state
    arduino_upload(arduino_neutral)
    #Leave LEDs in neutral state for the interStimTime
    time.sleep(interStimTime-uploadTime)

  #--------------------------------------------------------------------------------
  
  #Add time stamp for the end of the experiment
  add_timestamp('end_Exp_' + arduino_off)

  print('\n ================== ')
  print('\nExperiment complete')

  #Stops recording from cameras
  os.system("python ~/ros/flydra-kinetic/src/ros_flydra/scripts/simple-console.py --stop-saving-hdf5")

  #Stop the flow and close the flow controllers
  change_setpoint(CleanAirCtr, 0.0)
  change_setpoint(CO2_Ctr, 0.0)
  change_setpoint(OdorAirCtr, 0.0)
  CleanAirCtr.close()
  CO2_Ctr.close()
  OdorAirCtr.close()
  print('\nFlow controllers stopped at %s'%time_stamps[-1][1])

  #Calculate experimental duration
  duration = time_stamps[-1][2] - time_stamps[2][2]
  print('Total experimental duration: %s seconds'%duration)

  #All LEDs set to off
  print('\nSetting Arduino LED controller so all LEDs are off')
  arduino_upload(arduino_off)

  #save timestamps to output file
  save_timestamps()

  print('\n ================== \n')


except Exception as e:
  print(" * ERROR: %s"%e)


