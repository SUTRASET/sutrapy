#__author__ = 'sbginger'
import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess
import pandas as pd
import re
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
def open_read_well_file(fname):
    wellnodes = []
#    print type(wellnodes), wellnodes
    # First read the inp17 file to determine the number and index of all the lens well nodes
    inp17 = open(fname, 'r')
    while True:
        theline = inp17.readline()
        if len(theline) == 0:              # If there are no more lines
            break                          #     leave the loop
            #    theline.strip()
    #    print theline
        if theline.__contains__("Lens"):
            theline.strip()
            t = theline.split()
    #        print(t)
            wellnodes.extend([abs(int(t[0]))])
    #print wellnodes
    print ("Total # of well nodes =",wellnodes.__len__())
    inp17.close()
    return wellnodes

def create_timestep_list(mynewhandle):
    while True:
        theline = mynewhandle.readline()
        theline.strip()
    #    print theline
        if theline.startswith('## 2-D, REGULAR MESH'):
            z=theline.split()
            y_list=re.findall(r"\d+\.?\d*",z[5])
            y_element=int(y_list[0])
            x_list=re.findall(r"\d+\.?\d*",z[6])
            x_element=int(x_list[0])
           
        if theline.startswith('## NODEWISE'):
            theline.strip()
            t = theline.split()
    #        print(t)
            totaltimesteps = int(t[3])
            print ("Total # of time steps =", totaltimesteps-1)
            interim = [[0 for x in range(9)] for x in range(totaltimesteps)]
            timesteps = [[0 for x in range(2)] for x in range(totaltimesteps)]
    #
        if theline.startswith('##   --------------'):
            for num in range (0,totaltimesteps):
                theline = mynewhandle.readline()
                theline.strip()
    #            print theline
                interim[num] = theline.split()
                timesteps[num][0] = int(interim[num][1])
                timesteps[num][1] = float(interim[num][2])
                #print(timesteps[num])
            break
    return x_element,y_element,timesteps,t

# def resultreader( file_data, wellnodes, numberofnodes):
def resultreader( file_data, numberofnodes):
    
    allnodes       = []
    allnodes_data  = []
    while True:
        theline = file_data.readline()
        theline.strip()
#        print timestep_id, theline
        
        if theline.startswith('## TIME STEP'):
            #print theline
            theline.strip()
            theline = file_data.readline()
            theline = file_data.readline()
            
            # if theline.startswith('##   Node              X'):
            element=theline.split()
                #theline = file_data.readline()
            allnodes= [file_data.readline().strip().split() for num in range (0, numberofnodes)]
            allnodes_data=np.array([np.array(allnodes)])
            break
                    #allnodes[num] = [float(y) for lst in (theline.split() for x in allnodes[num]) for y in lst]
#                print allnodes[152520]
#                print allnodes[152521]
#                print allnodes[152522]
                # wellnodes = [[float(y) for y in nodeline] for i, nodeline in enumerate(allnodes) if i+1 in wellnodes]
#                print wellnodes[0]
#                print wellnodes[1]
#                print wellnodes[2]                
                # break

    # return wellnodes
    return allnodes,allnodes_data,element


def concentration_averager(all_wellnode_data,weighted_c_ts,homogeneous):
    for ts_data  in all_wellnode_data:
        x, y, z, p, c, s = zip(*ts_data)
    # The following section calculates the weighting based on the permeability distribution used in the Roi model
        x_array = np.sin(np.array(x)/100.*np.pi)
        y_array = np.sin(np.array(y)/100.*np.pi)
        z_array = np.sin(np.array(z)/4.*np.pi)
        c_array = np.array(c)
        LogRange = 3
        sin_func = 10**(LogRange*(1-abs((x_array*y_array*z_array))))
        Estimated_k = (0.00000000027/(10**LogRange))*sin_func
        weight = Estimated_k / np.mean(Estimated_k)
        for ii in range(len(c_array)):
            if c_array[ii] < 0.000055:
                print (c_array[ii])
                c_array[ii] = 0.000055
                print (c_array[ii])
        if homogeneous:
            c_weight = c_array * 1
        else:
            c_weight = c_array*weight
    # This line converts SUTRA concentration to equivalent Chloride concentration
#        print np.mean(c)/.0357*19600,np.mean(c_weight)/.0357*19600

        weighted_c_ts = np.mean(c_weight)/.0357*19600

    return weighted_c_ts



# fname_wellnodes = "file.inp17"
# wellnodes = open_read_well_file(fname_wellnodes)

# numberofnodes is a hardwired value that should be changed depending on the mesh
# numberofnodes = 491647

# fname1720 = "\\output\\Washover\\roi_sourcewashover_dec08.nod"

# fname1820 = "\\output\\Recovery_calibrated\\output\\roi_sourcewashover_2yr.nod"
# fname1821 = "\\output\\Recovery_wout_withdrawal\\roi_sourcewashover_2yr.nod"
# fname1822 = "\\output\\No_washover\\roi_sourcewashover_2yr.nod"
# fname1823 = "\\output\\Recovery_wout_artificial_recharge\\roi_sourcewashover_2yr.nod"
# fname1824 = "\\output\\Undisturbed_recovery\\roi_sourcewashover_2yr.nod"
# fname1825 = "\\output\\No_wash_No_art\\roi_sourcewashover_2yr.nod"
# fname1826 = "\\output\\Undisturbed_recovery\\roi_sourcewashover_2yr.nod"
# fname1827 = "\\output\\Average_stress\\roi_sourcewashover_2yr.nod"
# fname1828 = "\\output\\Recovery_after_wetseason\\roi_sourcewashover_2yr.nod"
fname="C:/Project/MDBA/data_deliverable/modelling_benchmark_submit/PART1.nod"
# fname="C:/SUTRA/SutraLab-master/example/SUTRA_examples/2D/Island2D/Island2D.nod"
# fname="C:/SUTRA/SutraLab-master/example/SUTRA_examples/2D/Henry/Henry2D.nod"

homogeneous = True #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# file_data = open(fname1820, 'r')     #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
file_data = open(fname, 'r')     #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
x_element,y_element,timesteps,t= create_timestep_list(file_data)
numberofnodes = x_element*y_element
ts_array = np.array(timesteps)
t_start = 0.0E+08-(0*24*3600.)-(0*24*3600)
ts_array_zeroed = (ts_array[:,1] - t_start)/3600/24

# all_wellnode_data = []
time_series_average = []
allnodes_data_timeseries= []
weighted_c_ts = []
#The following line must be changed to include all of the timesteps desired
i=0
for timestep_id, timestep_clock in timesteps:
    allnodes,allnodes_data,element = resultreader(file_data, numberofnodes)
    allnodes_data_timeseries.append(allnodes)
    i=i+1
    print(i)
allnodes_data_timeseries_array=np.array(allnodes_data_timeseries)
allnodes_data_timeseries_array=allnodes_data_timeseries_array.astype(np.float)
element=element[1:]

X_number=element.index('X')
Y_number=element.index('Y')
X=allnodes_data_timeseries_array[0,:,X_number]
Y=allnodes_data_timeseries_array[0,:,Y_number]
X_min=min(allnodes_data_timeseries_array[1,:,X_number])
X_max=max(allnodes_data_timeseries_array[1,:,X_number])
Y_min=min(allnodes_data_timeseries_array[1,:,Y_number])
Y_max=max(allnodes_data_timeseries_array[1,:,Y_number])
Pressure_number=element.index('Pressure')
Saturation_number=element.index('Saturation')
Concentration_number=element.index('Concentration')

sx=np.linspace(X_min,X_max,x_element)
sy=np.linspace(Y_min,Y_max,y_element)
# sx=np.unique(allnodes_data_timeseries_array[0,:,X_number])
# sy=np.unique(allnodes_data_timeseries_array[0,:,Y_number])
xv,yv=np.meshgrid(sx,sy)
for i in range(len(timesteps)):
    Saturation=allnodes_data_timeseries_array[i,:,Saturation_number]
    Saturation_interpolated=griddata((X,Y),Saturation,(xv,yv))   
    # Saturation2d=Saturation.reshape(len(sy),len(sx),order='F')
    plt.contourf(xv,yv,Saturation_interpolated)
    plt.title('Saturation')
    plt.pause(1)
    # all_wellnode_data.append(wellnodes_data)
    # average_c = concentration_averager(all_wellnode_data,weighted_c_ts,homogeneous)
    # time_series_average.append(average_c)
    
for i in range(len(timesteps)):
    Pressure=allnodes_data_timeseries_array[i,:,Pressure_number]
    Pressure_interpolated=griddata((X,Y),Pressure,(xv,yv))   
    # Saturation2d=Saturation.reshape(len(sy),len(sx),order='F')
    plt.contourf(xv,yv,Pressure_interpolated)
    plt.title('Pressure')    
    plt.pause(1)

for i in range(len(timesteps)):
    Concentration=allnodes_data_timeseries_array[i,:,Concentration_number]
    Concentration_interpolated=griddata((X,Y),Concentration,(xv,yv))   
    # Saturation2d=Saturation.reshape(len(sy),len(sx),order='F')
    plt.contourf(xv,yv,Concentration_interpolated)
    plt.title('Concentration')    
    plt.pause(1)

# allnodes_data_array=np.array(allnodes_data)
# df = pd.DataFrame(data=numpy_data, index=["row1", "row2"], columns=element)

# print homogeneous, "homogeneous flag"
# print len(all_wellnode_data), 'number of TS'

# #print results to a file

# out1720 = "\\output\\Washover\\cl_time.dat"

# out1820 = "\\output\\Recovery_calibrated\\output\\cl_time.dat"
# out1821 = "\\output\\Recovery_wout_withdrawal\\cl_time.dat"
# out1822 = "\\output\\No_washover\\cl_time.dat"
# out1823 = "\\output\\Recovery_wout_artificial_recharge\\cl_time.dat"
# out1824 = "\\output\\Undisturbed_recovery\\cl_time.dat"
# out1825 = "\\output\\No_wash_No_art\\cl_time.dat"
# out1826 = "\\output\\Undisturbed_recovery\\cl_time.dat"
# out1827 = "\\output\\Average_stress\\cl_time.dat"
# out1828 = "\\output\\Recovery_after_wetseason\\cl_time.dat"
# out="C:/Project/MDBA/data_deliverable/modelling_benchmark_submit/PART1.dat"
# # activeout = out1820      #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# activeout = out    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# outfile = open(activeout, 'w')


# for ii in range(len(time_series_average)):
#     outfile.write(str(ts_array_zeroed[ii])+'\t'+str(time_series_average[ii])+'\n') #
# outfile.close

# #read in observed chloride data
# fname = "E:\\Users\\sbginger\\Roi-Namur\\Python\\roi_wash_rec.csv"
# file_data = open(fname, 'r')
# observed_c = []
# while True:
#         theline = file_data.readline()
#         if len(theline) == 0:              # If there are no more lines
#             break                          #     leave the loop
#         theline = theline.strip()
# #        print theline
#         t = theline.split(',')
# #        print t
#         t = [np.float(x) if x!='' else np.float('Nan') for x in t]
#         observed_c.append(t)
#         observed_c_array = np.array(observed_c)
# #        print(observed_c_array)
# #        print (observed_c)
#         #print observed_c_array

# file_data.close()
# observed_c_array = zip(*observed_c_array)
# #print len(observed_c_array), len(observed_c_array[0])
# #start at time zero and convert to days


# fig1 = plt.figure(1, figsize=(6,6))
# ax1 = fig1.add_axes([0.25, 0.15, 0.65, 0.65])
# ax1.set_title('Chloride Washover and Recovery')
# ax1.set_xlabel('Time, in days')
# ax1.set_ylabel('Chloride Concentration, in mg/L')
# #ax1.plot(time, y, label='Data', linestyle='none', marker='o', mfc='black', mew=0)
# for yy in observed_c_array[1:]:
#     ax1.plot(observed_c_array[0], yy, linestyle='none', marker='o', mfc='DarkGray', mew=0)
#     ax1.plot(ts_array_zeroed, time_series_average, linestyle='solid', color='tomato')
# plt.show()
# #cmd = "uedit64.exe"
# #subprocess.call(['uedit64.exe', {out526}])
# #os.system(cmd)
# subprocess.Popen(['uedit64.exe', {activeout}])
