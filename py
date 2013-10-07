from scitools.std import *
import numpy as np

def extract(filename):
    infile = open(filename,'r')
    column1=[]
    column2=[]
    column3=[]
    column4=[]

    for line in infile:
        if len(line)<2:
            break
        values= line.split()
        
        column1.append(float(values[0]))
        column2.append(float(values[1]))
        column3.append(float(values[2]))
 	column4.append(float(values[3]))
    infile.close()
    return array(column1), array(column2), array(column3),array(column4)


fil='001.txt'
omega_r=0.01
rho,v1,v2,v3 = extract(fil)
n=len(rho)
plot(rho,v1,rho,v2,rho,v3,title=('omega_r=%g'% omega_r),legend=['v1','v2','v3'],xlabel='rho',ylabel='psi',hardcopy= '../../../../../../utdanning/compphys/prosjekt2/plot001.png')
