#!/usr/bin/env python
# -*- coding: utf-8 -*-
# *************************************************
#
# Monte Carlo calibration algorithm for ECH2O
#
# -------
# Routine: Subroutines for parameters & obs manip & outputs
# -------
# Author: S. Kuppel
# Created on 10/2016
# -------------------------------------------------

import shutil
import scipy.io as spio
import time, os, glob, sys, copy, subprocess
import csv
from time import localtime, strftime
from calendar import monthrange as lmon
import random
import numpy as np
import datetime
import pandas as pd
from pyDOE import *
from pcraster import *
from datetime import timedelta, datetime

# ----------------------------------------------------------------------------
# -- Generate random set of parameters values of write it


def saveToASCII(data, fileName, savePath_catchinfo, commandLine, format, mask):
    if not os.path.exists(savePath_catchinfo):
        os.makedirs(savePath_catchinfo)
    exec(commandLine)
    data = data.astype(np.double)
    data[np.isnan(mask)] = -9999.0
    line1 = 'ncols         22'
    line2 = 'nrows         30'
    line3 = 'xllcorner     442449.229'
    line4 = 'yllcorner     5798066.25'
    line5 = 'cellsize      500'
    line6 = 'NODATA_value  -9999'
    with open(savePath_catchinfo + fileName + '.asc', 'w') as f:
        f.write(line1+'\n'+line2+'\n'+line3+'\n'+line4+'\n'+line5+'\n'+line6+'\n')

    with open(savePath_catchinfo + fileName + '.asc', 'ab') as f:
        if format=='int':
            np.savetxt(f, data, fmt='%d')
        else:
            np.savetxt(f, data)

def sortNO3inputs(path, param):
    path1 = path + '/outputs.1/'
    if os.path.exists(path+'/catchinfo'):
        shutil.rmtree(path+'/catchinfo')
    shutil.copytree('/home/wusongj/dmc/forNitrate/catchinfo', path+'/catchinfo')

    mask = pcr2numpy(readmap(path1+'Spatial/DEM.map'), np.nan)

    tmp1 = pcr2numpy(readmap(path1+'Spatial/soildepth.L1.map'), np.nan)
    tmp2 = pcr2numpy(readmap(path1+'Spatial/soildepth.L2.map'), np.nan)
    tmp3 = pcr2numpy(readmap(path1+'Spatial/soildepth.map'), np.nan)

    depth_1 = tmp1
    depth_2 = tmp2
    depth_3 = tmp3 - tmp2 - tmp1


    saveToASCII(depth_1, 'depth_soilL1', path+'/catchinfo/', '', 'double', mask)
    saveToASCII(depth_2, 'depth_soilL2', path+'/catchinfo/', '', 'double', mask)
    saveToASCII(depth_3, 'depth_soilL3', path+'/catchinfo/', '', 'double', mask)

    saveToASCII(pcr2numpy(readmap(path1+'tmp/WiltingPoint.map'), np.nan), 'wp_0', path+'/catchinfo/', '', 'double', mask)
    saveToASCII(pcr2numpy(readmap(path1+'tmp/WiltingPoint.map'), np.nan), 'wp_1', path+'/catchinfo/', '', 'double', mask)
    saveToASCII(pcr2numpy(readmap(path1+'tmp/WiltingPoint.map'), np.nan), 'wp_2', path+'/catchinfo/', '', 'double', mask)
    saveToASCII(pcr2numpy(readmap(path1+'tmp/FCap0.map'), np.nan), 'sat_0', path+'/catchinfo/', '', 'double', mask)
    saveToASCII(pcr2numpy(readmap(path1+'tmp/FCap1.map'), np.nan), 'sat_1', path+'/catchinfo/', '', 'double', mask)
    saveToASCII(pcr2numpy(readmap(path1+'tmp/FCap2.map'), np.nan), 'sat_2', path+'/catchinfo/', '', 'double', mask)
    saveToASCII(pcr2numpy(readmap(path1+'tmp/poros_0.map'), np.nan), 'poros_0', path+'/catchinfo/', '', 'double', mask)
    saveToASCII(pcr2numpy(readmap(path1+'tmp/poros_1.map'), np.nan), 'poros_1', path+'/catchinfo/', '', 'double', mask) 
    saveToASCII(pcr2numpy(readmap(path1+'tmp/poros_2.map'), np.nan), 'poros_2', path+'/catchinfo/', '', 'double', mask)
    #saveToASCII(pcr2numpy(readmap(path1+'Spatial/srfovf_ratio.map'), np.nan), 'srfovf_ratio', path+'/catchinfo/', '', 'double', mask)

    len_tmp = len(np.fromfile(path1 + 'Outputs_tmp/ChnS.bin'))
    np.full(len_tmp, 0.0).tofile(path1 + 'Outputs_tmp/zeros.bin')

    conc0 = param[-7]
    conc1 = 2
    conc2 = 1
    conc3 = 1
    p0 = pcr2numpy(readmap(path1+'Spatial/p_0.map'), np.nan)
    p1 = pcr2numpy(readmap(path1+'Spatial/p_1.map'), np.nan)
    p2 = pcr2numpy(readmap(path1+'Spatial/p_2.map'), np.nan)
    p3 = pcr2numpy(readmap(path1+'Spatial/p_3.map'), np.nan)
    tmp = conc0*p0+conc1*p1+conc2*p2+conc3*p3
    saveToASCII(tmp, 'init_concIN', path+'/catchinfo/', '', 'double', mask)
    saveToASCII(tmp, 'init_concON', path+'/catchinfo/', '', 'double', mask)
    
    """"""
    srfMix_ratio = param[-6]
    with open(path + '/wqm_config.nml', 'r') as f:
        lines = f.readlines()
    for i in range(len(lines)):
        if 'SrfOVFmixratio =  ' in lines[i]:
            lines[i] = 'SrfOVFmixratio =  ' + str(srfMix_ratio) + '\n'
    with open(path + '/wqm_config.nml', 'w') as f:
        f.writelines(lines)
    

    
    crop_fertilisation = param[-5]
    crop_fertd1 = int(param[-2])
    wetland_fertd1 = int(param[-1])
    with open(path+'/catchinfo/cropdata.txt', 'r') as f:
        lines = f.readlines()
    for i in range(len(lines)):
        if 'crop_pre2000' in lines[i]:
            lines[i] = 'crop_pre2000	1	' + str(crop_fertilisation) + '	'+str(crop_fertd1)+'	0.4	0	0	0	15	110	0.3	0	0	0	0.5	60	20	260	0.3	0.8	30	10	0.5	0.05	0.3	0.01	85	0	220	0	0	0' +'\n'
        elif 'peatland' in lines[i]:
            lines[i] = 'peatland	5	40	'+str(wetland_fertd1)+'	0.2	0	0	0	30	120	0.08	8.6	220	0.1	0.5	60	59	290	0.3	0.8	1	10	0.8	0.03	0.4	0.001	120	0	350	0	0	0'
    with open(path+'/catchinfo/cropdata.txt', 'w') as f:
        f.writelines(lines)

    
    exp_depth_adentri = param[-4]
    frac_activeGW = param[-3]
    with open(path+'/wqm_config.nml', 'r') as f:
        lines = f.readlines()
    for i in range(len(lines)):
        if 'exp_depth_adentri' in lines[i]:
            lines[i] = 'exp_depth_adentri = ' + str(exp_depth_adentri) +'\n'
        if 'frac_activeGW' in lines[i]:
            lines[i] = 'frac_activeGW = ' + str(frac_activeGW) +'\n'
    with open(path+'/wqm_config.nml', 'w') as f:
        f.writelines(lines)



def run_NO3Module(mode):
    path_old = os.getcwd()
    path = os.path.abspath(os.path.join(path_old, "../.."))
    path1 = path + '/outputs.1/'
    path2 = path + '/outputs.1/tmp_no3/'

    if os.path.exists(path2):
        shutil.rmtree(path2)
    os.mkdir(path2)
    os.chdir(path)

    pvalue_new = np.loadtxt('./Input_ensemble/0.0_bestParams_no3.txt')
    
    createParam_no3(pvalue_new)

    sortNO3inputs(path, pvalue_new)

    subprocess.run('./nitrateModule', stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    #subprocess.run('./nitrateModule')

    for file in os.listdir(path2):
        arr = np.fromfile(path2+file)
        with open(path1+file,'a') as f_out:
            arr.tofile(f_out)

    os.chdir(path_old)

def createParam_no3(param):
    file = './wqm_parameter.nml'
    nlanduseType = 4
    with open(file, 'r') as f:
        lines = f.readlines()
    for i in range(len(lines)):   
        if '!> instream_denitrification' in lines[i]:
            idx = 0
            lines[i+10] = 'Nparameters('+str(idx+1)+')%pvalue = ' + ','.join([str(a) for a in param[idx*nlanduseType:(1+idx)*nlanduseType]]) + '\n'
        if '!> autotrophicuptk_rate' in lines[i]:
            idx = 1
            lines[i+10] = 'Nparameters('+str(idx+1)+')%pvalue = ' + ','.join([str(a) for a in param[idx*nlanduseType:(1+idx)*nlanduseType]]) + '\n'
        if '!> netprimaryprod_rate' in lines[i]:
            idx = 2
            lines[i+10] = 'Nparameters('+str(idx+1)+')%pvalue = ' + ','.join([str(a) for a in param[idx*nlanduseType:(1+idx)*nlanduseType]]) + '\n'
        if '!> denitrification_soil' in lines[i]:
            idx = 3
            lines[i+10] = 'Nparameters('+str(idx+1)+')%pvalue = ' + ','.join([str(a) for a in param[idx*nlanduseType:(1+idx)*nlanduseType]]) + '\n'
        if '!> mineralisationN_rate' in lines[i]:
            idx = 4
            lines[i+10] = 'Nparameters('+str(idx+1)+')%pvalue = ' + ','.join([str(a) for a in param[idx*nlanduseType:(1+idx)*nlanduseType]]) + '\n'
        if '!> degradationN_rate' in lines[i]:
            idx = 5
            lines[i+10] = 'Nparameters('+str(idx+1)+')%pvalue = ' + ','.join([str(a) for a in param[idx*nlanduseType:(1+idx)*nlanduseType]]) + '\n'
        if '!> dissolutionN_rate' in lines[i]:
            idx = 6
            lines[i+10] = 'Nparameters('+str(idx+1)+')%pvalue = ' + ','.join([str(a) for a in param[idx*nlanduseType:(1+idx)*nlanduseType]]) + '\n'
    with open(file, 'w') as f:
        f.writelines(lines)

def gen_paras(Opti, Config):

    #Opti.xtot = np.arange(1,Opti.nsamptot+1)

    # Latin Hypercube sampling
    if Config.sampling in ['LHS','LHS_m','LHS_r']:

        print('...using a latin hypercube sampling...')

        # 'Normal' LHS
        if Config.sampling=='LHS':
            # First, get the hypercube samples, ranging from 0 to 1
            mat = np.transpose(lhs(Opti.nvar,samples=Opti.nsamptot))

        # LHS with additional criterion: maixmim distane between samples
        elif Config.sampling=='LHS_m':
            print('...with maximin criterion -- it will take longer and a lot of memory')
            # First, get the hypercube samples, ranging from 0 to 1
            mat = np.transpose(lhs(Opti.nvar,samples=Opti.nsamptot,criterion='m'))

        # LHS with additional criterion: maixmim distane between samples
        elif Config.sampling=='LHS_r':
            print('...with correlation criterion -- it will take longer and a lot of memory')
            # First, get the hypercube samples, ranging from 0 to 1
            mat = np.transpose(lhs(Opti.nvar,samples=Opti.nsamptot,criterion='corr'))


        print('...LHS matrix generated...')
        
        # Second, scale with the actual range
        for i in range(Opti.nvar):
            # Log transform where needed
            if Opti.log[i]==1:                
                tmp = 10**(mat[i]*np.log10(Opti.max[i]/Opti.min[i])+np.log10(Opti.min[i]))
            else:
                tmp = mat[i]*(Opti.max[i]-Opti.min[i]) + Opti.min[i]
                
            
            if i==0:
                Opti.xtot = copy.copy(tmp)
            else:
                Opti.xtot = np.vstack((Opti.xtot,tmp))

    #-- Uniform distribution
    elif Config.sampling=='uniform':

        print('...using uniform distributions...')

        for i in range(Opti.nvar):
            # Log transform where needed
            if Opti.log[i]==1:
                if i==0:
                    Opti.xtot = 10**(np.random.uniform(np.log10(Opti.min[i]),np.log10(Opti.max[i]),Opti.nsamptot))
                else:
                    Opti.xtot = np.vstack((Opti.xtot,10**(np.random.uniform(np.log10(Opti.min[i]),
                                                                            np.log10(Opti.max[i]),
                                                                            Opti.nsamptot))))
            else:
                if i==0:
                    Opti.xtot = np.random.uniform(Opti.min[i],Opti.max[i],Opti.nsamptot)
                else:
                    Opti.xtot = np.vstack((Opti.xtot,np.random.uniform(Opti.min[i],Opti.max[i],Opti.nsamptot)))
        #print i, Opti.names[i], Opti.x[i]

    else:
        sys.exit('No proper sampling method selected!')


    print('Parameters sampling done!')

    # -- Reshape for output
    #Opti.xtot = np.transpose(Opti.xtot)#.shape((Opti.nvar,Opti)
    print(Opti.xtot.shape)

    print('')
    print('Writing in '+Config.ncpu+' files...('+str(Opti.nit)+' sets each)')

    # -- Write one file per parallel job
    for i in range(int(Config.ncpu)):
        f_out = Config.FILE_PAR+str(i+1)+'.txt'
        k = i*Opti.nit
        #print str(k)+','+str(k+Opti.nit)
        with open(f_out,'w') as fw:
            for j in range(Opti.nvar):
                #print i
                #print j
                #print k
                tmp=[a for a in Opti.xtot[j][k:k+Opti.nit]]
                #print len(tmp)
                fw.write(Opti.names[j]+','+','.join([str(a) for a in tmp])+'\n')

    # Write on file giving parameters range, log...(for later plots)
    f_out = Config.FILE_PAR+'char.txt'
    with open(f_out,'w') as fw:
        # fw.write('Sample,'+','.join(Opti.names)+'\n')
        fw.write('Names,Min,Max,Log\n')
        for i in range(Opti.nvar):
            fw.write(','.join([Opti.names[i],str(Opti.min[i]),str(Opti.max[i]),str(Opti.log[i])])+'\n')

    print('')
    
    #np.savetxt(f_out,np.transpose(Opti.xtot))#,header= 

# ----------------------------------------------------------------------------
# -- Write parameters values file

def get_par(Opti,Config):
  
    # Open one file for all samples
    f_in = Config.FILE_PAR+Config.numsim+'.txt'
    print(f_in)
    Opti.xpar = np.genfromtxt(f_in,delimiter=',',unpack=True)[1::]
    print(Opti.xpar.shape)

# ----------------------------------------------------------------------------
# -- Write parameters values file

def output_par(Opti,Config, it):
  
    # Open one file for all samples
    if Config.initpar==0:
        Config.f_par = Config.PATH_OUT+'/Parameters.txt'
        if Config.restart == 0:
            with open(Config.f_par,'w') as f_in:
                f_in.write('Iteration,'+','.join(Opti.names)+'\n')
        Config.initpar = 1

    with open(Config.f_par,'a') as f_in:
        f_in.write(str(it+1)+','+','.join([str(x) for x in Opti.x])+'\n')

# ----------------------------------------------------------------------------
# -- Updating inputs for ECH2O

def create_inputs(Opti, Paras, Site, Config, it):

    # Small switch not to read the vegetation params file every time
    readveg=1
    
    # -- Get the parameter sample of the current iteration
    Opti.x = Opti.xpar[it]

    outmapM = Site.bmaps['chanmask']*0

    outmapR = Site.bmaps['chanmask']*0

    for pname in Paras.names:

        #print pname

        ## - Mapped parameters
        if Paras.ref[pname]['veg']==0:
            # Soil unit dependence
            #if (Paras.ref[pname]['soil']==1 and pname!='Khoriz') or (Paras.ref[pname]['soil']==1 and Config.mode==2):
            if Paras.ref[pname]['soil']==1:
                #print 'Soil dependent !!'
                outmap = Site.bmaps['unit']*1
#                outmap = Site.bmaps['unit']*0
                # Read each soil map unit and apply param value                
                for im in range(Site.ns):
#                    outmap+= Site.bmaps[Site.soils[im]]*Opti.x[Paras.ind[pname]+int(im)]
                    outmap*= Opti.x[Paras.ind[pname][im]]**Site.bmaps[Site.soils[im]]

                if Opti.simRock == 1:
                    # Taking into account rock/scree: micro-topsoil, low poros and fixed anisotropy
                    if pname=='HLayer1':
                        outmap = outmap*(Site.bmaps['unit']-Site.bmaps['rock'])+Site.bmaps['rock']*0.001 
                    # if pname=='Porosity':
                    #    outmap = outmap*(Site.bmaps['unit']-Site.bmaps['rock'])+Site.bmaps['rock']*0.25 
                    if pname=='Khoriz':
                        outmap = outmap*(Site.bmaps['unit']-Site.bmaps['rock'])+Site.bmaps['rock']*0.000001 
                    if pname=='Anisotropy':
                        outmap = outmap*(Site.bmaps['unit']-Site.bmaps['rock'])+Site.bmaps['rock']*0.1 

                report(outmap,Config.PATH_SPA+'/'+Paras.ref[pname]['file']+'.map')

            # No spatial/veg dependence, but channel stuff
            else:
                #print 'Not dependent !!'
                if Paras.ref[pname]['file'] in ['chanmanningn']:
                    if pname=='manningRiv_all':
                        outmapM += (Site.bmaps['chanmask'] - Site.bmaps['chanmask_wetland'])*Opti.x[Paras.ind[pname]]
                    else:
                        outmapM += Site.bmaps['chanmask_wetland']*Opti.x[Paras.ind[pname]]
        
                elif Paras.ref[pname]['file'] in ['chanwidth']:
                    outmap = Site.bmaps['chanmask']*Opti.x[Paras.ind[pname]]
                elif Paras.ref[pname]['file'] in ['chanrough']:
                    if Opti.wetland == 1:
                        if pname=='chanrough_all':
                            outmapR += (Site.bmaps['chanmask'] - Site.bmaps['chanmask_wetland'])*Opti.x[Paras.ind[pname]]
                        else:
                            outmapR += Site.bmaps['chanmask_wetland']*Opti.x[Paras.ind[pname]]
                    else:
                        if pname=='chanrough_all':
                            outmapR += Site.bmaps['chanmask']*Opti.x[Paras.ind[pname]]
                elif Paras.ref[pname]['file'] == 'chanparam':
                    outmap = Site.bmaps['chanmask']*Opti.x[Paras.ind[pname]]
                else:
                    outmap = Site.bmaps['unit']*Opti.x[Paras.ind[pname]]

                if Paras.ref[pname]['file'] in ['chanmanningn']:
                    report(outmapM,Config.PATH_SPA+'/'+Paras.ref[pname]['file']+'.map')
                elif Paras.ref[pname]['file'] in ['chanrough']:
                    report(outmapR,Config.PATH_SPA+'/'+Paras.ref[pname]['file']+'.map')
                else:
                    report(outmap,Config.PATH_SPA+'/'+Paras.ref[pname]['file']+'.map')

        ## - Vegetation parameters
        else:
            # Change the value based on param name correspondance
            vegnew = copy.copy(Opti.vref)
#            print(Opti.vref)
#            print(Site.nv)
            for iv in range(Site.nv):
#                print(iv)
#                print(pname)
#                vegnew[iv][vegnew['name'].index(pname)] = str(Opti.x[Paras.ind[pname]+int(iv)])

                vegnew[iv][vegnew['name'].index(pname)] = str(Opti.x[Paras.ind[pname][iv]])


    # ------------------------------------------------------------------------------
    # Finalizing the preps....
    # Check for initial condition/other parameter dependence
    # Back up relevant maps
                
    # Initial soil water content: needs porosity profile and depth to have consistent
    # initial soil moisture for each layer

    # print(Paras.names)
    if 'Porosity0' in Paras.names:
        poros = readmap(Config.PATH_SPA+'/'+Paras.ref['Porosity0']['file']+'.map')
    elif 'Porosity' in Paras.names:
        poros = readmap(Config.PATH_SPA+'/'+Paras.ref['Porosity']['file']+'.map')
        #print(Config.PATH_SPA+'/'+Paras.ref['Porosity']['file']+'.map')
    else:
        poros = readmap(Config.PATH_SPA+'/poros.map')

    if 'kPorosity' in Paras.names:
#        print("kPorosity activated")
        kporos = readmap(Config.PATH_SPA+'/'+Paras.ref['kPorosity']['file']+'.map')

        if 'HLayer1' in Paras.names:
            dL1 = readmap(Config.PATH_SPA+'/'+Paras.ref['HLayer1']['file']+'.map')
        else:
            dL1 = readmap(Config.PATH_SPA+'/soildepth.L1.map')
        if 'HLayer2' in Paras.names:
            dL2 = readmap(Config.PATH_SPA+'/'+Paras.ref['HLayer2']['file']+'.map')
        else:
            dL2 = readmap(Config.PATH_SPA+'/soildepth.L2.map')
        if 'Depth' in Paras.names:
            dTot = readmap(Config.PATH_SPA+'/'+Paras.ref['Depth']['file']+'.map')
        else:
            dTot = readmap(Config.PATH_SPA+'/soildepth.map')

        # Layer-integrated values from profile
        porosL1 = kporos*poros*(1-exp(-dL1/kporos))/dL1
        porosL2 = kporos*poros*(exp(-dL1/kporos)-exp(-(dL1+dL2)/kporos))/dL2
        porosL3 = kporos*poros*(exp(-(dL1+dL2)/kporos)-exp(-dTot/kporos))/(dTot-dL1-dL2)
        #print("Checked porosity")
    else:

        porosL1 = poros
        porosL2 = poros
        porosL3 = poros

    #report(porosL1*dL1,Config.PATH_SPA+'/storageL1.map')
    
    # -- Use a fraction of these porosities as initial soil moisture
    report(porosL1*0.95,Config.PATH_SPA+'/SWC.L1.map')
    report(porosL2*0.9,Config.PATH_SPA+'/SWC.L2.map')
    report(porosL3*0.9,Config.PATH_SPA+'/SWC.L3.map')       
            
    ## - Finalizing soil parameterization
    # Check that initial soil moisture is not smaller residual soil 
    # tbd... for now just pick the porosity and thetar reange wisely enough

    ## - Finalizing the vegetation parameterization
    if Paras.isveg > 0:
    # Equalize leaf turnover and additional turnover due to water and/or temperature stress
        #for iv in range(Site.nv):
        #    vegnew[iv][vegnew['name'].index('TurnovL_MWS')] = copy.copy(vegnew[iv][vegnew['name'].index('TurnovL')])
        #    vegnew[iv][vegnew['name'].index('TurnovL_MTS')] = copy.copy(vegnew[iv][vegnew['name'].index('TurnovL')])
    # Write the vegetation params file (if there's at least one veg param)
        vegfile = open(Config.PATH_SPA+'/'+Site.vfile,'w')
        vegfile.write('\t'.join(Opti.vref['header'])+'\n')
        for iv in range(Site.nv):
            vegfile.write('\t'.join(vegnew[iv]) +'\n')
        vegfile.write('\t'.join(Opti.vref['footer'])+'\n')
        vegfile.write('\t'.join(Opti.vref['name'])+'\n')
        vegfile.close()

# ----------------------------------------------------------------------------
# -- Check is ECH2O ran properly

def runOK(Data, Opti, Config):

    isOK = 1

    # 1. Check it ran
    if len(glob.glob(Config.PATH_EXEC+'/BasinSummary.txt')) == 0 or os.stat(Config.PATH_EXEC+'/BasinSummary.txt').st_size == 0 :
        #print("Something went wrong, "+Config.PATH_EXEC+"BasinSummary.txt is missing/empty...")
        isOK = 0
    
    else:
        for oname in Data.names:
            if Data.obs[oname]['type']=='Ts' or Data.obs[oname]['type']=='Total':
                f_test = Config.PATH_EXEC+'/'+Data.obs[oname]['sim_file']
                if len(glob.glob(f_test)) == 0:
                    print("Something went wrong, no output for "+oname+" !!")
                    print('(i.e., '+f_test+' is missing...)')
                    isOK = 0    
                    # print 'Outputs are there...'
    
        f_test = Config.PATH_EXEC+'/BasinSummary.txt'
        try:
            tmp = np.genfromtxt(f_test,skip_header=1,unpack=True)
        except ValueError:
            isOK = 0
        else:
            if tmp!=Data.lsim:
                isOK = 0
                print("Something went wrong, output does not match the supposed sim length!")

        """
        f_test = Config.PATH_EXEC+'/BasinSummary.txt'
        try:
            tmp = np.genfromtxt(f_test,skip_header=1,unpack=True)[0]
        except ValueError:
            isOK = 0
        else:
            # print str(len(tmp2))
            # print str(Data.lsim)
            if type(tmp)==float or type(tmp)==np.float64 or type(tmp)==np.float32:
                isOK = 0
                print("Something went wrong, output of length 1 !")
                
            if len(tmp)!=Data.lsim: # & len(tmp2)==Data.lsim:
                # print 'It ran until the end !'
                isOK = 0
                print("Something went wrong, output does not match the supposed sim length!")
                print('Output: '+str(len(tmp))+' , supposed to be: '+str(Data.lsim))
        """

    return isOK
        
# ----------------------------------------------------------------------------
# -- Manage outputs        

def manage_outputs(Data, Opti, Config, it):

    # -- Report the full BasinSummary.txt files?
    #if Config.repBS == 1:
    #    os.system('mv '+Config.PATH_EXEC+'/BasinSummary.txt '+Config.PATH_OUT+'/BasinSummary_run'+str(it+1)+'.txt')
            
    # -- Group the output files in one across simulations, 
    #    separating by observations points and veg type where it applies
    for oname in Data.names:
        if (Data.obs[oname]['type']!='map' or Data.obs[oname]['type']!='mapTs') and (it==0 or Opti.begfail==1):
            # Historic time series file names
            if Config.restart == 0:
                Data.obs[oname]['sim_hist'] = Config.PATH_OUT+'/'+oname+'.tab'
            # Header of files
            # Exclude header. No need for binary file. Modified by Songjun
            open(Data.obs[oname]['sim_hist'],'w')

            #with open(Data.obs[oname]['sim_hist'],'w') as f_out:
            #    f_out.write('Sample,'+','.join([str(i+1) for i in range(Config.trimL)])+'\n')

    # Reinit begfail (otherwise will never write all!)
    Opti.begfail = 0

    # Save current run outputs (and delete the file to relieve Maxwell...)
    for oname in Data.names:
        startDate = (Data.simbeg + timedelta(days=Config.trimB)).date()
        endDate   = startDate + timedelta(days=Config.trimL-1)
        dateList = (pd.date_range(startDate,endDate,freq='M')).date
        dateList = np.append(dateList, endDate)
        dateList = np.append(startDate, dateList)
        idxList = [(d-startDate).days+int(Config.trimB) for d in dateList]
        startList = idxList[:-1]
        endList = idxList[1:]




        # Integrated variables (in BasinSummary.txt)
        if Data.obs[oname]['type']=='Total':
            idx = Data.obs[oname]['sim_pts']-1

            tmp = np.genfromtxt(Data.obs[oname]['sim_file'],delimiter='\t',
                                skip_header=1,unpack=True)[idx]*Data.obs[oname]['conv']

            # Shave off the transient part (if any)
            if Config.trimB > 1:
                tmp = tmp[Config.trimB-1:Config.trimB-1+Config.trimL]
                if len(tmp) != Config.trimL:
                    sys.exit("ERROR -> Problem with output trim: we've got "+str(len(tmp))+
                             ' instead of '+str(Config.trimL))

            with open(Data.obs[oname]['sim_hist'],'a') as f_out:
                f_out.write(str(it+1)+','+','.join([str(j) for j in list(tmp)])+'\n')

        # Time series for calibration:
        if Data.obs[oname]['type']=='Ts_cali':
            hskip = Data.nts+3
            # export multiple points, modified by Songjun    
            #idx = np.argsort(np.array(Data.sim_order))[Data.obs[oname]['sim_pts']-1]+1
            #idx = np.array(Data.obs[oname]['sim_pts'])
            idx = np.argsort(np.array(Data.sim_order))[Data.obs[oname]['sim_pts']-1] + 1                    
            tmp = np.genfromtxt(Data.obs[oname]['sim_file'],delimiter='\t',
                                skip_header=hskip,unpack=True)[idx]*Data.obs[oname]['conv']


            # Shave off the transient part (if any)
            if Config.trimB > 1:
                tmp = tmp[:,Config.trimB-1:Config.trimB-1+Config.trimL]

            valid_idx = np.fromfile('/home/wusongj/dmc/obs/'+oname+'_mask.bin', dtype=np.int64)
            tmp = tmp.flatten()[valid_idx]
            # export multiple points, modified by Songjun 
            #with open(Data.obs[oname]['sim_hist'],'a') as f_out:
            #    f_out.write(str(it+1)+','+','.join([str(j) for j in list(tmp)])+'\n')
            with open(Data.obs[oname]['sim_hist'],'a') as f_out:
                tmp.tofile(f_out)

        # Time series
        if Data.obs[oname]['type']=='Ts':
            hskip = Data.nts+3
            # export multiple points, modified by Songjun    
            #idx = np.argsort(np.array(Data.sim_order))[Data.obs[oname]['sim_pts']-1]+1
            #idx = np.array(Data.obs[oname]['sim_pts'])
            idx = np.argsort(np.array(Data.sim_order))[Data.obs[oname]['sim_pts']-1] + 1                    
            tmp = np.genfromtxt(Data.obs[oname]['sim_file'],delimiter='\t',
                                skip_header=hskip,unpack=True)[idx]*Data.obs[oname]['conv']


            # Shave off the transient part (if any)
            if Config.trimB > 1:
                tmp = tmp[:,Config.trimB-1:Config.trimB-1+Config.trimL]
            # export multiple points, modified by Songjun 
            #with open(Data.obs[oname]['sim_hist'],'a') as f_out:
            #    f_out.write(str(it+1)+','+','.join([str(j) for j in list(tmp)])+'\n')
            with open(Data.obs[oname]['sim_hist'],'a') as f_out:
                tmp.tofile(f_out)
          
        if Data.obs[oname]['type']=='map':   # renewed by Songjun, lets not use netCDF but PCraster
            f_out = open(Data.obs[oname]['sim_hist'],'a')
            for kk in range(len(startList)):
                for it in np.arange(startList[kk], endList[kk], 1):
                    mapName = Data.obs[oname]['sim_file']+str(it+1)
                    if it==startList[kk]:
                        arr  = pcr2numpy(readmap(mapName), np.nan)
                    else:
                        arr += pcr2numpy(readmap(mapName), np.nan)
                arr /= (endList[kk] - startList[kk])
                arr = arr.astype(np.float64)   
                arr.tofile(f_out)         # shape(iteration, months, row, col)
            f_out.close()



        if Data.obs[oname]['type']=='mapTs':   # renewed by Songjun, lets not use netCDF but PCraster
            arr = np.array([])
            for Ts in Data.obs[oname]['Ts']:
                mapName = Data.obs[oname]['sim_file']+str(Ts)
                arr = np.append(arr, pcr2numpy(readmap(mapName), np.nan))

            with open(Data.obs[oname]['sim_hist'],'a') as f_out:
                arr.tofile(f_out)         # shape(iteration, year, month, row, col)
                
            

# ----------------------------------------------------------------------------
# -- Restart: trim outputs files to match the specified restart iteration 

def restart(Config, Opti, Data):

    # -- Get the last iteration that worked
        
    # File names for grouped simulations
    for oname in Data.names:
        Data.obs[oname]['sim_hist'] = Config.PATH_OUT+'/'+oname+'_all.tab'

    # Read one and get the last iteration written
    # (take one before, just to be sure writing was not ongoing there when it stopped)
    idx = np.genfromtxt(Data.obs[Data.names[0]]['sim_hist'],delimiter=',',skip_header=1,
                        unpack=True,usecols=(0))
    
    Config.itres = int(idx[::-1][0])      # Removed the "-1" since want to re-run whichever was the last run being written / completed
    # In case some iterations failed the max number of rows to read is smaller than Config.itres-1!!
    
    mxRow = len(idx)-1      # Still keep -1 here since want to discard the last run, since this is being re-run
    #print Config.itres-1
    #print mxRow


    # Some cleaning of the outputs
    for oname in Data.names:
        # -- Read grouped simulations
        tmp = np.genfromtxt(Data.obs[oname]['sim_hist'],delimiter=',',skip_header=1,
                            max_rows=mxRow)[::,1::]
        
        #print tmp.shape
        # -- Take out the iterations from/beyond restarting point
        #tmp = tmp[::,1::]
        #print Config.itres
        #print Config.trim
        #tmp = tmp[0:Config.itres-1,Config.trim+1::] ## TEMPORARY!!
        #print tmp.shape
        # -- Overwrite
        # Header
        with open(Data.obs[oname]['sim_hist'],'w') as f_out:
            f_out.write('Sample,'+','.join([str(i+1) for i in range(len(tmp[0]))])+'\n')
        # Iterations before starting point
        with open(Data.obs[oname]['sim_hist'],'a') as f_out:
            for i in range(mxRow):
                f_out.write(str(idx[i])+','+','.join([str(j) for j in list(tmp[i])])+'\n')

    # -- Some cleaning for parameters
    Config.f_par = Config.PATH_OUT+'/Parameters.txt'
    # Read grouped simulations
    tmp = np.genfromtxt(Config.f_par,delimiter=',',skip_header=1,max_rows=mxRow)[::,1::]
    # Take out the iteration from/beyond the restarting point
    ###tmp = tmp[0:Config.itres,1::]
    
    # Write
    with open(Config.f_par,'w') as f_out:
        f_out.write('Iteration,'+','.join(Opti.names)+'\n')
    with open(Config.f_par,'a') as f_out:
        for i in range(mxRow):
            f_out.write(str(idx[i])+','+','.join([str(x) for x in list(tmp[i])])+'\n')


# ----------------------------------------------------------------------------
# -- Spinup + clean up/edit after

def spinup(Config):

    # Run spinup (1Y-run)
    os.system(Config.spin_ech2o+' > ech2o_spin.log')

    # Check if it ran properly
    # 1. Check it ran
    if len(glob.glob(Config.PATH_EXEC+'/BasinSummary.txt')) == 0:
        sys.exit("Something went wrong in the spinup, BasinSummary.txt is missing...")
    # 2. Check it ran until the end
    # get the last day number
    if Config.leap == 1:
        lspin = 366*Config.spinup
    else:
        lspin = 365*Config.spinup
    #format for ECH2O map outputs
    if lspin < 1000:
        espin = '0.'+format(lspin,'03')
    if lspin>=1000 :
        espin = '1.'+format(lspin-1000,'03')
    if lspin>=2000 :
        espin = '2.'+format(lspin-2000,'03')
    if lspin>=3000 :
        espin = '3.'+format(lspin-3000,'03')
    if lspin>=4000 :
        espin = '4.'+format(lspin-4000,'03')
    if lspin>=5000 :
        espin = '4.'+format(lspin-5000,'03')
    if lspin>=6000 :
        espin = '4.'+format(lspin-6000,'03')
    if lspin>=7000 :
        espin = '4.'+format(lspin-7000,'03')

    if len(glob.glob(Config.PATH_EXEC+'/Q000000'+espin)) == 0:
        sys.exit("Something went wrong in the spinup, the last-step-map Q000000"+espin+" is missing...")

    # Copy last-step maps to Spatial directory
    os.system('cp Q000000'+espin+' '+Config.PATH_SPA+'/Q.map')      # Initial streamflow
    os.system('cp SWE0000'+espin+' '+Config.PATH_SPA+'/SWE.map')    # Initial snow water equiv
    os.system('cp SWC1_00'+espin+' '+Config.PATH_SPA+'/SWC.L1.map') # Initial soil moisture L1
    os.system('cp SWC2_00'+espin+' '+Config.PATH_SPA+'/SWC.L2.map') # Initial soil moisture L2
    os.system('cp SWC3_00'+espin+' '+Config.PATH_SPA+'/SWC.L3.map') # Initial Soil moisture L3
    os.system('cp Ts00000'+espin+' '+Config.PATH_SPA+'/Ts.map')     # Initial soil temperature

    if Config.isTrck == 1:
        # -- For initial isotopes, use measurement-derived knowledge as initial value
        # Snowpack, approximate generic value, may not have huge influence
        os.system('cp '+Config.PATH_SPA+'/dD_snowpack.map '+Config.PATH_SPA+'/dD.snowpack.map') 
        os.system('cp '+Config.PATH_SPA+'/d18O_snowpack.map '+Config.PATH_SPA+'/d18O.snowpack.map') 
        # Stream : extrapolated outlet value
        os.system('cp '+Config.PATH_SPA+'/dD.stream.20130221.map '+Config.PATH_SPA+'/dD.surface.map')
        os.system('cp '+Config.PATH_SPA+'/d18O.stream.20130221.map '+Config.PATH_SPA+'/d18O.surface.map')
        # Soil : Josie's measurements (peat/gley/podzol)
        # Podzol extrapolated to ranker
        # L3 is taken as L2 (measurements depth)
        os.system('cp '+Config.PATH_SPA+'/dD.s10cm.20130221.map '+Config.PATH_SPA+'/dD.L1.map')
        os.system('cp '+Config.PATH_SPA+'/d18O.s10cm.20130221.map '+Config.PATH_SPA+'/d18O.L1.map')
        os.system('cp '+Config.PATH_SPA+'/dD.s20cm.20130221.map '+Config.PATH_SPA+'/dD.L2.map')
        os.system('cp '+Config.PATH_SPA+'/d18O.s20cm.20130221.map '+Config.PATH_SPA+'/d18O.L2.map')
        os.system('cp '+Config.PATH_SPA+'/dD.s40cm.20130221.map '+Config.PATH_SPA+'/dD.L3.map')
        os.system('cp '+Config.PATH_SPA+'/d18O.s40cm.20130221.map '+Config.PATH_SPA+'/d18O.L3.map')
        # Groundwater: extrapolation of Bernhard's deep wells (DW1-4)
        # (although it's not the same year, as it is very stable)
        os.system('cp '+Config.PATH_SPA+'/dD.GW.DW201603.map '+Config.PATH_SPA+'/dD.GW.map')
        os.system('cp '+Config.PATH_SPA+'/d18O.GW.DW201603.map '+Config.PATH_SPA+'/d18O.GW.map')

        # -- For initial water age, use spinup
        os.system('cp Agesnw0'+espin+' '+Config.PATH_SPA+'/Age.snowpack.map')
        os.system('cp Agesrf0'+espin+' '+Config.PATH_SPA+'/Age.surface.map')
        os.system('cp AgesL10'+espin+' '+Config.PATH_SPA+'/Age.L1.map')
        os.system('cp AgesL20'+espin+' '+Config.PATH_SPA+'/Age.L2.map')
        os.system('cp AgesL30'+espin+' '+Config.PATH_SPA+'/Age.L3.map')
        os.system('cp Agegw00'+espin+' '+Config.PATH_SPA+'/Age.GW.map') 

    # Clean up
    os.system('rm -f *.txt *.tab *'+espin)

# --------------------------------------------------------------------------------------------
# -- Creation of the Morris trajectories
# Follows methodology and recommandations of
# Sohier, Farges, and Piet-Lahanier (2014), Improvements of the representativity of the Morris
# method for air-launch-to-orbit separation, Proc. 19th congress of IFAC.

def morris_trajs(Config, Opti):

    # Number of trajectories -> even number determined form ncpu
    # Opti.nr = np.round(int(Config.ncpu))
    # Number of levels -> equal nr
    Opti.nlev = copy.copy(Opti.nr)

    vals={}
    
    # Step: plus-minus 0.5 of the range of each parameter
    Opti.step = np.zeros((Opti.nvar),np.float32)

    # Construct B* matrix, for each repetition
    Opti.Bnorm = np.zeros((Opti.nvar,Opti.nvar+1,Opti.nr),np.float32) # the normalized
    Opti.Bstar = np.zeros((Opti.nvar,Opti.nvar+1,Opti.nr),np.float32) # the final used in runs

    # Starting point: latin hypercube sampling, maximizing 'distance' between point, and centered
    Opti.Bnorm[:,0,:] = np.transpose(lhs(Opti.nvar,samples=Opti.nr,criterion='cm'))

    # Construct samples
    for ir in range(Opti.nr):
        for iv in range(Opti.nvar):
            # Mode 1 : trajectories
            # Mode 2 : radial points
            # In both cases the other of one-at-a-time change is fixed: 1st param changes, 
            # then 2nd param, etc.
            # the randomness is assured by the initial LHS + the fixed step of +-0.5
            
            if(Config.MSspace=='trajectory'):
                # copy previous location
                Opti.Bnorm[:,iv+1,ir] = copy.copy(Opti.Bnorm[:,iv,ir])
            elif(Config.MSspace=='radial'):
                # alway start from initial
                Opti.Bnorm[:,iv+1,ir] = copy.copy(Opti.Bnorm[:,0,ir])
            else:
                sys.exit('Wrong option for the MS parameter space definition !')

            # Successive changes with + (resp. -) 0.5, depending on if they are 
            # above (resp. below) the mid-interval
            if Opti.Bnorm[iv,iv+1,ir] < 0.5:
                Opti.Bnorm[iv,iv+1,ir] += 0.5
            elif Opti.Bnorm[iv,iv+1,ir] > 0.5:
                Opti.Bnorm[iv,iv+1,ir] -= 0.5
            else:
                Opti.Bnorm[iv,iv+1,ir] = 0.2

            # print(Opti.Bnorm[iv,iv,ir],Opti.Bnorm[iv,iv+1,ir])
                
            # Check for error
            if Opti.Bnorm[iv,iv+1,ir] > 1-1/(2*Opti.nlev) or Opti.Bnorm[iv,iv+1,ir] < 1/(2*Opti.nlev):
                print('Error in the incrementation of the parameter '+Opti.names[iv])
                print(Opti.Bnorm[iv,iv,ir],Opti.Bnorm[iv,iv+1,ir])
                sys.exit()
                
    # Construct the actual Bstar, with non-normalized values
    for iv in range(Opti.nvar):
        if Opti.log[iv]==1:                
            Opti.Bstar[iv,:,:] = 10**(Opti.Bnorm[iv,:,:]*np.log10(Opti.max[iv]/Opti.min[iv])+np.log10(Opti.min[iv]))
            Opti.step[iv] = 0.5 * (np.log10(Opti.max[iv])-np.log10(Opti.min[iv]))
        else:
            Opti.Bstar[iv,:,:] = Opti.Bnorm[iv,:,:]*(Opti.max[iv]-Opti.min[iv]) + Opti.min[iv]
            Opti.step[iv] = 0.5 * (Opti.max[iv]*Opti.min[iv])

    # Check if outputs directory exists
    if len(glob.glob(Config.PATH_TRAJ)) == 0: os.system('mkdir '+ Config.PATH_TRAJ)
            
    # Write on file giving parameters range, log...(for later plots)    
    f_out = Config.FILE_TRAJ+'.Parameters_char.txt'
    with open(f_out,'w') as fw:
        # fw.write('Sample,'+','.join(Opti.names)+'\n')
        fw.write('Names,'+','.join(Opti.names)+'\n')
        fw.write('Min,'+','.join([str(Opti.min[x]) for x in range(Opti.nvar)])+'\n')
        fw.write('Max,'+','.join([str(Opti.max[x]) for x in range(Opti.nvar)])+'\n')
        fw.write('Log,'+','.join([str(Opti.log[x]) for x in range(Opti.nvar)])+'\n')
        fw.write('Step,'+','.join([str(Opti.step[x]) for x in range(Opti.nvar)])+'\n')
        fw.write('StepN,'+','.join([str(Opti.stepN[x]) for x in range(Opti.nvar)])+'\n')

    # Write Bstar for each trajectory
    for ir in range(Opti.nr):
        trajnb = str(ir+1)#'%02i' % int(ir+1)
        # print trajnb
        with open(Config.FILE_TRAJ+'.Bstar_traj'+trajnb+'.txt','w') as fw:
            csv_writer = csv.writer(fw)
#            print(Opti.names)
            csv_writer.writerow(Opti.names)
            for irun in range(Opti.nvar+1):
                csv_writer.writerow(Opti.Bstar[:,irun,ir])
        exit

# --------------------------------------------------------------------------------------------
# -- Calculation of elementary effects for Morris Sensitivity Analysis
# Applied to one trajectory / radial points with same origin
# Since time series are the outputs, metrics for d(outputs) are bias and RMSE

# !! Pandas unavailable so using numpy (unhandy) numpy for I/O+data-manip
def morris_ee(Config, Data, Opti):

    firstObs = 0
    numObs = 0
    outObs = []

    for oname in Data.names:
        
        # Only look into time series
        if Data.obs[oname]['type'] != 'map':

            # Read file
            f_in = Data.obs[oname]['sim_hist']
            # df_sim = pd.read_csv(f_in,skiprows=1)

            # # Diff between two sims
            # df_diff = df_sim.set_index('Sample').diff().loc[2:Opti.nvar+1,]

            # # Get bias
            # bias = df_diff.mean(axis=1)
            # # Get RMSE
            # RMSE = np.sqrt((df.diff**2).mean(axis=1))
            # # Get corresponding elementary effect
            # bias_ee = bias / Opti.stepN
            # RMSE_ee = RMSE / Opti.stepN
            # # Write the files
            # bias_ee.index = Opti.names
            # RMSE_ee.index = Opti.names

            sim = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True)[1:Config.trimL+1,:]

            # Take into account accumulated fluxes
            if Data.obs[oname]['type']=='Total' and Data.obs[oname]['sim_pts'] in [1,11,12,13,14,15,16,17,18,19,20]:
                sim = np.diff(sim,axis=0)
            # Diff between sims
            if(Config.MSspace=='trajectory'):
                simd = np.diff(sim)
            elif(Config.MSspace=='radial'):
                simd = sim[:,1::] - sim[:,0][...,None]
            print(Opti.nvar)
            # Elementary effect (keep direction of param change)
            bias_ee = np.zeros((Opti.nvar),np.float32)*np.nan
            RMSE_ee = np.zeros((Opti.nvar),np.float32)*np.nan
            for i in range(Opti.nvar):
                bias_ee[i] = np.mean(simd[:,i]) / Opti.dx[i,i]
                RMSE_ee[i] = np.sqrt(np.mean(simd[:,i]**2)) / Opti.dx[i,i]
            # bias_ee = bias / Opti.BnormstepN
            # RMSE_ee = RMSE / Opti.stepN

            # Add the overall data frame
            if(firstObs==0):
                bias_ee_tot = bias_ee[...,None] # Creates a (..,1) dimension
                RMSE_ee_tot = RMSE_ee[...,None] # Creates a (..,1) dimension
                # bias_ee_tot = pd.DataFrame(bias_ee).assign(oname=bias_ee)
                # RMSE_ee_tot = pd.DataFrame(RMSE_ee).assign(oname=RMSE_ee)
            else:
                bias_ee_tot = np.append(bias_ee_tot,bias_ee[...,None],1)
                RMSE_ee_tot = np.append(RMSE_ee_tot,RMSE_ee[...,None],1)
                # bias_ee_tot = bias_ee_tot.assign(oname=bias_ee)
                # RMSE_ee_tot = RMSE_ee_tot.assign(oname=RMSE_ee)

            # Update
            firstObs = 1
            # Increment number of obs actually evaluated
            numObs += 1
            # Append name of obs actually evaluated
            outObs = outObs + [oname]

    # # Ugly: drop the first column (named 0) that had to be left (duplicate of 2nd column)
    # bias_ee_tot.drop(0,axis=1,inplace=True)
    # RMSE_ee_tot.drop(0,axis=1,inplace=True)
    # Write outputs ------------------------------------------------------------------------

    # Check if directory exists
    if len(glob.glob(Config.PATH_EE)) == 0: os.system('mkdir '+ Config.PATH_EE)

    if(Config.MSspace=='trajectory'):
        # bias_ee_tot.to_csv(Config.FILE_EE+'.EE.Traj'+trajnb+'.bias.txt')
        # RMSE_ee_tot.to_csv(Config.FILE_EE+'.EE.Traj'+trajnb+'.RMSE.txt')
        with open(Config.FILE_EE+'.EE.Traj'+Config.numsim+'.bias.txt','w') as f_out:
            f_out.write('Parameter'+','+','.join([outObs[j] for j in range(numObs)])+'\n')
            for i in range(Opti.nvar):
                f_out.write(Opti.names[i]+','+','.join([str(bias_ee_tot[i,j]) for j in range(numObs)])+'\n')
        exit
        with open(Config.FILE_EE+'.EE.Traj'+Config.numsim+'.RMSE.txt','w') as f_out:
            f_out.write('Parameter'+','+','.join([outObs[j] for j in range(numObs)])+'\n')
            for i in range(Opti.nvar):
                f_out.write(Opti.names[i]+','+','.join([str(RMSE_ee_tot[i,j]) for j in range(numObs)])+'\n')
        exit

    if(Config.MSspace=='radial'):
        # bias_ee_tot.to_csv(Config.FILE_EE+'.EE.RadP'+trajnb+'.bias.txt')
        # RMSE_ee_tot.to_csv(Config.FILE_EE+'.EE.RadP'+trajnb+'.RMSE.txt')
        with open(Config.FILE_EE+'.EE.RadP'+Config.numsim+'.bias.txt','w') as f_out:
            f_out.write('Parameter'+','+','.join([outObs[j] for j in range(numObs)])+'\n')
            for i in range(Opti.nvar):
                f_out.write(Opti.names[i]+','+','.join([str(bias_ee_tot[i,j]) for j in range(numObs)])+'\n')
        exit
        with open(Config.FILE_EE+'.EE.RadP'+Config.numsim+'.RMSE.txt','w') as f_out:
            f_out.write('Parameter'+','+','.join([outObs[j] for j in range(numObs)])+'\n')
            for i in range(Opti.nvar):
                f_out.write(Opti.names[i]+','+','.join([str(RMSE_ee_tot[i,j]) for j in range(numObs)])+'\n')
        exit
    
