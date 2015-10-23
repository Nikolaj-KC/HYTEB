# -*- coding: utf-8 -*-
#@author: Nikolaj

import os
import numpy as np
import math as mat
import load_write
import subprocess
import shutil

# *****************************************************************************
# *****************************************************************************
class hydro_module:
    """
    Class containing hydrology
    """
    
    
    def __init__(self,model_dis ='model_dis.ini', hydro_ini = 'hydro_ini.ini',inversion_ini='inversion_setup.ini'):

        self.write = load_write.write_file()
        
        self.load = load_write.load_file()
        self.load.ini_model_dis(model_dis)
        self.load.ini_hydro(hydro_ini)
        self.load.ini_inversion(inversion_ini)

#==============================================================================
    def head2folder(self,model_num,modflow_nam,well):
#==============================================================================
        """
        Copy true hydraulic head measurements and head prediction measurements 
        to folder.
        
        Parameters 
        ----------
        model_num :         int
            Model realization number
        modflow_nam :       str
            Name of modlfow name file
        well :              str
            on  --> run modflow with active well, see modflow-nam file \n
            off --> run modflow without active well, see modflow-nam file
            
       
        """
        

        print 'head2folder      ',
        
        
        try:
            subprocess.call('mf2k.exe '+modflow_nam+'.nam')
        except:
            print 'WARNING IN HEAD2FODLER'
            print 'MODFLOW run fail'
        
        if well == 'off':

        #-------------------------------
        #--- run GWM WITHOUT PUMPING ---        
        #-------------------------------
                
#        try:
#            if run_type == 'true':
#                subprocess.call('mf2k.exe '+self.load.modflow_file_name+'_true.nam')               
#            if run_type == 'calib':
#                subprocess.call('mf2k.exe '+self.load.modflow_file_name+'_calib.nam')
#            else:
#                print 'WARNING, wrong input value for run_type '
#                
#        except:
#            print 'Warning in run_prediction'
#            print 'MODFLOW fail in run WITHOUT pumping'
        
            
    
            head_obs_file = [self.load.head_obs_points,self.load.head_pred_points]
            
            # --- write in files for mod2smp ---
            for i in head_obs_file:
                output_file         = i.split('.')
                try:
                    
                    output_file_name    = 'mod2smp_'+output_file[0]+modflow_nam+'.in'
                    self.write.mod2smp_in(output_file_name,self.load.modflow_file_name,self.load.model_spc,i)
                    print '     write:               '+output_file_name
                    
                    print '     run mod2smp:        ',
                    myinput = open(output_file_name,'r')
                    tmp     = open('null.dat','w')
                    subprocess.call('mod2smp.exe',stdin=myinput,stdout=tmp)
                    myinput.close()
                    tmp.close()  
                    
                    src = output_file[0]
                    dst = self.load.path_results+'/'+src
                    shutil.copyfile(src+'.smp',dst+'.smp'+'.'+str(model_num)) 
                                   
                    print 'done'  
                except:
                    print 'head2folder FAIL:   well='+well+'  ',output_file
                            
        # --------------------------------------------------------------------
        # --- run GWM WITH PUMPING ---        
        # --------------------------------------------------------------------              
        elif well == 'on':
#            try:
#                if run_type == 'true':
#                    subprocess.call('mf2k.exe '+self.load.modflow_file_name+'_pump.nam')               
#                if run_type == 'calib':
#                    subprocess.call('mf2k.exe '+self.load.modflow_file_name+'_calib_pump.nam')
#                else:
#                    print 'WARNING, wrong input value for run_type '
#            except:
#                print 'Warning in run_prediction'
#                print 'MODFLOW fail in run WITH pumping'
#                
                   
              
            head_obs_file = [self.load.head_obs_points,self.load.head_pred_points]
            
            # --- write in files for mod2smp ---
            for i in head_obs_file:
                print i
                output_file         = i.split('.')
                try:
                    output_file         = i.split('.')
                    
                    output_file_name    = 'mod2smp_'+output_file[0]+modflow_nam+'.in'
                    self.write.mod2smp_in(output_file_name,self.load.modflow_file_name,self.load.model_spc,i)
                    print '     write:               '+output_file_name
                    
                    print '     run mod2smp:        ',
                    myinput = open(output_file_name,'r')
                    tmp     = open('null.dat','w')
                    subprocess.call('mod2smp.exe',stdin=myinput,stdout=tmp)
                    myinput.close()
                    tmp.close()  
                    
                    src = output_file[0]
                    dst = self.load.path_results+'/'+src+'_pump'
                    shutil.copyfile(src+'.smp',dst+'.smp'+'.'+str(model_num)) 
                                   
                    print 'done'  
                except:
                    print 'head2folder FAIL:   ',output_file
                    
                
                
#==============================================================================
    def ibound(self,x_top,x_bot,valley_north,fixed_head):
#==============================================================================
        """
        Ibound files for modflow

        Parameters
        ----------
        x_top :             float
            Width of the burried valley top 
        x_bot :             float
            Width of the burried valley top
        valley_north :      float
            Length of valley from southern bondary
        fixed_head :        int
            Row number where head is constant in modflow simulatons. 
            
        Returns
        -------
        #.inf
          
            
        """
        
        print 'ibound for modflow      ',        
        
        # --- node position for southern and northern extension of valley --- 
        y_south = self.load.ny-1
        y_north = (self.load.ny*self.load.dx-valley_north)/self.load.dx+self.load.n_cap
        
        # --- depth of valley ---
        depth   = self.load.dz*(self.load.n_layer-self.load.n_cap)   
        
        # --- slope of valley --- 
        angel   = mat.degrees(mat.atan(depth/((x_top/2.0)-(x_bot/2.0))))
        print '\t slope of valley sides:   ', np.round(angel,2),
        
        # --- valley going from node nx1 to nx2 in east-west direction ---
        nx1      = (self.load.dx*self.load.nx/2.0-x_top/2.0)/float(self.load.dx)
        nx2      = (self.load.dx*self.load.nx/2.0+x_top/2.0)/float(self.load.dx)

        # --- allocate numpy array, and fill array  ---
        ibound   = 99999*np.ones((self.load.ny,self.load.nx,self.load.n_layer))
        for layer in range(self.load.n_layer-1,-1,-1):    
            # --- for capping part of groundwater model ---            
            if layer < self.load.n_cap:
                for i in range(self.load.ny):
                    for j in range(self.load.nx):
                        if i == fixed_head-1:
                            ibound[i,j,layer] = -1
                        else:
                            ibound[i,j,layer] = 1 
            
            # for valley part of groundwater model ---
            if layer >= self.load.n_cap:
                cells = ((layer-5)*self.load.dz/mat.tan(mat.radians(angel)))/self.load.dx
                cell = np.ceil(cells)
                for i in range(self.load.ny):
                    for j in range(self.load.nx):
                        if i > y_north+layer and i < y_south-layer:
                            if j > nx1+cell and j < nx2-cell:
                                ibound[i,j,layer] = 1
                            else:
                                ibound[i,j,layer] = 0
                        else:
                            ibound[i,j,layer] = 0

        # --- write ibound file. One for each layer ---
        for layer in range(self.load.n_layer):
            fid1 = open(self.load.ibound+str(layer+1)+'.'+self.load.dtype_ibound,'w')
            for i in range(self.load.ny):
                for j in range(self.load.nx):
                    fid1.write(str(int(ibound[i,j,layer]))+'\t')
                fid1.write('\n')
            fid1.close()
        
        print 'done'
        
#==============================================================================
    def int2zone(self,sufix=''):        
#==============================================================================
        """
        write modflow zone file from integer file 
        
        Parameters
        ----------
        sufix : str
                optional; putting extension to standart modflow project name
        
        Returns
        -------
        "modflow project name".zns        
        
        """
        
        # --- allocate numpy array ---
        data = np.zeros(((self.load.ny,self.load.nx,self.load.n_layer)),dtype=int)

        # --- load zone array ---
        for layer in range(self.load.n_layer):
            fid1 = open('zone_int_array_'+str(layer+1)+'.inf','r')
            data[:,:,layer] = np.array(fid1.read().split()).reshape(self.load.ny,self.load.nx)
            fid1.close()
        
        # --- write modflow zone file ---
        fid1 = open(self.load.modflow_file_name+sufix+'.zns','w')
        fid1.write(str(self.load.n_layer)+'\n')
        for layer in range(self.load.n_layer):
            fid1.write('ZLAY_'+str(layer+1)+'\n')
            fid1.write('INTERNAL 1 (FREE) -1'+'\n')
            for i in range(self.load.ny):
                for j in range(self.load.nx):
                    fid1.write(str(data[i,j,layer])+' ')
                fid1.write('\n')
        fid1.close()

#==============================================================================
    def integer_array_file(self):         
#==============================================================================
        """
        Write intger array file         
        have to run prior to pilot_point_file

        Returns
        -------
        zone_mat :          array_like, float
            
            zone_int_array_L(layer).inf
        
        
        """
        

        if self.load.parametrization == -1:
            # --- one zone for modflow model ----
            zone_mat = np.ones((self.load.ny,self.load.nx,self.load.n_layer),dtype=int)
            
        elif self.load.parametrization == -2 or self.load.parametrization == -3:
            # --- load zone array from modflow in to matrice ---
            zone_mat = self.load.zns2mat(self.load.modflow_zone_file,self.load)            
       
        else:
            'Error in interger array!'

        # --- write integer array file for each layer ---
        self.write.integer_array_file(self.load.integer_array_name,zone_mat)
        
        return zone_mat
        
        

#==============================================================================
    def model2mat(self,model_type='true'):        
#==============================================================================
        """
        load 3D model into numpy array
        
        Parameters
        ----------
        model_type :        logical
            true  --> load "true" model,  (see HYTEB-ini-file)
            calib --> load "calib" model, (see HYTEB-ini-file)
        
        Returns
        -------
        mat :               numpy.array
        
        """
        
        

        mat = np.zeros((self.load.ny,self.load.nx,self.load.n_layer),dtype=float)
        for i in range(self.load.n_layer):
            if model_type == 'true':
                fname = self.load.model_name+str(i+1)+'.'+self.load.dtype_model
            elif model_type == 'calib':
                fname = self.load.model_name_calib+str(i+1)+'.'+self.load.dtype_model
            else:
                print 'Wrong input for MODEL_TYPE!!!'
            fid1 = open(fname,'r')
            mat[:,:,i] = np.array(fid1.read().split(),dtype=float).reshape(self.load.ny,self.load.nx)
            fid1.close()
            
        return mat                
                    
#==============================================================================
    def modflow_adj(self,sufix='',operation='obsall',iuparobs=49,iadjsu=43,
                    iadjfm=12,iadjscl=0,iadjpc=0,iuobsmap=0,
                    iusolv=48,iadjxdu=44,iaform=3,discut=.000001,
                    pvalueout=1,ipdisp=1,
                    resetp='noreset',parnam='hk_',ln=0):
#==============================================================================
        """
        Input for adjoint based sensitivity package
        
        Parameters
        ----------
        sufix :             str
            sufix to modflow nam-file, see HYTEB-ini-file
        operation :         str
            default('obsall'), see modflow_adjoint doc for more information
        iuparobs :          int
            defalut(49), see modflow_adjoint doc for more information 
        iadjsu :            int
            default(43), see modflow_adjoint doc for more information
        iadjfm              int
            deafault(12), see modflow_adjoint doc for more information
        iadjscl             int
            dafault(0), see modflow_adjoint doc for more information
        iadjpc              int
            default(0), see modflow_adjoint doc for more information
        iuobsmap            int
            default(0), see modflow_adjoint doc for more information
        iusolv              int
            default(48), see modflow_adjoint doc for more information
        iadjxdu             int
            default(44), see modflow_adjoint doc for more information
        iaform              int
            dafualt(3), see modflow_adjoint doc for more information
        discut              float
            default(.000001), see modflow_adjoint doc for more information
        pvalueout           int
            dafualt(1), see modflow_adjoint doc for more information
        ipdisp              int
            default(1), see modflow_adjoint doc for more information
        resetp              str
            dafault('noreset'), see modflow_adjoint doc for more information
        parnam              str
            default('hk\_'), see modflow_adjoint doc for more information
        ln                  int
            default(0), see modflow_adjoint doc for more information

        
        """
        
        print 'modflow ADJ file:          ',
        
        npe     = self.load.n_layer
        npardis = self.load.n_layer       
        
        fid1 = open(self.load.modflow_file_name+'.adj'+sufix,'w')
        fid1.write('# Input for adjoint based sensitivity package'+'\n')
        fid1.write(' '+operation+'  '+str(iuparobs)+'\n')
        fid1.write(' '+str(iadjsu)+'  '+str(iadjfm)+'  '+str(iadjscl)+'  '+str(iadjpc)+'\n')
        fid1.write(' '+str(iuobsmap)+'\n')
        fid1.write(' '+str(iusolv)+'\n')
        fid1.write(' '+str(npe)+'  '+str(npardis)+'\n')
        fid1.write(' '+str(iadjxdu)+'  '+str(iaform)+' '+str(discut)+'  '+str(pvalueout)+'\n')
        for i in range(npe):           
            fid1.write(' '+str(ipdisp+i))
        fid1.write('\n')
        fid1.write(' '+str(resetp)+'\n')
        for i in range(npe):
            fid1.write(' '+str(parnam)+str(i+1)+' '+str(ln)+'\n')
        fid1.close()

        
        print 'done'
#==============================================================================
    def modflow_bas(self,sufix='',external=int(0)):
#==============================================================================            
        """
        making of MODFLOW BAS file

        Parameters
        ----------
        sufix :             str, optional
            sufix to modflow nam-file, see HYTEB-ini-file
        external:           int
            external > 0, --> specify modflow.hds file as initial head valus
        

        """
        print 'make MODFLOW BAS file:       ',

        # --- write modflow bas file ---        
        fid1 = open(self.load.modflow_file_name+sufix+'.bas','w')
        fid1.write('# MODFLOW2000 Basic Package'+'\n')
        fid1.write('FREE'+'\n')
        for layer in range(self.load.n_layer):
            fid1.write("OPEN/CLOSE '"+self.load.ibound+str(layer+1)+'.'+self.load.dtype_ibound+"'  1    "+"'(FREE)'"+'  -1    ibound layer '+str(layer+1)+'\n')
        fid1.write(' 1e30'+'\n')
        
        # --- homogen head field as inital starting point ---
        if external == 0:
            for layer in range(self.load.n_layer):
                fid1.write('         0 0.000e+00(10e12.4)                   -1'+'\n')
        # --- use external head file as intial starting point ---
        if external > 0:
            for layer in range(self.load.n_layer):
                fid1.write('EXTERNAL	'+str(external)+'	1	(BINARY)	-1'+'\n')
        fid1.close()  
        
        print 'done'        

#==============================================================================
    def modflow_dis(self):
#==============================================================================
        """
        write MODFLOW2000 Discretization File
        
        
        """
        print 'make MODFLOW DIS file:       ',
        
        # --- write dis file   ---      
        fid1 = open(self.load.modflow_file_name+'.dis','w')
        fid1.write('# MODFLOW2000 Discretization File'+'\n')
        fid1.write(' '+str(self.load.n_layer)+'  '+str(self.load.ny)+'  '+str(self.load.nx)+'  1  1  2'+'\n')
        fid1.write(' ')
        for layer in range(self.load.n_layer):
            fid1.write(str(0)+' ')
        fid1.write('\n')
        fid1.write('         0    '+str(format(self.load.dx,'.3f'))+'(10E12.4)                    0'+'\n')
        fid1.write('         0    '+str(format(self.load.dy,'.3f'))+'(10E12.4)                    0'+'\n')
        fid1.write('         00.0000e+00(10e12.4)                   -1'+'\n')
        for layer in range(self.load.n_layer):
            fid1.write('         0   '+str(format(-self.load.dz*(layer+1),'.3f'))+'(10e12.4)                   -1'+'\n')
        fid1.write('  1.000000  1  1.000000  SS')
        fid1.close()
        
        print 'done'

#==============================================================================
    def modflow_gmg(self,rclose=1.0, iiter=100,hclose=1e-3,mxiter=200,damp=0.3,
                    iadamp=1,ioutgmg=1,ism=0,isc=1,relax=1):
#==============================================================================
        """
        modflow solver module GMG
        
        Parameters
        ----------
        rclose :            float, optional
            see MODFLOW-2000 manual for more information
        iiter :             int, optional
            see MODFLOW-2000 manual for more information
        hclose :            float, optional
            see MODFLOW-2000 manual for more information
        mxiter :            int, optional
            see MODFLOW-2000 manual for more information
        damp :              float, optional
            see MODFLOW-2000 manual for more information
        iadamp :            int, optional
            see MODFLOW-2000 manual for more information
        ioutgmg :           int, optional
            see MODFLOW-2000 manual for more information
        ism :               int, optional
            see MODFLOW-2000 manual for more information
        isc :               int, optional
            see MODFLOW-2000 manual for more information
        relax :             int, optional
            see MODFLOW-2000 manual for more information
        
        """
        
        print 'MODFLOW gmg file:              ',
        fid1 = open(self.load.modflow_file_name+'.gmg','w')
        fid1.write(str(rclose)+'  '+str(iiter)+'  '+str(hclose)+'  '+str(mxiter)+'     RCLOSE IITER HCLOSE MXITER \n')
        fid1.write(str(damp)+'  '+str(iadamp)+'  '+str(ioutgmg)+'            DAMP IADAMP IOUTGMG \n')
        fid1.write(str(ism)+'  '+str(isc)+'                          ISM ISC \n')
        fid1.write(str(relax)+'                           RELAX')
        fid1.close()
        
        print 'done'
        

#==============================================================================
    def modflow_kfield(self,data,sufix='',log=None):        
#==============================================================================
        """
        Write mulitplier K-fields for modflow adjoint
        
        Parameters
        ----------
        data :              numpy.array(row,col)
            2D array with hydraulic conductivity field
        sufix :             str, optional
            sufix to modflow nam-file, see HYTEB-ini-file
        
        """
        

        [ix,iy,iz] = np.shape(data)

        for z in range(iz):
            fid1 = open(self.load.modflow_file_name+sufix+str(z+1)+'.dat','w')
            for i in range(ix):
                for j in range(iy): 
                    if log != None:   
                        fid1.write('%14.6e'%(np.log10(data[i,j,z])))
                    else:
                        fid1.write('%14.6e'%(data[i,j,z]))
                fid1.write('\n')
            
        fid1.close()
        
#==============================================================================
    def modflow_lpf(self,sufix='',run_type='true'):
#==============================================================================
        """
        
        MODFLOW2000 Layer Property Flow (LPF) Package
        
        Parameters
        ----------
        sufix :             str, optional
            sufix to modflow nam-file, see HYTEB-ini-file
        run_type :          str
            true    --> run modflow with true model
            calib   --> run modflow, GWM discretized by pilot points (real array
            files)
            zone    --> run modflow, GWM discretized by zones
            adjoint --> run modflow, GWM discretized by pilot poins (real array
            files)
      

        """
        print 'make MODFLOW LPF file:       ',
        print run_type,

        count = 0
        while count < 1:
                # --- True model LFP file ---
                if run_type == 'true':            
                    fid1 = open(self.load.modflow_file_name+sufix+'.lpf','w')
                    fid1.write('# MODFLOW2000 Layer Property Flow (LPF) Package'+'\n')
                    fid1.write('        50  -1.000000e+030         0'+'\n')
                    for layer in range(self.load.n_layer):
                        fid1.write(' 0 ')                
                    fid1.write('\n')
                    for layer in range(self.load.n_layer):
                        fid1.write(' 0 ')
                    fid1.write('\n ')
                    for layer in range(self.load.n_layer):
                        fid1.write(format(1,'.3e'))
                        fid1.write(str(' '))
                    fid1.write('\n')
                    for layer in range(self.load.n_layer):
                        fid1.write(' 0 ')                
                    fid1.write('\n')
                    for layer in range(self.load.n_layer):
                        fid1.write(' 0 ')
                    fid1.write('\n')
                    for layer in range(self.load.n_layer):
                        fid1.write(" OPEN/CLOSE '"+self.load.model_name+str(layer+1)+"."+self.load.dtype_model+"'  1.0  '(FREE)'  -1                 hk layer"+str(layer+1)+'\n')
                        fid1.write(" OPEN/CLOSE '"+self.load.model_name+str(layer+1)+"."+self.load.dtype_model+"'  1.0  '(FREE)'  -1                 vk layer"+str(layer+1)+'\n')
                    fid1.close()
                    
                # --- pilot point dis LPF file ---
                elif run_type == 'calib':
                    fid1 = open(self.load.modflow_file_name+sufix+'.lpf','w')
                    fid1.write('# MODFLOW2000 Layer Property Flow (LPF) Package'+'\n')
                    fid1.write('        50  -1.000000e+030         0'+'\n')
                    for layer in range(self.load.n_layer):
                        fid1.write(' 0 ')                
                    fid1.write('\n')
                    for layer in range(self.load.n_layer):
                        fid1.write(' 0 ')
                    fid1.write('\n ')
                    for layer in range(self.load.n_layer):
                        fid1.write(format(1,'.3e'))
                        fid1.write(str(' '))
                    fid1.write('\n')
                    for layer in range(self.load.n_layer):
                        fid1.write(' 0 ')                
                    fid1.write('\n')
                    for layer in range(self.load.n_layer):
                        fid1.write(' 0 ')
                    fid1.write('\n')
                    for layer in range(self.load.n_layer):
                        fid1.write(" OPEN/CLOSE '"+self.load.model_name_calib+str(layer+1)+"."+self.load.dtype_model+"'  1.0  '(FREE)'  -1                 hk layer"+str(layer+1)+'\n')
                        fid1.write(" OPEN/CLOSE '"+self.load.model_name_calib+str(layer+1)+"."+self.load.dtype_model+"'  1.0  '(FREE)'  -1                 vk layer"+str(layer+1)+'\n')
                    fid1.close()
                
                # --- zone dis LPF file ---
                elif run_type == 'zone':
                    fid1 = open(self.load.modflow_file_name+sufix+'.lpf','w')
                    fid1.write('# MODFLOW2000 Layer Property Flow (LPF) Package'+'\n')
                    fid1.write('        50  -1.000000e+030         '+str(2*self.load.n_zone_hk)+'\n')
                    for layer in range(self.load.n_layer):
                        fid1.write(' 0 ')                
                    fid1.write('\n')
                    for layer in range(self.load.n_layer):
                        fid1.write(' 0 ')
                    fid1.write('\n ')
                    for layer in range(self.load.n_layer):
                        fid1.write(format(1,'.3e'))
                        fid1.write(' ')
                    fid1.write('\n')
                    for layer in range(self.load.n_layer):
                        fid1.write(' 1 ')                
                    fid1.write('\n')
                    for layer in range(self.load.n_layer):
                        fid1.write(' 0 ')
                    fid1.write('\n')
                    
                    hz = self.load.hydro_zones('hydro_zones.ini',self.load.n_zone_hk,4)
                    
                    for i in range(self.load.n_zone_hk):
                            fid1.write(' '+str(hz[i,0])+'_'+str(hz[i,1])+' '+str(hz[i,0])+' '+str(hz[i,2])+' '+str(self.load.n_layer)+'\n')
                            for layer in range(self.load.n_layer):
                                fid1.write(' '+str(layer+1)+' NONE ZLAY_'+str(layer+1)+' '+str(hz[i,3])+'\n')
                    for i in range(self.load.n_zone_hk,2*self.load.n_zone_hk,1):
                            fid1.write(' '+str(hz[i,0])+'_'+str(hz[i,1])+' VANI '+str(1.0)+' '+str(self.load.n_layer)+'\n')
                            for layer in range(self.load.n_layer):
                                fid1.write(' '+str(layer+1)+' NONE ZLAY_'+str(layer+1)+' '+str(hz[i,3])+'\n')
                    for i in range(2*self.load.n_layer):
                            fid1.write(' '+str(-1)+'\n')
                    fid1.close()
                    
                # --- ---
                elif run_type == 'adjoint':
                    fid1 = open(self.load.modflow_file_name+sufix+'.lpf','w')
                    fid1.write('# MODFLOW2000 Layer Property Flow (LPF) Package'+'\n')
                    fid1.write('        0  -1.000000e+030         '+str(self.load.n_layer)+'\n') #1\n')
                    for layer in range(self.load.n_layer):
                        fid1.write(' 0 ')                
                    fid1.write('\n')
                    for layer in range(self.load.n_layer):
                        fid1.write(' 0 ')
                    fid1.write('\n ')
                    for layer in range(self.load.n_layer):
                        fid1.write(format(1,'.3e'))
                        fid1.write(' ')
                    fid1.write('\n')
                    for layer in range(self.load.n_layer):
                        fid1.write(' 0 ')                
                    fid1.write('\n')
                    for layer in range(self.load.n_layer):
                        fid1.write(' 0 ')
                    fid1.write('\n')
                     
                    for layer in range(self.load.n_layer):
                        fid1.write('HK_'+str(layer+1)+'  HK   1  1 \n')#+str(self.load.n_layer)+' \n')
                        fid1.write(str(layer+1)+' HK'+str(layer+1)+' ZLAY_'+str(layer+1)+' 1'+' \n') 
                    for layer in range(self.load.n_layer):
                        fid1.write('         2                                              HK of layer 1 print code'+' \n')
                        fid1.write('         0         1(20G14.0)                    1  12. VK of layer 1'+' \n')
                else:
                     print 'Wrong model_name in MODLFOW_LPF!'                             
                count += 1
                
        print 'done'
    
#==============================================================================
    def modflow_pcg(self,sufix='',MXITER=1000,ITER1=30,NPCOND=1,HCLOSE=0.001,RCLOSE=0.001,RELAX=1,NBPOL=1,IPRPCG=1,MUTPCG=0,DAMP=1):
#==============================================================================
        """
        MODFLOW PCG for the Preconditioned Conjugate-Gradient Package
        See fx MODFLOW-2000 documentation for varivable explanations
        
        Parameter
        ---------
        sufix :             str, optional
            sufix to modflow nam-file, see HYTEB-ini-file
        MXITER:             int
            See fx MODFLOW-2000 documentation for varivable explanations
        ITER1:              int
            See fx MODFLOW-2000 documentation for varivable explanations
        NPCOND:             int
            See fx MODFLOW-2000 documentation for varivable explanations
        HCLOSE:             float
            See fx MODFLOW-2000 documentation for varivable explanations
        RCLOSE:             float
            See fx MODFLOW-2000 documentation for varivable explanations
        RELAX:              int
            See fx MODFLOW-2000 documentation for varivable explanations
        NBPOL:              int
            See fx MODFLOW-2000 documentation for varivable explanations
        IPRPCG:             int
            See fx MODFLOW-2000 documentation for varivable explanations
        MUTPCG:             int
            See fx MODFLOW-2000 documentation for varivable explanations
        DAMP:               int
            See fx MODFLOW-2000 documentation for varivable explanations
        
        """
        
        print 'MODFLOW pcg file:            ',
        fid1 = open(self.load.modflow_file_name+sufix+'.pcg','w')
        fid1.write(str(MXITER)+' '+str(ITER1)+' '+str(NPCOND) +'\n')
        fid1.write(str(HCLOSE)+' '+str(RCLOSE )+' '+str(RELAX )+' '+str(NBPOL )+' '+str(IPRPCG )+' '+str(MUTPCG )+' '+str(DAMP))
        fid1.close()
        
        print 'done'
              

#==============================================================================
    def modflow_mpn(self,run_type='true',well ='on',recharge ='true',sufix=''):
#==============================================================================
        """
        Write MODPATH MPN file
         
        Parameters
        ----------
        run_type :          str
            true    --> run modflow with true model
            calib   --> run modflow, GWM discretized by pilot points (real array
            files)
            
        """
        print 'MODFLOW MPN file:           ',
        
        
        fid1 = open(self.load.modflow_file_name+sufix+'.mpn','w')
        fid1.write('main          11      valley.mpt'+'\n')
        if recharge == 'true':
            fid1.write('rch           12      valley.rch'+'\n')
        if recharge == 'calib':
            fid1.write('rch           12      valley_calib.rch'+'\n')
        if well == 'on':
            fid1.write('wel           13      valley.wel'+'\n')
        if run_type == 'calib':
            fid1.write('budget        40      valley_calib.cbb'+'\n')
            fid1.write('head(binary)  30      valley_calib.hds'+'\n')
            fid1.write('dis           15      valley.dis'+'\n')
        if run_type == 'true':
            fid1.write('budget        40      valley.cbb'+'\n')
            fid1.write('head(binary)  30      valley.hds'+'\n')
            fid1.write('dis           15      valley.dis'+'\n')
        if run_type != 'true' and run_type !='calib':
            print "WARNING, wrong input for run_type!"
        fid1.close()
        
        print 'done'

    
#==============================================================================
    def modflow_mpt(self):         
#==============================================================================
        """ 
        Write MODPATH MPT file 
        """

        print 'MODFLOW MPT file          '
        
        fid1 = open(self.load.modflow_file_name+'.mpt','w')
        fid1.write(str(self.load.nx)+' '+str(self.load.ny)+' '+str(self.load.n_layer)+' 0   0    1    0     1.0e35  1.0e35  10000    NCOL NROW NLAY NCBL IGRID NPER MAXSIZ HNOFLO HDRY NPART'+'\n')
        fid1.write('\n')
        for layer in range(self.load.n_layer):
            fid1.write(' 0')
        fid1.write('\n')                
        for layer in range(self.load.n_layer):
            fid1.write(" OPEN/CLOSE '"+self.load.ibound+str(layer+1)+"."+self.load.dtype_ibound+"'  1  '(FREE)'  -1                ibound layer "+str(layer+1)+'\n')
        for layer in range(self.load.n_layer):
            fid1.write(" OPEN/CLOSE '"+self.load.por_name+str(layer+1)+"."+self.load.dtype_por+"'  1  '(FREE)'  -1             porosity layer "+str(layer+1)+'\n')
        fid1.close()
        
        print 'done'

#==============================================================================
    def modflow_mult(self,num=100): 
#==============================================================================
        """
        MODFLOW Multiplier File.
        Input to define multiplier arrays is read from the file that is 
        specified with "MULT" as the file type.
        
        Parameters
        ----------
        num :               int
            identity(index) number in modflow nam-file
        """
        
        print 'MODFLOW MULT file:           ',
        fid1 = open(self.load.modflow_file_name+'.mult','w')
        fid1.write('# Multiplier file for MODFLOW-2000, UNC1NLI1 Report Example \n ')
        fid1.write(str(2*self.load.n_layer)+ ' \n') 
#        fid1.write(str(1)+ '\n')
        for i in range(self.load.n_layer):
            fid1.write('HK'+str(i+1)+ '\n')
            fid1.write('          '+str(num+i)+'         1(FREE)                      -1  4. HK Multiplier array for layer 1'+ '\n')
            fid1.write('VK'+str(i+1)+ '\n')
            fid1.write('          0           1(20G14.0)                   -1  4. VK Multiplier array for layer 1'+ '\n')
        fid1.close()
      
        
        
        print 'done'
    
#==============================================================================
    def modflow_nam(self,output_file,run_type='true',zone='off',well='off',river='off',initial_head='off',solver='GMG',hob='off',adj='off'):         
#==============================================================================
        """
        Write modflow input name file
        
        Parameters
        ----------
        output_file :       str
            output name of file with extension "nam"
        run_type :          str
            true    --> run modflow with true model
            calib   --> run modflow, GWM discretized by pilot points (real array
            files)
        zone :              str
            optional('on',off')
        well :              str
            optional('on',off')
        river :             str
            optional('on',off')
        initial_head :      str
            optional('on',off')
        solver :            str
            GMG -->
            PCG -->
        hob :               str
            optional('on',off')
        adj :                   str
            optional('on',off')
            Name of modflow adjoint input file (*.adj)
                     
             
        """
        print 'modflow_nam:               ',
        fid1 = open(output_file+'.nam','w')
        if run_type == 'true':
            fid1.write('%1.20s %10.0i %5.30s\n'%('GLOBAL',6,self.load.modflow_file_name+'.glo'))            
            fid1.write('%1.20s %10.0i %5.30s\n'%('LIST',7,self.load.modflow_file_name+'.lst'))
            fid1.write('%1.20s %10.0i %5.30s\n'%('BAS6',1,self.load.modflow_file_name+'.bas'))
            fid1.write('%1.20s %10.0i %5.30s\n'%('DIS',95,self.load.modflow_file_name+'.dis'))
            fid1.write('%1.20s %10.0i %5.30s\n'%('LPF',31,self.load.modflow_file_name+'.lpf'))
            fid1.write('%1.20s %10.0i %5.30s\n'%('RCH',18,self.load.modflow_file_name+'.rch'))
            if zone == 'on':
                fid1.write('%1.20s %10.0i %5.30s\n'%('ZONE',32,self.load.modflow_zone_file))
            if well == 'on':
                fid1.write('%1.20s %10.0i %5.30s\n'%('WEL',16,self.load.modflow_file_name+'.wel'))
            fid1.write('%1.20s %10.0i %5.30s\n'%('OC',22,self.load.modflow_file_name+'.oc'))
            fid1.write('%1.20s %10.0i %5.30s\n'%(solver,23,self.load.modflow_file_name+'.'+solver))
            if river == 'on':
                fid1.write('%1.20s %10.0i %5.30s\n'%('RIV',14,self.load.modflow_file_name+'.riv'))
            fid1.write('%1.20s %10.0i %5.30s\n'%('DATA(Binary)',50,self.load.modflow_file_name+'.cbb'))
            fid1.write('%1.20s %10.0i %5.30s\n'%('DATA(Binary)',51,self.load.modflow_file_name+'.hds'))
            fid1.write('%1.20s %10.0i %5.30s\n'%('DATA(Binary)',52,self.load.modflow_file_name+'.ddn'))
            if initial_head == 'on':
                fid1.write('%1.20s %10.0i %5.30s\n'%('DATA(Binary)',75,self.load.modflow_file_name+'.hds'))        
        
        if run_type == 'calib':
            fid1.write('%1.20s %10.0i %5.30s\n'%('LIST',7,self.load.modflow_file_name+'_calib'+'.lst'))
#            fid1.write('%1.20s %10.0i %5.30s\n'%('GLOBAL',61,self.load.modflow_file_name+'_calib'+'.glo'))            
            fid1.write('%1.20s %10.0i %5.30s\n'%('BAS6',1,self.load.modflow_file_name+'_calib'+'.bas'))
            fid1.write('%1.20s %10.0i %5.30s\n'%('DIS',95,self.load.modflow_file_name+'.dis'))
            fid1.write('%1.20s %10.0i %5.30s\n'%('LPF',31,self.load.modflow_file_name+'_calib'+'.lpf'))
            if self.load.rch_est < 0:
                fid1.write('%1.20s %10.0i %5.30s\n'%('RCH',18,self.load.modflow_file_name+'_calib'+'.rch'))
            if self.load.rch_est == 0:
                fid1.write('%1.20s %10.0i %5.30s\n'%('RCH',18,self.load.modflow_file_name+'.rch'))
            if zone == 'on':
                fid1.write('%1.20s %10.0i %5.30s\n'%('ZONE',32,self.load.modflow_zone_file))
            if well == 'on':
                fid1.write('%1.20s %10.0i %5.30s\n'%('WEL',16,self.load.modflow_file_name+'.wel'))
            fid1.write('%1.20s %10.0i %5.30s\n'%('OC',22,self.load.modflow_file_name+'.oc'))
            fid1.write('%1.20s %10.0i %5.30s\n'%(solver,23,self.load.modflow_file_name+'.'+solver))
            if hob == 'on':
                fid1.write('%1.20s %10.0i %5.30s\n'%('HOB',15,self.load.modflow_file_name+'.hob'))
            if river == 'on':
                fid1.write('%1.20s %10.0i %5.30s\n'%('RIV',14,self.load.modflow_file_name+'.riv'))
                
            if adj != 'off':
                fid1.write('%1.20s %10.0i %5.30s\n'%('MULT',94,self.load.modflow_file_name+'.mult'))
                fid1.write('%1.20s %10.0i %5.30s\n'%('ADJ',42,self.load.modflow_file_name+'.'+adj))
                fid1.write('\n')
                fid1.write('%1.20s %10.0i %5.30s\n'%('DATA',43,'lump_par.sen'))
                fid1.write('%1.20s %10.0i %5.30s\n'%('DATA',44,'dist_par.sen'))
                fid1.write('\n')
                fid1.write('%1.20s %10.0i %5.30s\n'%('DATA',38,self.load.modflow_file_name+'_calib_'+'head.hds'))
                fid1.write('%1.20s %10.0i %5.30s\n'%('DATA',39,self.load.modflow_file_name+'_calib_'+'heads_well.hds'))
                fid1.write('%1.20s %10.0i %5.30s\n'%('DATA',48,self.load.modflow_file_name+'_adj.pcg'))
                fid1.write('%1.20s %10.0i %5.30s\n'%('DATA',49,'obs_no.dat'))
                for i in range(self.load.n_layer):
                    fid1.write('%1.20s %10.0i %5.30s\n'%('DATA',100+i,self.load.modflow_file_name+'_kfield_L'+str(i+1)+'.dat'))
            else:
                fid1.write('%1.20s %10.0i %5.30s\n'%('DATA(Binary)',50,self.load.modflow_file_name+'_calib'+'.cbb'))
                fid1.write('%1.20s %10.0i %5.30s\n'%('DATA(Binary)',51,self.load.modflow_file_name+'_calib'+'.hds'))
                fid1.write('%1.20s %10.0i %5.30s\n'%('DATA(Binary)',52,self.load.modflow_file_name+'_calib'+'.ddn'))
                if initial_head == 'on':
                    fid1.write('%1.20s %10.0i %5.30s\n'%('DATA(Binary)',75,self.load.modflow_file_name+'_initial_head.hds'))
        fid1.close()
        
        print 'done'
#==============================================================================
    def modflow_oc(self,head=51,drawdown=52):        
#==============================================================================
        """
        write modflow OC file
        
        Parameters
        ----------
        head :              int
            deafault(51)
        drawdown            int
            default(52)

        """
        print 'modflow_oc:                  ', 
        
        fid1    = open(self.load.modflow_file_name+'.oc','w')
        fid1.write('HEAD SAVE UNIT '+str(head)+'\n')
        fid1.write('HEAD PRINT FORMAT 0\n')
        fid1.write('DRAWDOWN SAVE UNIT '+str(drawdown)+'\n')
        fid1.write('DRAWDOWN PRINT FORMAT 0\n')
        fid1.write('COMPACT BUDGET AUX\n')
        fid1.write('PERIOD 1  STEP 1\n')
        fid1.write('SAVE HEAD\n')
        fid1.write('SAVE DRAWDOWN\n')
        fid1.write('SAVE BUDGET\n')
        fid1.close()
        
        print 'done'
        
#==============================================================================
    def modflow_rch(self,run_type='true',sufix=''):
#==============================================================================
        """
        Write MODFLOW recharge module input
        
        Parameters
        ----------
        run_type :          str
            true    --> run modflow with true model
            calib   --> run modflow, GWM discretized by pilot points (real array
            files)
        sufix :             str, optional
            sufix to modflow nam-file, see HYTEB-ini-file
        """
        
        print 'modflow_rch:                 ',
        if run_type == 'true':
            fid1 = open(self.load.modflow_file_name+sufix+'.rch','w')
            fid1.write('# MODFLOW2000 Recharge Package\n')
            fid1.write('PARAMETER  0\n')
            fid1.write('    3        50\n')
            fid1.write('    1         1\n')
            fid1.write("OPEN/CLOSE '"+self.load.recharge+"."+self.load.dtype_recharge+"' 1.0  '(FREE) ' -1                 recharge\n")
            fid1.close()
        elif run_type == 'calib':
            fid1 = open(self.load.modflow_file_name+'_calib.rch','w')        
            fid1.write('# MODFLOW2000 Recharge Package\n')
            fid1.write('PARAMETER  0\n')
            fid1.write('    3        50\n')
            fid1.write('    1         1\n')
            tmp = self.load.rch_est_calc.split('.')
            fid1.write("OPEN/CLOSE '"+tmp[0]+"."+self.load.dtype_recharge+"' 1.0  '(FREE) ' -1                 recharge\n")
            fid1.close()
        
        elif run_type == 'zone':
            hz = self.load.hydro_zones('hydro_zones.ini',3,4)

            fid1 = open(self.load.modflow_file_name+sufix+'.rch','w')
            fid1.write('# MODFLOW2000 Recharge Package \n')
            fid1.write('PARAMETER '+str(self.load.n_zone_rch)+'\n')
            fid1.write('    3        50'+'\n')
            for i in range(2*self.load.n_zone_rch,3*self.load.n_zone_rch,1):

                fid1.write(str(hz[i,0])+'_'+str(hz[i,1])+' '+str(hz[i,0])+' '+str(hz[i,2])+' 1'+'\n')
                fid1.write('NONE ZLAY_1 '+str(hz[i-6,3])+'\n') 
            fid1.write(str(self.load.n_zone_rch)+'\n')
            
            for i in range(2*self.load.n_zone_rch,3*self.load.n_zone_rch,1):
                fid1.write(str(hz[i,0])+'_'+str(hz[i,1])+'\n')
            fid1.close()
            
        else:
            print 'WARNING MODFLOW_RCH: wrong input for, run_type!'             
            

        print 'done'
#==============================================================================
    def modflow_river(self,r_south= 0,r_north= 3500,r_east= 2500,hk_vertical=5e-3,stage=3.0,rbot=1.0,head_low=1.0):        
#==============================================================================
        """
        Calulates the position of the river.
        The river is setup as a "Gaining" River.
        The river is running straight from North to South.
        
        Parameters
        ----------
        r_south :           float
            Starting postion of river (southern boundary)
        r_north :           float
            Ending postion of river (Northern boundary)
        r_east :            float
            East-West location of river
        hk_vertical :       float
            vertical hk of riverbed
        stage :             float
            head in river below the initial GWM head 
        rbot :              float
            bottom of river, in meters below river_head
        head_low :          float
            stop criteria to substrating stage from initial gwm head
        """
        import subprocess
        
        inputfile_1 = 'river_obs_points'
#        inputfile_2 = 'heads_river_pos.smp'
        inputfile_3 = 'mod2smp_river.in'
        

        y_dim       = self.load.y
        dx          = self.load.dx
        dy          = self.load.dy
        
        # --------------------------------------------------
        #   --- write inputfile_2 with y coordinates  (x fixed) ---
        # ---------------------------------------------------
        y = np.arange(r_south+dy/2.0,r_north+dy,dy)
        y = np.flipud(y)
        x = r_east-dx/2.0
        
        fid1 = open(inputfile_1,'w')        
        for i in range(len(y)):
            fid1.write('Riv_'+str(i+1)+'    '+str(x)+'     '+str(y[i])+'     '+str(1)+'     '+str(0)+'\n')      
        fid1.close()
        
        #-----------------------------------------------------
        # --- run gwm without river to find location of river buttom
        #-----------------------------------------------------

        subprocess.call('mf2k.exe '+self.load.modflow_file_name+'_river.nam')
        
        self.write.mod2smp_in(inputfile_3,self.load.modflow_file_name,self.load.model_spc,inputfile_1)        
        
        myinput = open(inputfile_3,'r')
        tmp     = open('null.dat','w')
        subprocess.call('mod2smp.exe',stdin=myinput,stdout=tmp)
        myinput.close()
        tmp.close()
        
        # -------------------------------------------------------------
        #   --- import head measurements (smp-file) for 1 layer in initial gwm ---
        # -------------------------------------------------------------

        fid1 = open(inputfile_1+'.smp','r')
        text = np.array(fid1.read().split()).reshape(len(y),4)
        fid1.close()

        head = map(float,text[:,3])
        
        # --------------------------------------------------
        #      --- input parameters for RIVER file ---
        # --------------------------------------------------
       
        MXACTR = len(y)-1
        IRIVCB = 50
        ITMP   = len(y)-1
        NP     = 0
        LAYER  = 1
        cond   = hk_vertical
        
        fid2 = open(self.load.modflow_file_name+'.riv','w')
        fid2.write('#MODFLOW2000 RIVER PACKAGE INPUTFILE WRITTEN'+'\n')
        fid2.write(str(MXACTR)+'\t'+str(IRIVCB)+'\n')
        fid2.write(str(ITMP)+'\t'+str(NP)+'\n')
        
        node_north = (y_dim-r_north)/dy  
        node_south = (y_dim-r_south)/dy
        node_east  = r_east/dx 
                     
        # --------------------------------------------------------------------
        # --- calc where head in river is lower than sea level (head[i]=0) ---
        # --- when river-head is below head_low ( 1 meter)
        # --------------------------------------------------------------------

        
        count = 0 
        for i in range(len(y)):
            if head[i]-stage <= head_low:
                hp = head[i]-stage
                break
            count += 1
        
        h_stop = count
        step = hp/(len(y)-float(h_stop))
        
        count = 0 
        count2 =1
        for i in range(node_north,node_south+1,1):
            if i <= node_south-(len(y)-h_stop):
                tmp =[]                
                head_river = head[count]-stage
                tmp = head_river
                # ------------------------------------------------>         layer Row          Colum      Stage      Cond       Rbot
                fid2.write('%1.0f %10.0f %10.0f %10.3f %10.1e %10.3f \n'%(LAYER,i,node_east,head_river,cond,head_river-rbot))
            if i > node_south-(len(y)-h_stop):
                head_river = hp-step*count2
                # ------------------------------------------------>         layer Row          Colum      Stage      Cond       Rbot
                fid2.write('%1.0f %10.0f %10.0f %10.3f %10.1e %10.3f \n'%(LAYER,i,node_east,head_river,cond,head_river-rbot))   
                count2 += 1
            count +=1
        fid2.close()
        
        print 'done'  

#==============================================================================
    def modflow_rsp(self,name_file,start_loc,sufix=''): 
#==============================================================================
        """
        Write MODPATH rsp file.
        
        Parameters
        ----------
        name_file           str
            Modpath name file.            
        start_loc           str
            Name of data file containing startinglocations.      
        sufix               str
            sufix to modflow nam-file, see HYTEB-ini-file.    
        
        """
        
        fid2 = open(self.load.modflow_file_name+'_'+'modpath'+sufix+'.rsp','w')
        fid2.write(' @[MODPATH 5.0]'+'\n')
        fid2.write('* ENTER THE NAME FILE:'+'\n')
        fid2.write('@RESPONSE: '+'\n')
        fid2.write(name_file+'\n')
        fid2.write('* DO YOU WANT TO STOP COMPUTING PATHS AFTER A SPECIFIED LENGTH OF TIME ?'+'\n')
        fid2.write('@RESPONSE: '+'\n')         
        fid2.write('n'+'\n')
        fid2.write('* SELECT THE OUTPUT MODE:'+'\n')
        fid2.write('*     1 = ENDPOINTS'+'\n')
        fid2.write('*     2 = PATHLINE'+'\n')
        fid2.write('*     3 = TIME SERIES'+'\n')
        fid2.write('@RESPONSE:     '+'\n')     
        fid2.write('2'+'\n')
        fid2.write('* DO YOU WANT TO COMPUTE LOCATIONS AT SPECIFIC POINTS IN TIME ?'+'\n')
        fid2.write('@RESPONSE:        '+'\n')  
        fid2.write('n'+'\n')
        fid2.write('* HOW ARE STARTING LOCATIONS TO BE ENTERED?'+'\n')
        fid2.write('*     1 = FROM AN EXISTING DATA FILE'+'\n')
        fid2.write('*     2 = ARRAYS OF PARTICLES WILL BE GENERATED INTERNALLY'+'\n')
        fid2.write('@RESPONSE:   '+'\n')       
        fid2.write('1'+'\n')
        fid2.write('* ENTER NAME OF DATA FILE CONTAINING STARTING LOCATIONS:'+'\n')
        fid2.write('@RESPONSE:          '+'\n')
        fid2.write(start_loc+'\n')
        fid2.write('* IN WHICH DIRECTION SHOULD PARTICLES BE TRACKED?'+'\n')
        fid2.write('*     1 = FORWARD IN THE DIRECTION OF FLOW'+'\n')
        fid2.write('*     2 = BACKWARDS TOWARD RECHARGE LOCATIONS'+'\n')
        fid2.write('@RESPONSE:          '+'\n')
        fid2.write('1'+'\n')
        fid2.write('* HOW SHOULD PARTICLES BE TREATED WHEN THEY ENTER CELLS WITH INTERNAL SINKS ?'+'\n')
        fid2.write('*     1 = PASS THROUGH WEAK SINK CELLS'+'\n')
        fid2.write('*     2 = STOP AT WEAK SINK CELLS'+'\n')
        fid2.write('*     3 = STOP AT WEAK SINK CELLS THAT EXCEED A SPECIFIED STRENGTH'+'\n')
        fid2.write('@RESPONSE:          '+'\n')
        fid2.write('2'+'\n')
        fid2.write('* DO YOU WANT TO STOP PARTICLES WHENEVER THEY ENTER ONE SPECIFIC ZONE ?'+'\n')
        fid2.write('@RESPONSE:          '+'\n')
        fid2.write('n'+'\n')
        fid2.write('* DO YOU WANT TO COMPUTE VOLUMETRIC BUDGETS FOR ALL CELLS ?'+'\n')
        fid2.write('@RESPONSE:          '+'\n')
        fid2.write('n'+'\n')
        fid2.write('*  DO YOU WANT TO CHECK DATA CELL BY CELL ?'+'\n')
        fid2.write('@RESPONSE:          '+'\n')
        fid2.write('n'+'\n')
        fid2.write('* SUMMARIZE FINAL STATUS OF PARTICLES IN SUMMARY.PTH FILE ?'+'\n')
        fid2.write('@RESPONSE:       '+'\n')   
        fid2.write('n'+'\n')
        fid2.close()
#==============================================================================
    def modflow_rvob(self): 
#==============================================================================
        """ 
        
        Write modflow rvob file.
 
        """
        
        # --- load modflow river file ---
        fid1    = open('valley.riv','r')
        riv1     = fid1.readlines()
        fid1.close()
        fid1    = open('bud2smp_river.smp','r')
        riv2     = fid1.readline().split()
        fid1.close()
        
        
        NQRV        = 1
        NQCRV       = len(riv1)-4
        NQTRV       = 1
        IURVOBSV    = 1 
        TOMULTRV    = 1
        EVFRV       = 1 
        IOWTQRV     = 1
        NQOBRV      = 1
        NQCLRV      = 1 
        OBSNAM      = 'flow'
        IREFSP      = 1
        TOFFSET     = 0.0
        HOBS        = float(riv2[3])
        std         = 0.20
        STATISTIC   = 1/(HOBS*std)**2
        STAT_FLAG   = 1
        PLOT_SYMBOL = 1
        
        # --- load file with xy positions and layer ---       
        
        fid1 = open('valley.rvob','w')
        fid1.write('# River package package - RVOB'+'\n')
        fid1.write(str(NQRV)+'\t'+str(NQCRV)+'\t'+str(NQTRV)+'\t'+str(IURVOBSV)+'\n')
        fid1.write(str(TOMULTRV)+'\t'+str(EVFRV)+'\t'+str(IOWTQRV)+'\n')
        fid1.write(str(NQOBRV)+'\t'+str(NQCLRV)+'\n')
        #fid1.write('FLOW  1  0.0  1  '+str(riv2[3])+'\n')
        fid1.write(str(OBSNAM)+'\t'+str(IREFSP)+'\t'+str(TOFFSET)+'\t'+str(HOBS)+'\t'+str(STATISTIC)+'\t'+str(STAT_FLAG)+'\t'+str(PLOT_SYMBOL)+'\n')
        for i in range(len(riv1)-4):
            temp = riv1[4+i].split()
            fid1.write(str(temp[0])+'\t'+str(temp[1])+'\t'+str(temp[2])+'\t'+str(1.0)+'\n')
        fid1.close()
               
#==============================================================================
    def modflow_wel(self):         
#==============================================================================
        """
        Write MODFLOW2000 Well Package.
        
        """
        print 'make MODFLOW MPT file:       ',

        fid1 = open(self.load.modflow_file_name+'.wel','w')
        fid1.write('# MODFLOW2000 Well Package'+'\n')
        fid1.write('PARAMETER  0  0'+'\n')                    
        fid1.write('   1         0' +'\n')
        fid1.write('   1         0' +'\n')
        fid1.write(' '+str(self.load.n_layer)+'  204       100  -0.01500        -1'+'\n' )
        fid1.close()
        
        print 'done'

#==============================================================================
    def model_specification(self):
#==============================================================================
        """
        Write model_specification file for PEST groundwater utilities
        """
        
        print 'model_specification:         ',
        
        fid1 = open(self.load.model_spc,'w')
        fid1.write(str(self.load.ny)+' '+str(self.load.nx)+'\n')
        fid1.write('0.0'+' '+str(self.load.y)+' '+'0.0'+'\n')
        fid1.write(str(self.load.nx)+'*'+str(self.load.dx)+'\n')
        fid1.write(str(self.load.ny)+'*'+str(self.load.dy)+'\n')
        fid1.close()
        
        print 'done'
#==============================================================================
    def noisemodel2file(self,mean=0,std=1,number=40):
#==============================================================================
        """
        Writes a noise ini file with a noise model for fwr2tem.
        
        Parameters
        ----------
        mean :              float
            gaussian mean    
        std :               float
            1 STD
        number :            int
            length of noise-vector
        
        """
        
        print 'noise2model              ',
        
        # --- noise model as numpy vector ---
        G = np.random.normal(mean,std,number)
        
        # --- write G vector to file ---
        fid1 = open(self.load.hydro_noise_model_file,'w')
        for i in range(number):
            fid1.write('%9.11f \n'%(G[i]))
        fid1.close()
        
    
        print 'done'

#==============================================================================
    def obs2hob(self,TYPE='obs',MOBS=0,MAXM=0,IUHOBSV=42,HOBDRY= '',
                TOMULTH=1.0,IREFSP=1,TOFFSET=0.0,STAT_FLAG=1,
                PLOT_SYMBOL=1):
#==============================================================================
        """
        Rewrites PEST format bore file to MODFLOW HOB format
        
        Parameters
        ----------
        TYPE :              str
            Obs --> head observation with noise.
            Smp --> head observations without noise. 
        MOBS :              int
            See MODFLOW doc.
        MAXM :              int
            See MODFLOW doc.
        IUHOBSV :           int
            See MODFLOW doc.
        HOBDRY :            str
            See MODFLOW doc.
        TOMULTH :           float
            See MODFLOW doc.
        IREFSP :            int
            See MODFLOW doc.
        TOFFSET :           float
            See MODFLOW doc.
        STAT\_FLAG :         int
            See MODFLOW doc.
        PLOT\_SYMBOL :       int
            See MODFLOW doc.
        
        """
        
        # --- input data ---
        NH          = self.load.n_well
        STATISTIC   = 1.0/self.load.std_well**2
        
        
        # --- load file with xy positions and layer from "PEST" borehole format ---       
        fid1    = open(self.load.head_obs_points,'r')
        hop     = np.array(fid1.read().split()).reshape(self.load.n_well,5)
        fid1.close()
        
        # --- load bore obsdervation file with noise ---
        tmp = self.load.head_obs_points.split('.')
        fid1    = open(tmp[0]+'.'+TYPE,'r')
        head    = np.array(fid1.read().split()).reshape(self.load.n_well,4)
        fid1.close()
        
        # --- write MODFLOW HOB file ---
        fid1 = open(self.load.modflow_file_name+'.hob','w')
        fid1.write('#Head-Observation Package - HOB'+'\n')

        
        fid1.write(str(NH)+'\t'+str(MOBS)+'\t'+str(MAXM)+'\t'+str(IUHOBSV)+'\t'+str(HOBDRY)+'\n')
        fid1.write(str(TOMULTH)+'\n')
        for i in range(self.load.n_well):
            row_1   = float(hop[i,2])/self.load.dx
            roff    = np.mod(row_1,1)
            if roff >=0.5:
                row     = int(row_1-roff)+1
                row_1   = int(np.ceil(row_1))
                roff    = roff-1
            if roff <= 0.5:
                row     = int(row_1-roff)-1
                row_1   = int(np.floor(row_1))
                roff    = roff
                        
            col_1   = float(hop[i,1])/self.load.dy
            coff    = np.mod(col_1,1)
            if coff >=0.5:
                col     = int(col_1-coff)+1
                col_1   = int(np.ceil(col_1))
                coff    = coff-1
            if coff <= 0.5:
                col     = int(col_1-coff)-1
                col_1   = int(np.floor(col_1))
                coff    = coff
                        

            fid1.write(str(hop[i,0])+'\t\t'+str(hop[i,3])+'\t'+str(row_1)+'\t'+str(col_1)+'\t'+str(IREFSP)+'\t'+str(TOFFSET)+'\t')
#                          OBSNAM               LAYER              ROW             COLUMN          IREFSP           TOFFSET 
            fid1.write(str(roff)+'\t'+str(coff)+'\t'+str( head[i,3])+'\t'+str(STATISTIC)+'\t'+str(STAT_FLAG)+'\t'+str(PLOT_SYMBOL)+'\n')
#                          ROFF           COFF            HOBS                STATISTIC           STAT-FLAG           PLOT-SYMBOL

        fid1.close()

#==============================================================================
    def pred2obs(self,noise='on'):        
#==============================================================================
        """
        module to adding gaussian noise to hydraulic observations (head prediction)
        outfile has extension *.obs
        
        Parameters
        ----------
        noise :             str
            optinonal('on','off')
            
        """
        # --- noise model from file or automatic updated ---        
        if self.load.hydro_noise_model == -1:
            G = np.random.normal(0,1,self.load.n_well)
        elif self.load.hydro_noise_model == 0:
            G = self.load.noisemodel2vec(self.load.hydro_noise_model_file)
        else:
            print 'SMP2OBS -warning, '
            print 'wrong number for noise_model in hydro.ini'

            
        # ---open head observation file with true values ---
        predfile = self.load.head_pred_points.split('.')
        fid1 = open(predfile[0]+'.smp','r')
        fid2 = open(predfile[0]+'.obs','w')
        count = 0
        
        
        while count < self.load.n_pred_well:
            text = fid1.readline().split()
            if self.load.std_well == 0.0:
                print ' no noise added to well observation ' , 
                well_data   = float(float(text[3]))
                fid2.write(text[0]+'  '+text[1]+'  '+text[2]+'\t'+str(well_data)+'\n')
                count += 1
            else:
                if noise == 'on':
                    noise       = G[count]                    # adding noise to data
                else:
                    noise = 0.0
                    
                well_data   = float(float(text[3])+noise)
                fid2.write(text[0]+'  '+text[1]+'  '+text[2]+'\t'+str(well_data)+'\n')
                count += 1 
            
        print 'done'       
        
#==============================================================================
    def run_hydro_obs(self,modpath='off'):        
#==============================================================================
        """
        writing input file for mod2smp and bud2smp
        Used to generate hydarulic observation without noise.
        True head values at the head prediction points is also generatede and saved. 
        outputfile has extension *.smp
        
        Parameters
        ----------
        modpath :           str
            optinonal('on','off')
            
        """
        import subprocess
        print 'run_hydro_obs    '
        
        head_obs_file = [self.load.head_obs_points,self.load.head_pred_points]
        
        # --- write in files for mod2smp ---
        for i in head_obs_file:
            output_file         = i.split('.')
            output_file_name    = 'mod2smp_'+output_file[0]+'.in'
            self.write.mod2smp_in(output_file_name,self.load.modflow_file_name,self.load.model_spc,i)
            print '     write:               '+output_file_name
            
            print '     run mod2smp:        ',
            myinput = open(output_file_name,'r')
            tmp     = open('null.dat','w')
            subprocess.call('mod2smp.exe',stdin=myinput,stdout=tmp)
            myinput.close()
            tmp.close()       
            print 'done'            
        
        output_file = 'bud2smp_river'
        print '     write:               '+output_file
        self.write.bud2smp_in(output_file,self.load.modflow_file_name,self.load.model_spc,self.load.n_layer,self.load.ibound)
        print '     run bud2smp:        ',
        myinput = open(output_file+'.in','r')
        tmp     = open('null.dat','w')
        subprocess.call('bud2smp.exe',stdin=myinput,stdout=tmp)
        myinput.close()
        tmp.close() 
        
        if modpath == 'on':
            try:                
                fid1 = open('Mpathr5_0.bat','w')
                fid1.write('Mpathr5_0.exe < mpathr5_0_pumping.ini \n')
                fid1.close()
#                myinput = open('valley_modpath_pumping.rsp','r')
#                myinput = 'valley_modpath_pumping.rsp'
#                tmp     = open('null.dat','w')
#                subprocess.call('Mpathr5_0.exe',stdin=myinput)#,stdout=tmp)     
                myinput.close()
#                tmp.close() 
                subprocess.call('Mpathr5_0.bat')#,stdout=tmp)
                
            except:
                print 'Warning16 in run_prediction'
                print 'MPATHR_5_0 fail in run with pumping'
                                      
            try:
                self.write.sortbyept('endpoint',self.load.nx,self.load.ny,1,100,204,self.load.n_layer)
                src = 'rchsum_true.dat'
                
                subprocess.call('sortbyept.exe')
                fid1 = open(src,'w')
                subprocess.call('rchsum.exe',stdout=fid1)
                fid1.close()
                
            except:
                print 'WARNING17, problems in sortbyept.exe or rchsum.exe'

        elif modpath == 'off':
            pass
        else:
            'Warning, wrong input for modpath!!!'
        
        print 'done' 
        
#==============================================================================
    def run_mf2k(self,input_file,save_initial_head=None):        
#==============================================================================
        """
        Run MODFLOW-2000
        
        Parameters
        ----------
        input_file :        str
            MODFLOW nam-file
        save_initial_head : logical
            optional(True/None)
            
        Returns
        -------
        file :              str
            if save_inital_head --> True
            MODFLOW-2000 filename + '_initial_head.hds'
            
        """
        import shutil
        import subprocess
        
        #-----------------------------------------------------
        # --- run gwm 
        #-----------------------------------------------------

        subprocess.call('mf2k.exe '+input_file+'.nam')
        if save_initial_head == 'on':            
            shutil.copyfile(self.load.modflow_file_name+'.hds',self.load.modflow_file_name+'_initial_head.hds')
        
#==============================================================================
    def run_prediction(self,model_num,run_type,modpath='on'):
#==============================================================================
        """
        Script to groundwater model ini forward mode for 1000 realizations.
        The GWM is runned WITHOUT pumping well and some predictions are made
        Then the GWM is runned WITH pumping well and some other predictions are made
        
        Parameters
        ----------
        model_num :         int
            geological reference structure number (realization number)
        run_type :          str
            true    --> run modflow with true model
            calib   --> run modflow, GWM discretized by pilot points (real array
            files)        
        modpath :           str
            optional('on','off')
        
        
        """
        print run_type 
        import shutil
        import subprocess
        
        # --- copy PEST output files ---
        dtype = ['jco','par','pst','rec','rei','res','sen','seo','svd']
        for i in dtype:
            try:                
                src = self.load.PEST_name+'.'+i
                dst = self.load.path_results+'/'+src
                shutil.copyfile(src,dst+'.'+str(model_num))  
                del src, dst
            except:
                print ' PEST file:   *.', i,' not copied'
                pass 
         
        
        # --- copy calibrated hydraulic conductivity files ---                
        for i in range(self.load.n_layer):
            try:
                src = self.load.pp_file_name+str(i+1)+'.'+self.load.dtype_model
                dst = self.load.path_results+'/'+self.load.pp_file_name+str(i+1)+'.'+self.load.dtype_model
                shutil.copyfile(src,dst+'.'+str(model_num))          
            except:
                print ' Hydraulic conductivty not copied!'  

        # --- copy calibrated recharge file ---  
        if self.load.rch_est == -1:                      
            try:
                src = 'recharge_est.ref'
                dst = self.load.path_results+'/'+src
                shutil.copyfile(src,dst+'.'+str(model_num))          
            except:
                print ' Recharge_est not copied!'
        
                  
        #--------------------------------------------------------------
        #--- run GWM WITHOUT PUMPING ---        
        #--------------------------------------------------------------
                
        try:
            if run_type == 'true':
                subprocess.call('mf2k.exe '+self.load.modflow_file_name+'.nam') 
            elif run_type == 'calib':
                subprocess.call('mf2k.exe '+self.load.modflow_file_name+'_calib.nam')
            else:
                print 'WARNING1, wrong input value for run_type '
        except:
            print 'Warning2 in run_prediction'
            print 'MODFLOW fail in run WITHOUT pumping'

        try:                
            fid1 = open('Mpathr5_0.bat','w')
            fid1.write('Mpathr5_0.exe < mpathr5_0.ini' +'\n')
            fid1.close()
            subprocess.call('Mpathr5_0.bat')     
            
        except:
            print 'Warning3 in run_prediction'
            print 'MPATHR_5_0 fail in run with pumping'
                        
        # --- copy MODFLOW output files ---
        dtype = ['hds','lst']
        for i in dtype:
            try:
                if run_type == 'true':                
                    src = self.load.modflow_file_name+'.'+i
                elif run_type == 'calib':                
                    src = self.load.modflow_file_name+'_calib.'+i
                else:
                    print 'WARNING4, wrong input value for run_type!'
                    
                dst = self.load.path_results+'/'+src
                shutil.copyfile(src,dst+'.'+str(model_num))  
                del src, dst
            except:
                print ' Warning5 in run_prediction'
                print ' MODFLOW file:   *.', i,' not copied'
                pass           
        
        # --- calc head prediction ---         
        try:            
            if run_type == 'true':
                self.write.mod2obs_in('mod2obs_pred',self.load.modflow_file_name,self.load.model_spc,self.load.head_pred_points,self.load.n_layer)                 
            elif run_type == 'calib':    
                self.write.mod2obs_in('mod2obs_pred',self.load.modflow_file_name+'_calib',self.load.model_spc,self.load.head_pred_points,self.load.n_layer)                 
            myinput = open('mod2obs_pred.in','r')
            tmp     = open('null.dat','w')
            subprocess.call('mod2obs.exe',stdin=myinput,stdout=tmp)
            myinput.close()
            tmp.close()            
        except:
            print 'Warning6 in run_prediction:   mod2obs!'
        
        # --- copy MODFLOW head predictions ---
        try:
            tmp = self.load.head_pred_points.split('.')
            src = tmp[0]+'_model.smp'
            dst = self.load.path_results+'/'+src
            shutil.copyfile(src,dst+'.'+str(model_num))             
            del dst
        except:
                print ' Warning6 in run_prediction'
                print ' MODFLOW file:   *.', src,' not found'
                del src
                pass
            
        # --- read prediction from modpath file and print to file --- 
        try:
            os.system('python finaltime.py')
            src = 'time.dat'
            dst = self.load.path_results+'/'+src
            shutil.copyfile(src,dst+'.'+str(model_num))
        except:
            print 'Warning7 in run_prediction:       finaltime.py!'
            
        # --------------------------------------------------------------------
        # --- delete some intemediate files---
        # --------------------------------------------------------------------    
        try:
            tmp = self.load.head_pred_points.split('.')
            src = tmp[0]+'_model.smp'
            os.remove(src)
            os.remove('pathline')
            os.remove('endpoint')
        except:
            print 'WARNING8; cannot remove *.smp, pathline or endpoint'
        
        dtype = ['hds','lst','ddn']    
        for i in dtype:     
            try:
                if run_type == 'true':                
                        src = self.load.modflow_file_name+'.'+i
                        os.remove(src)
                elif run_type == 'calib':                
                        src = self.load.modflow_file_name+'_calib.'+i
                        os.remove(src)
                else:
                    print 'WARNING9, wrong input value for run_type '
            except:
                print 'fail to remove  *.',i,' or not found!'
                
        # --------------------------------------------------------------------
        # --- run GWM WITH PUMPING ---        
        # --------------------------------------------------------------------              
        try:
            if run_type == 'true':
                subprocess.call('mf2k.exe '+self.load.modflow_file_name+'_pumping.nam')               
            elif run_type == 'calib':
                subprocess.call('mf2k.exe '+self.load.modflow_file_name+'_calib_pumping.nam')
            else:
                print 'WARNING10, wrong input value for run_type '
        except:
            print 'Warning11 in run_prediction'
            print 'MODFLOW fail in run WITH pumping'
               
               
        # --- copy MODFLOW output files ---            
        dtype = ['hds','lst']
        for i in dtype:
            try:
                if run_type == 'true':                
                    src = self.load.modflow_file_name+'.'+str(i)
                    dst = self.load.path_results+'/'+self.load.modflow_file_name+'_pump.'+str(i)
                elif run_type == 'calib':                
                    src = self.load.modflow_file_name+'_calib.'+str(i)
                    dst = self.load.path_results+'/'+self.load.modflow_file_name+'_calib_pump.'+str(i)
                else:
                    print 'WARNING12, wrong input value for run_type! '
                    
                shutil.copyfile(src,dst+'.'+str(model_num))  
                del src, dst
            except:
                print ' Warning13 in run_prediction'
                print ' MODFLOW file:   *.', i,' not copied'
                pass           
        
        # --- calc head prediction ---         
        try:
            if run_type == 'true':
                self.write.mod2obs_in('mod2obs_pred_pump',self.load.modflow_file_name,self.load.model_spc,self.load.head_pred_points,self.load.n_layer)                 
            if run_type == 'calib':    
                self.write.mod2obs_in('mod2obs_pred_pump',self.load.modflow_file_name+'_calib',self.load.model_spc,self.load.head_pred_points,self.load.n_layer)                 
            myinput = open('mod2obs_pred_pump.in','r')
            tmp     = open('null.dat','w')
            subprocess.call('mod2obs.exe',stdin=myinput,stdout=tmp)
            myinput.close()
            tmp.close()            
        except:
            print 'Warning14 in run_prediction:   mod2obs!'    
        
        # --- copy MODFLOW head predictions ---
        try:
            tmp = self.load.head_pred_points.split('.')
            if run_type == 'true':
                src = tmp[0]+'_model.smp'
            if run_type == 'calib':
                src = tmp[0]+'_model.smp'
                
            dst = self.load.path_results+'/'+tmp[0]+'_model_pump.smp'
            shutil.copyfile(src,dst+'.'+str(model_num))             
            del dst
        except:
                print ' Warning15 in run_prediction'
                print ' MODFLOW file:   *.', src,' not found'
                del src
                pass
            
        if modpath == 'on':
            try:                
                fid1 = open('Mpathr5_0.bat','w')
                fid1.write('Mpathr5_0.exe < mpathr5_0_pumping.ini \n')
                fid1.close()
                subprocess.call('Mpathr5_0.bat')     
                
            except:
                print 'Warning16 in run_prediction'
                print 'MPATHR_5_0 fail in run with pumping'
                          
            
            try:
                self.write.sortbyept('endpoint',self.load.nx,self.load.ny,1,100,204,self.load.n_layer)
                src = 'rchsum.dat'
                dst = self.load.path_results+'/'+src
                
                subprocess.call('sortbyept.exe')
                fid1 = open(src,'w')
                subprocess.call('rchsum.exe',stdout=fid1)
                fid1.close()
                
                shutil.copyfile(src,dst+'.'+str(model_num))
            except:
                print 'WARNING17, problems in sortbyept.exe or rchsum.exe'
            
            


#==============================================================================
    def smp2obs(self):        
#==============================================================================
        """
        Module to adding gaussian noise to hydraulic observations (head and river)
        outfile has extension *.obs
        """
        # --- noise model from file or automatic updated ---
        if self.load.hydro_noise_model == -1:
            G = np.random.normal(0,1,self.load.n_well)
        elif self.load.hydro_noise_model == 0:
            G = self.load.noisemodel2vec(self.load.hydro_noise_model_file)
        else:
            print 'SMP2OBS -warning, '
            print 'wrong number for noise_model in hydro.ini'
            
            
        # ---open head observation file with true values ---
        obsfile = self.load.head_obs_points.split('.')
        fid1 = open(obsfile[0]+'.smp','r')
        fid2 = open(obsfile[0]+'.obs','w')
        count = 0
        
        
        while count < self.load.n_well:
            text = fid1.readline().split()
            if self.load.std_well == 0.0:
                print ' no noise added to well observation ' , 
                well_data   = float(float(text[3]))
                fid2.write(text[0]+'  '+text[1]+'  '+text[2]+'\t'+str(well_data)+'\n')
                count += 1
            else:
                noise       = G[count]                    # adding noise to data
                well_data   = float(float(text[3])+noise)
                fid2.write(text[0]+'  '+text[1]+'  '+text[2]+'\t'+str(well_data)+'\n')
                count += 1 
            
        # ---open river flux observation file with true values ---

        fid1 = open('bud2smp_river'+'.smp','r')
        fid2 = open('bud2smp_river'+'.obs','w')
        count = 0
        while count < self.load.n_river:
            text = fid1.readline().split()
            if self.load.noise_river == 0.0:
                river_data   = float(float(text[3]))
                fid2.write(text[0]+'  '+text[1]+'  '+text[2]+'\t'+str(river_data)+'\n')
                count += 1    
            else:
                std_river   = self.load.noise_river*abs(float(text[3]))/100.0
                noise       = np.random.normal(0,std_river,1)                    # adding noise to data
                river_data   = float(float(text[3])+noise)
                fid2.write(text[0]+'  '+text[1]+'  '+text[2]+'\t'+str(river_data)+'\n')
                count += 1                   
        fid1.close()
        fid2.close()
        print 'done'
        
#==============================================================================
    def strat2zns(self,fname):
#==============================================================================
        """
        Load (TPROGS) start-files and write MODFLOW zone file
        
        Parameters
        ----------
        fname :             str
            strat file name


        Returns
        -------
        *.zns  
        
        mat :      array_like, int
             numpy array with lithologies in 3D (node_y,node_x,node_z)

        
        """
        
        print 'strat2zns:           ',
        
        # --- load strat file --- 
        mat = self.load.strat(fname,self.load)

        # --- write modflow zone file ---
        self.write.modflow_zone(self.load.modflow_file_name,mat,self.load)
        
        print 'done'
        
        return mat
