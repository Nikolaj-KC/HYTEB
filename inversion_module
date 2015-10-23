# -*- coding: utf-8 -*-
"""
Created on Tue Aug 12 12:06:16 2014

@author: Nikolaj
"""

import os
import numpy as np
import math as mat
import subprocess
import shutil

import load_write



class inversion_module:
    """
    class containing modules to setup a PEST inversion in the HYTEB enviroment
   

    """   
    
    
#==============================================================================
    def __init__(self,inversion_ini='inversion_setup.ini', model_dis ='model_dis.ini', hydro_ini = 'hydro_ini.ini',geophys_ini='geophys_TEM.ini'):
#==============================================================================
        self.write = load_write.write_file()              
        self.load = load_write.load_file()        
        
        self.load.ini_model_dis(model_dis)
        self.load.ini_hydro(hydro_ini)
        self.load.ini_inversion(inversion_ini)
        self.load.ini_geophys(geophys_ini)
        
#        # --- load from inversion_setup.ini ---
#        self.load.pp_file_name       = load.pp_file_name
#        self.load.pp_name            = load.pp_name
#
#
#        # --- load from model_dis.ini ---        
#        self.load.dx                 = load.dx
#        self.load.dy                 = load.dy
#        self.load.dz                 = load.dz
#        self.load.n_layer            = load.n_layer
#        self.load.n_cap              = load.n_cap
#        self.load.x                  = load.x
#        self.load.y                  = load.y
#        self.load.nx                 = load.nx
#        self.load.ny                 = load.ny
#        self.load.n_pp_x             = load.n_pp_x
#        self.load.n_pp_y             = load.n_pp_y
#        
#        # --- load from hydro.ini.ini ---
#        self.load.model_name         = load.model_name
#        self.load.dtype_model        = load.dtype_model
#        self.load.ibound             = load.ibound
#        self.load.dtype_ibound       = load.dtype_ibound
#        self.load.por_name           = load.por_name
#        self.load.dtype_por          = load.dtype_por
#        self.load.recharge           = load.recharge
#        self.load.dtype_recharge     = load.dtype_recharge
#        self.load.model_name_calib   = load.model_name_calib
#        self.load.modflow_file_name  = load.modflow_file_name
#==============================================================================
    def batch_all_models(self,mf_nam,pp_x,pp_y,grid=1,pp_file_name=[]):
#==============================================================================
        """
        Write batch file for running both the geophysical and hydrolodical models
        
        Parameters
        ----------
        mf_nam :            str
            string with name of modflow nam file for running the calibration      
        pp_x :              List(int)  
            python list - node positions of pilot points in column direction
        pp_y :              List(int)
            python list - node positions of pilot points in row direction
        grid :              int
            grid number used when running model, see ini file for discretization        
        pp_file_name :      Optional(List(str))
            Python list with sufix to fac2real-ini files 

        """        
        
        print 'batch_allmodels:           '
        
        # --- open file for writing ---
        fid1 = open('allmodels.bat','w')  
        
        fid1.write('@echo off'+'\n')
        fid1.write('echo ALLMODELS'+'\n')
        fid1.write('Rem ###################################'+'\n')
        fid1.write('Rem DELETE MODFLOW OUTPUT FILES'+'\n')
        fid1.write('Rem ###################################'+'\n')
        fid1.write('\n')
        fid1.write('del '+self.load.modflow_file_name+'_calib.hds'+'\n')
        fid1.write('del '+self.load.modflow_file_name+'_calib.cbb'+'\n')
        fid1.write('del '+self.load.modflow_file_name+'_calib.dnn'+'\n')
        if self.load.rch_est != 0:
            tmp = self.load.rch_est_calc.split('.')
            fid1.write('del '+tmp[0]+'.ref'+'\n')
        fid1.write('del hk_diff_L*.dat'+'\n')
        fid1.write('del '+self.load.model_name_calib+'*.ref'+'\n')    
        fid1.write('\n')
        fid1.write('\n')
        fid1.write('Rem ###################################'+'\n')
        fid1.write('Rem Now the actual model is run'+'\n')
        fid1.write('Rem ###################################'+'\n')
        fid1.write('\n')
        fid1.write('for %%k in (')
        
        # --- run all geophysical measurements ---
        for i in range(len(pp_y)):
            fid1.write(str(pp_y[i])+' ')
        fid1.write( ') do ('+'\n')
        fid1.write('    for %%j in (')
        for j in range(len(pp_x)):
            fid1.write(str(pp_x[j])+' ')
        fid1.write(') do ('+'\n')        
        fid1.write('        '+self.load.geophys_software+'.exe  '+self.load.mod_file+'%%k_%%j.mod > nul'+'\n')
        fid1.write('   )'+'\n')
        fid1.write(')'+'\n')
        fid1.write('\n')
        
        # --- run PEST utility to interpolate parameter values from pilot point 
        # --- file to grid cells 
        fid1.write('set /A %F = '+str(self.load.n_layer)+'\n')
        fid1.write('\n')        
        fid1.write('for /L %%k in (1,1,%F%) do ('+'\n')
        
        # --- use grid defind in HYTEB ini ---
        if grid == 1:             
                fid1.write('   	fac2real.exe < fac2real_L%%k.in > nul'+'\n')
        if grid == 2:
            for pp_name in pp_file_name:
                fid1.write('   	fac2real.exe < fac2real_L%%k'+pp_name+'.in > nul'+'\n')
        fid1.write(')'+'\n')
        fid1.write('\n')
        
        #--- petrophysical parameter estimation with pilot points ---
        if self.load.petro_est ==-3:
            fid1.write(self.load.petro_est_calc+'\n') 
            
        # --- calc recharge qua linear log log relatioship ---    
        if self.load.rch_est != 0:
            fid1.write(self.load.rch_est_calc+'\n')             
        fid1.write('\n') 
        
        # --- run regularization ---
        fid1.write(self.load.reg_calc+'\n')
        
        # --- run modflow-2000 ---
        fid1.write('mf2k.exe '+mf_nam+'.nam > nul'+'\n')
        
        # --- run PEST utilities to read binary output data into acis files ---
        fid1.write('mod2obs.exe < mod2obs.in > nul'+'\n')
        fid1.write('bud2smp.exe < bud2smp_model.in > nul'+'\n')
        fid1.write('\n')      

        fid1.close()

        
        
#==============================================================================
    def batch_beopest(self,run_inversion ='off'):
#==============================================================================
        """
        Writes batch file for running beopest
        
        Parameter
        ---------
        run_inversion :     str
            off --> do not run PEST 
            on  --> run PEST 
        
        """

        print 'batch_beopest:            ',
        
        
        # --- write batch file for running beopest ---
        self.write.batch_beopest(self.load.IP_addr,self.load)
        
        # --- run beopest ---
        if run_inversion == 'on':
            if self.load.PEST_niter != 0:
                tmp     = open('run_beopest_output.dat','w')
                subprocess.call('run_beopest.bat',stdout=tmp)
                tmp.close()
            # --- if iteration is zero call PEST instead of BeoPEST ---
            if self.load.PEST_niter == 0:
                program = 'pest.exe'
                input1  = self.load.PEST_name+'.pst'
                dos_tex = [program,input1]
                outputtext= subprocess.call(dos_tex)
        
        
        
#==============================================================================
    def batch_gp_model(self):        
#==============================================================================
        """
        Makes batch stript to run one geophysical model. Using batcht input
        parameter to call geophysical forward response. 
        
        
        """
        print 'batch_gp_model:            ',
        
        fid1 = open('gp_model.bat','w')
        fid1.write('@echo off'+'\n')
        fid1.write('echo TEM_MODEL'+'\n')
        fid1.write('Rem ###################################'+'\n')
        fid1.write('Rem Now the actual model is run'+'\n')
        fid1.write('Rem ###################################'+'\n')
        fid1.write('\n')
        fid1.write('\n')        
        
        # --- running python script for calculating regularization ---- 
        fid1.write(self.load.reg_calc+'\n')
        
        # --- call geophysical software for running one forward model --- 
        fid1.write(self.load.geophys_software+'.exe  '+self.load.mod_file+'%1_%2.mod > run_info_Aarhus_Inv.dat'+'\n')
        fid1.close()
        
        print 'done'   
        
#==============================================================================
    def batch_gw_model(self,mf_nam,pp_x,pp_y,grid=1,pp_file_name=[],run_fac2real=True):
#==============================================================================
        """
        Makes batch stript to run groundwater model.
        
        Parameters
        ----------
        mf_nam :            str
            name of modflow nam file   
        grid :              int
            grid number used when running model, see ini file for discretization
        pp_file_name :      Optional(List((str))
            Python list with sufix to fac2real-ini files 
        run_fac2real :      logical
            True --> run fac2real to generate hk-field files. One for each layer
                
        """
        
        print 'batch_gw_model:             ',
        
        fid1 = open('gw_model.bat','w')
        fid1.write('@echo off'+'\n')
        fid1.write('echo GW_MODELS'+'\n')
        fid1.write('Rem ###################################'+'\n')
        fid1.write('Rem DELETE MODFLOW OUTPUT FILES'+'\n')
        fid1.write('Rem ###################################'+'\n')
        fid1.write('\n')
        fid1.write('del '+self.load.modflow_file_name+'_calib.hds'+'\n')
        fid1.write('del '+self.load.modflow_file_name+'_calib.cbb'+'\n')
        fid1.write('del '+self.load.modflow_file_name+'_calib.ddn'+'\n')
        
        # --- if gwm is running in prediction mode delete modpath output files ----        
        if self.load.pred2calibdata < 0:
            if self.load.pred2calibdata == -1 or self.load.pred2calibdata == -3: 
                fid1.write('del endpoint'+'\n')
                fid1.write('del pathline'+'\n')
                fid1.write('\n')
        
        # --- calc recharge qua linear log log relatioship ---  
        if self.load.parametrization == -1:
            if self.load.rch_est != 0:
                if run_fac2real == True:
                    tmp = self.load.rch_est_calc.split('.')
                    fid1.write('del '+tmp[0]+'.ref'+'\n')
            fid1.write('del hk_diff_L*.dat'+'\n')
            if run_fac2real == True:
                fid1.write('del '+self.load.model_name_calib+'*.ref'+'\n')    
            
        fid1.write('\n')
        fid1.write('\n')
        fid1.write('Rem ###################################'+'\n')
        fid1.write('Rem Now the actual model is run'+'\n')
        fid1.write('Rem ###################################'+'\n')
        fid1.write('\n')
        fid1.write('\n')
        
        # --- run PEST utility to interpolate parameter values from pilot point 
        # --- file to grid cells
        if run_fac2real == True:
            if self.load.parametrization == -1 or self.load.parametrization ==-3:
                fid1.write('set /A %F = '+str(self.load.n_layer)+'\n')
                fid1.write('\n')        
                fid1.write('for /L %%k in (1,1,%F%) do ('+'\n')
                
                # --- use grid defind in HYTEB ini ---
                if grid == 1:             
                    fid1.write('   	fac2real.exe < fac2real_L%%k.in > nul'+'\n')
                if grid == 2:
                    for pp_name in pp_file_name:
                        fid1.write('   	fac2real.exe < fac2real_L%%k'+pp_name+'.in > nul'+'\n')
                        
                fid1.write(')'+'\n')
                fid1.write('\n')
            


            #--- petrophysical parameter estimation with pilot points ---
            if self.load.petro_est ==-3:
                fid1.write(self.load.petro_est_calc+'\n') 
                
            # --- calc recharge qua linear log log relatioship ---    
            if self.load.rch_est == -1:
                fid1.write(self.load.rch_est_calc+'\n')
                fid1.write('\n')
                    
            # --- run regularization ---   
            if self.load.parametrization == -1 or self.load.parametrization ==-3: 
                if self.load.inv_type != 'hydrologic':
                    if self.load.petro_est == 0 or self.load.petro_est == -1:
                        fid1.write(self.load.reg_calc+'\n')
                    elif self.load.petro_est == -2:
                        tmp1 = self.load.reg_calc.split('.')
                        
                        fid1.write('set /A %F = '+str(self.load.n_layer)+'\n')      
                        fid1.write('for /L %%k in (1,1,%F%) do ('+'\n')
                        fid1.write('   	'+tmp1[0]+'_L%%k.'+tmp1[1]+'\n')
                        fid1.write(')'+'\n')
                        fid1.write('\n')
        
        # --- run modflow-2000 ---                        
        fid1.write('\n')    
        fid1.write('mf2k.exe '+mf_nam+'.nam >> mf2k.dat'+'\n')
        
        # --- run PEST utilities to read binary output data into acis files -- 
        fid1.write('mod2obs.exe < mod2obs.in > nul'+'\n')
        fid1.write('bud2smp.exe < bud2smp_model.in > nul'+'\n')
        fid1.write('\n')
        
        
        # --- add predictions to calibration dataset ---
        if self.load.pred2calibdata < 0:
            if self.load.pred2calibdata == -1 or self.load.pred2calibdata == -3: 
                fid1.write('Mpathr5_0.exe < mpathr5_0_pumping.ini \n')
                fid1.write('sortbyept.exe'+'\n')
                fid1.write('rchsum.exe > rchsum.dat'+'\n')
                fid1.write('\n')
                
            if self.load.pred2calibdata == -1 or self.load.pred2calibdata == -2:
                fid1.write('mod2obs.exe < mod2obs_pred_pump.in > nul'+'\n')
            
                self.write.mod2obs_in('mod2obs_pred_pump',self.load.modflow_file_name+'_calib',self.load.model_spc,self.load.head_pred_points,self.load.n_layer)                 
        
        fid1.close()        
        
        # --- write ini files for mod2obs.exe and bud2smp.exe ---
        self.write.bud2smp_in('bud2smp_model',self.load.modflow_file_name+'_calib',self.load.model_spc,self.load.n_layer,self.load.ibound)
        self.write.mod2obs_in('mod2obs',self.load.modflow_file_name+'_calib',self.load.model_spc,self.load.head_obs_points,self.load.n_layer)        
        
        print 'done'
        
#==============================================================================
    def copy2slaves(self,pp_x,pp_y,run_allmodels='on',tmp_numb=''):
#==============================================================================
        """
        Creates n slave dir.

        Make n new slave subdirectories and copy pest and model files to folder

        Parameters
        ----------
        pp_x :              List(int)  
            python list - node positions of pilot points in column direction
        pp_y :              List(int)
            python list - node positions of pilot points in row direction
        run_allmodel :      str
            run batch file "allmodels" before copying files into slave dir
        tmp_nump :          int
            
        
        """
        import errno
        
        print 'copy2slave:'
        
        cwd = os.getcwd()
        path ='C:/temp/HYTEB'+str(tmp_numb)
        
        for i in range(self.load.n_slaves):
            try:
                shutil.rmtree(path+str(tmp_numb))
            except: 
#                ishutil.Error as e:
                print 'copy2slaves: folder does not excist!!! '
        
        # --- delete temporay geophys files ---  
        for i in pp_y:
            for j in pp_x:
                try:                
                    os.remove(self.load.tem_name_raw+str(i)+'_'+str(j)+'.tem')
                except:
                    if IOError:
                        pass                    
                try:                        
                    os.remove(self.load.mod_file_temp+str(i)+'_'+str(j)+'.mod')
                except:
                    if IOError:
                        pass
        
        # --- delete modpath output files ---                
        try:
            os.remove('pathline')
        except:
                pass
         
        # --- run all_models ---
        if run_allmodels == 'on':
            if self.load.inv_type == 'joint':
                subprocess.call('allmodels.bat')
            if self.load.inv_type == 'sequential':
                subprocess.call('gw_model.bat')
                       
        print 'done'                                           
        # --- copy work folde to slave# and move folders back into work folder
        for i in range(self.load.n_slaves):
            # --- copy to dir ---            
            print str(i+1)+'/'+str(self.load.n_slaves)
            try: 
                shutil.copytree(cwd,path+'/slave'+str(i+1))
            # --- Directories are the same
            except shutil.Error as e:
                print 'Error in copy2slave:'
                print('Directory not copied. Error: %s' % e)
            # --- Directories already exist, delete and try copy again
            except OSError as e:
                try:
                    shutil.rmtree(path+'/slave'+str(i+1))
                except shutil.Error as e:
                    print 'Error in copy2slave:'
                    print('Directory not copied. Error: %s' % e)    
                try:
                    shutil.copytree(cwd,path+'/slave'+str(i+1))
                except shutil.Error as e:
                    print 'Error in copy2slave:'
                    print('Directory not copied. Error: %s' % e)
                # --- Directories already exist, delete and try copy again
                except OSError as e:
                    print 'Error in copy2slaves:'
                    print('Directory not copied. Error: %s' % e)
                    
        
        for i in range(self.load.n_slaves):
            try:
                shutil.move(path+'/slave'+str(i+1),cwd)
            # --- Directories are the same
            except shutil.Error as e:
                print 'Error in copy2slave:'
                print('Directory not MOVED. Error: %s' % e)
                break
            except OSError as e:
                if e.errno != errno.EEXIST:
                    print 'OSError in copyfiles'
                    break
                    
                else:
                    print 'Error in copy2slave:'
                    print "\nBE CAREFUL! Directory %s already exists." % path
                

        
#==============================================================================
    def delete_files(self,del_slave='off',
                     ext_PEST   =('.1','.2','.3','.4','.5','.6','.7','.9','.bat','.jac','.jco','.jst','.mtt','.par','.prf','.pst','.rec','.rmr','rei','.res','.rst','.stp','.svd','.sen','.seo','.sen'),
                     ext_modflow=('.inf','.spc','.pth','.mpt','.mpn','.nam','.hds','.cbb','.ddn','.glo','.lst','.bas','.dis','.dis','.lpf','.rch','.wel','.oc','.gmg','.riv'),
                     ext_other  = ('.emm','.emo','.err','.dat','.ite','.log','.ref','.tpl','.ins','.ite','.fwr','.tem','.mod','.pts','.dat','.in','.smp','.obs','.vtk','.voo')):             
#==============================================================================
        """
        Module to cleanup HYTEB for "old" output files  
        
        Parameters
        ----------
        del_slave :         str
            delete PEST slave folders
        ext_PEST :          tuple(str)
            tuple with dtype of PEST file(s) to delete
        ext_modflow :       tuple(str)
            tuple with dtype of modflow file(s) to delete
        ext_other :         tuple(str)
            tuple with dtype of AarhusINV file(s) to delete
        
        """

        
        print 'delete_fiels:'        
        cwd = os.getcwd()        
        
        # --- deleting foldes with slaves ---
        if del_slave == 'on':
            print '     deleting all slaves'
            for i in range(100):               
                shutil.rmtree(cwd+'\\'+'slave'+str(i+1),ignore_errors=True)
        if del_slave == 'off':
            print '     no slaves are deleted'
        if del_slave != 'on' and del_slave != 'off':
            print '     wrong input for del_slave!'
         
         
        # --- delete files with extensions ---
        def scandirs(path,exts):
            for root, dirs, files in os.walk(path):
                for currentFile in files:                    
                    if any(currentFile.lower().endswith(ext) for ext in exts):
                        print "deleting file: " + currentFile
                        os.remove(os.path.join(root, currentFile))
        
        scandirs(cwd,ext_other)
        scandirs(cwd,ext_PEST)
        scandirs(cwd,ext_modflow)

        # --- modpath outputfiles to delete ---
        try:
            os.remove('pathline')
        except:
            if IOError:
                print "     no pathline file"
        try:    
            os.remove('endpoint')        
        except:
            if IOError:
                print "     no endpoint file"

#==============================================================================
    def diff_hk_tpl(self,inifile='diff_hk',left=0,right=39,split='REM'):
#==============================================================================
        """
        Writes PEST templates for regularization files
        
        Parameters
        ----------
        inifile :           str
            file of ini-file for calc regularization
        left :              int
            PEST marker, read from this "left" position to "right"
        right :             int
            PEST marker, read to this "left" position from "left"
        split :             str
            string delimiter
        
        """

        print 'diff_hk_tpl:                ',
        
        # --- estimated petrophysical realtionship for all layers ---
        if self.load.petro_est == -1:
            fid1 = open(inifile+'.ini','r')
            text = fid1.readlines()
            fid1.close()
            fid1 = open(inifile+'.tpl','w')
            fid1.write('ptf # \n')   
    
            tex = text[0].split('REM')    
            tmp ='%'+str(left)+'.1s %0.3s %'+str(right-left-3-3)+'.1s %5.'+str(len(tex[1])+5)+'s'        
            fid1.write(tmp % ('#','p_a','# ',split+' '+tex[1]))
    
            tex = text[1].split('REM')        
            tmp ='%'+str(left)+'.1s %0.3s %'+str(right-left-3-3)+'.1s %5.'+str(len(tex[1])+5)+'s'
            fid1.write(tmp % ('#','p_b','# ',split+' '+tex[1]))
            for i in range(2,len(text),1):
                fid1.write(text[i])
                        
            fid1.close()
        
        # --- esimates one petrophysical relationship for each layer --- 
        if self.load.petro_est == -2:
            for layer in range(self.load.n_layer):
                
                fid1 = open(inifile+str(layer+1)+'.ini','r')
                text = fid1.readlines()
                fid1.close()
        
                fid1 = open(inifile+str(layer+1)+'.tpl','w')
                fid1.write('ptf # \n')   
        
                tex = text[0].split('REM')    
                tmp ='%'+str(left)+'.1s %0.9s %'+str(right-left-3-3)+'.1s %5.'+str(len(tex[1])+5)+'s'        
                fid1.write(tmp % ('#','p_a_'+str(layer+1),'# ',split+' '+tex[1]))
        
                tex = text[1].split('REM')        
                tmp ='%'+str(left)+'.1s %0.9s %'+str(right-left-3-3)+'.1s %5.'+str(len(tex[1])+5)+'s'
                fid1.write(tmp % ('#','p_b_'+str(layer+1),'# ',split+' '+tex[1]))
                for i in range(2,len(text),1):
                    fid1.write(text[i])
                            
                fid1.close()
        
        print 'done'   
        
#==============================================================================
    def het2realarray(self,mean_type='A',log_trans=None):
#==============================================================================
        """
        Introduce hetrogenity into zone estimated zones.
        
        mean_type A --> arthemetic mean, G --> geometric mean, H --> harmonic mean
        
        Parameters
        ----------
        mean_type :         str
            A --> arthemetic mean, 
            G --> geometric mean, 
            H --> harmonic mean,
        log_trans :         logical
            True/None
        
        """
        # --- import module to do statiscics ---
        from scipy import stats

        # --- load cat from hydrozones.ini---
        hydro_zones = self.load.hydro_zones(self.load.zone_par_file,self.load.n_zone_hk,4)
        
        # --- load estimated hk parameters ---
        line = [0,1,2]
        par_vec = []
        for i in range(1,self.load.blocksis_n_real+1,1):
            tmp = self.load.PEST_par(self.load.path_results+'\\'+'case.par.'+str(i),self.load,line=line)
            par_vec.append(tmp)
            
            
        # --- multiply estiated parameter values to n zone realizations  ---
        realarray    = np.zeros((self.load.blocksis_n_real,self.load.ny,self.load.nx,self.load.n_layer),dtype=float)
        print 'realizations: '
        for n_real in range(1,self.load.blocksis_n_real+1,1):
            print n_real,
            intarray = self.load.zns2mat(self.load.blocksis_path+'\\'+self.load.blocksis_name+'.'+str(n_real),self.load)
            
            for layer in range(self.load.n_layer):
                for i in range(self.load.ny):
                    for j in range(self.load.nx):
                        for k in range(len(hydro_zones[0:self.load.n_zone_hk,3])):
                            if intarray[i,j,layer] == int(hydro_zones[k,3]):
                                if log_trans == True:
                                    par = np.log10(par_vec[n_real-1][k])
                                elif log_trans == None:
                                    par = par_vec[n_real-1][k]
                                else:
                                    print 'Warning in het2realarray; log_trans != True/None!!!'
                                    break                                          
                                realarray[n_real-1,i,j,layer] = par
        
        # --- calc mean hk fields for n realizations ---
        if mean_type == 'A':
            x1 = np.mean(realarray,axis=0)
        elif mean_type == 'G':
            x1 = stats.gmean(realarray,axis=0)
        elif mean_type == 'H':
            x1 = stats.hmean(realarray,axis=0)
        else:
            print 'Warning in het2realarray; mean_type = A,G or H!!!'
   
        if log_trans == None:
            x2 = x1
        elif log_trans == True:
            x2 = 10**x1
        else:
            print 'Warning in het2realarray; log_trans != True/None!!!'
        
        
        
        # --- write realarray file ---
        # --- return x1,x2,realarray,np.log10(par_vec[n_real-1][k]),intarray
        for layer in range(5):
            fid1 = open(self.load.pp_file_name+str(layer+1)+'.ref','w')
            for i in range(self.load.ny):
                for j in range(self.load.nx):
                    fid1.write('%16.7e' % (x2[i,j,layer]))
                fid1.write('\n')
                
#==============================================================================
    def in_diff_hk(self):
#==============================================================================
        """
        1e-12                     REM a value of the transformation: hk = a*(1/rho)**b.
        4					REM b value.
        hk_L				REM hydraulic conductivity parameter file.
        mod\_				REM electrical resistivity parameter file (mod-file).
        5					REM n header lines in mod-fil (electrical resistivity parameter file).
        model_dis.ini			REM ini-file with model discretization.
        1 2 3 4 5				REM layers to calc hk diff.
        """
        
        print 'in_diff_hk',  

        # 
        if self.load.petro_est == 0 or self.load.petro_est == -1:
           fid1 = open('diff_hk.ini','w') 
           fid1.write(str(self.load.beta_1)+'                 REM a value of the transformation: hk = a*(1/rho)**b'+'\n')
           fid1.write(str(self.load.beta_2)+'                 REM b value'+'\n')
           fid1.write(self.load.model_name_calib+'            REM hydraulic conductivity parameter file'+'\n')
           fid1.write(self.load.mod_file+'                    REM electrical resistivity parameter file (mod-file)'+'\n')
           fid1.write(str(self.load.n_mod_header)+'           REM n header lines in mod-fil (electrical resistivity parameter file)'+'\n')
           fid1.write('model_dis.ini'+'                       REM ini-file with model discretization'+'\n')
           
           for layer in range(self.load.n_layer):
               fid1.write(str(layer+1)+' ')
           fid1.write('                 REM layers to calc hk diff'+'\n')
           
           fid1.close()
        
        
        # --- ini-files for each layer for estimating petro-relation ---
        elif self.load.petro_est == -2:
            for layer in range(self.load.n_layer):
                
                # --- copy diff_hk.py to scripts for calc diff for each layer ---
                tmp1 = self.load.reg_calc.split('.')
                shutil.copy(self.load.reg_calc,tmp1[0]+'_L'+str(layer+1)+'.'+tmp1[1])
                
                # --- write ini-files for diff_hk_l#.ini ---
                fid1 = open('diff_hk_L'+str(layer+1)+'.ini','w') 
                fid1.write(str(self.load.beta_1)+'                  REM a value of the transformation: hk = a*(1/rho)**b'+'\n')
                fid1.write(str(self.load.beta_2)+'                  REM b value'+'\n')
                fid1.write(self.load.model_name_calib+'             REM hydraulic conductivity parameter file'+'\n')
                fid1.write(self.load.mod_file+'                     REM electrical resistivity parameter file (mod-file)'+'\n')
                fid1.write(str(self.load.n_mod_header)+'            REM n header lines in mod-fil (electrical resistivity parameter file)'+'\n')
                fid1.write('model_dis.ini'+'                        REM ini-file with model discretization'+'\n')
                fid1.write(str(layer+1)+' ')
                fid1.write('                 REM layers to calc hk diff'+'\n')
           
                fid1.close()
        print 'done'
    
#==============================================================================
    def in_fac2real(self,grid=1,pp_file_name='',low_lim=1e-20,upper_lim=1e20):         
#==============================================================================
        """
        Writes fac2real ini-files ---
        
        Parameters
        ----------
        grid :              int
            grid number used when running model, see ini file for discretization
        pp_file_name :      Optional(List((str))
            Python list with sufix to fac2real-ini files         
        low_lim :           float
            lower value for interpolating values to grid cells. See PEST doc...
            for more details
        upper_lim :         float
            upper value for interpolating values to grid cells. See PEST doc...
            for more details
        """
        
        for i in range(self.load.n_layer):
            if i < self.load.n_cap:
                fid1 = open('fac2real_L'+str(i+1)+pp_file_name+'.in','w')
                fid1.write('factors_L'+str(i+1)+'.dat'+'\n')
                fid1.write('f'+'\n')
                
                if grid == 1:
                    fid1.write(self.load.model_name_calib+str(i+1)+'.pts'+'\n')
                if grid == 2:                    
                    fid1.write(pp_file_name+str(i+1)+'.pts'+'\n')
                    
                fid1.write('s'+'\n')
                fid1.write(str(low_lim)+'\n')
                fid1.write('s'+'\n')
                fid1.write(str(upper_lim)+'\n')
                if grid == 1:
                    fid1.write(self.load.model_name_calib+str(i+1)+'.ref'+'\n')
                if grid == 2:
                    fid1.write(pp_file_name+str(i+1)+'.ref'+'\n')
                fid1.write('1e35'+'\n')
                fid1.close()
            elif i >= self.load.n_cap:
                fid1 = open('fac2real_L'+str(i+1)+pp_file_name+'.in','w')
                fid1.write('factors_L'+str(i+1)+'.dat'+'\n')
                fid1.write('f'+'\n')
                if grid == 1:
                    fid1.write(self.load.model_name_calib+str(i+1)+'.pts'+'\n')
                if grid == 2:                    
                    fid1.write(pp_file_name+str(i+1)+'.pts'+'\n')
                fid1.write('s'+'\n')
                fid1.write('1e-20'+'\n')
                fid1.write('s'+'\n')
                fid1.write('1e20'+'\n')
                if grid == 1:
                    fid1.write(self.load.model_name_calib+str(i+1)+'.ref'+'\n')
                if grid == 2:
                    fid1.write(pp_file_name+str(i+1)+'.ref'+'\n')
                fid1.write('1e35'+'\n')
                fid1.close()
                
        print('in_fac2real:                done')

#==============================================================================
    def in_ppk2fac(self,output_file,model_spc = 'model.spc',struct_file='struct.ini',structure='structure1',grid=1):
#==============================================================================
        """  
        Writes ini-files for PEST utility ppk2fac.exe        
        
        Parameters
        ----------
        output_file :      str
            String with outputfile name
        model_spc :         str
            pest model specification file
        struct_file :       str
            pest geostatistical structure file        
        structure:          str
            geostatistical structure in struct file
        grid :              int
            pilot point discretization to use, read from model_dis.ini
        
        """       
        
        print 'in_ppk2fac:  ',
        for layer in range(self.load.n_layer):
            if layer < self.load.n_cap:
                fid1 = open(output_file+str(layer+1)+'.in','w')
                fid1.write(model_spc+'\n')
                if grid == 1:
                    fid1.write(self.load.model_name_calib+str(layer+1)+'.pts'+'\n')
                elif grid == 2:
                    fid1.write('a_L'+str(layer+1)+'.pts'+'\n')
                fid1.write('0'+'\n')

                if self.load.parametrization == -1 or self.load.parametrization == -3:
                    tmp = self.load.integer_array_name.split('.')
                    fid1.write(tmp[0]+'_'+str(layer+1)+'.'+tmp[-1]+'\n')
                fid1.write(struct_file+'\n')
                if self.load.parametrization == -1:
                    fid1.write(structure+'\n')
                    fid1.write('o'+'\n')
                    fid1.write('1e20'+'\n')
                    fid1.write('1'+'\n')
                    fid1.write('8'+'\n')
                if self.load.parametrization == -3:
                    for i in range(len(self.load.pp2zone.split(','))+len(self.load.pp2zone_tied.split(','))):                        
                        fid1.write(structure+'\n')
                        fid1.write('o'+'\n')
                        fid1.write('1e20'+'\n')
                        fid1.write('1'+'\n')
                        fid1.write('8'+'\n')
                fid1.write('factors_L'+str(layer+1)+'.dat'+'\n')
                fid1.write('f'+'\n')
                fid1.write('std_L'+str(layer+1)+'.ref'+'\n')
                fid1.write('regularisation_L'+str(layer+1)+'.dat'+'\n')
                fid1.close()      
                
            elif layer >= self.load.n_cap:
                fid1 = open(output_file+str(layer+1)+'.in','w')
                fid1.write(model_spc+'\n')
                if grid == 1:
                    fid1.write(self.load.model_name_calib+str(layer+1)+'.pts'+'\n')
                elif grid == 2:
                    fid1.write('a_L'+str(layer+1)+'.pts'+'\n')
                fid1.write('0'+'\n')

                if self.load.parametrization == -1 or self.load.parametrization == -3:
                    tmp = self.load.integer_array_name.split('.')
                    fid1.write(tmp[0]+'_'+str(layer+1)+'.'+tmp[-1]+'\n')
                fid1.write(struct_file+'\n')
                fid1.write(structure+'\n')
                fid1.write('o'+'\n')
                fid1.write('1e20'+'\n')
                fid1.write('1'+'\n')
                fid1.write('8'+'\n')
                fid1.write('factors_L'+str(layer+1)+'.dat'+'\n')
                fid1.write('f'+'\n')
                fid1.write('std_L'+str(layer+1)+'.ref'+'\n')
                fid1.write('regularisation_L'+str(layer+1)+'.dat'+'\n')
                fid1.close()                

        print '                done'     
        
#==============================================================================
    def ins_files_diff(self,output_file,pp_x,pp_y,x_valley,y_valley,left=39,right=63):       
#==============================================================================
        """
        Writes PEST instructions files for regularization outputfiles
        
        Parameters
        ----------
        output_file :           str
            Outputfilename
        pp_x :              List(int)  
            python list - node positions of pilot points in column direction
        pp_y :              List(int)
            python list - node positions of pilot points in row direction        
        x_valley :          List(int)
            python list - node positions of pilot points in column direction
            inside the buried valley
        y_valley :          List(int)
             python list - node positions of pilot points in row direction
             inside the buried valley
        left :              int
            PEST marker, read from this "left" position to "right" 
        right :             int
            PEST marker, read to this "left" position from "left"
        """         
        
        print 'ins_files_diff:            ',
        
         
        for layer in range(self.load.n_layer):                     
            fid1 = open(output_file+str(layer+1)    +'.ins','w')
            fid1.write('pif #\n')
            
            count1 = 1            
            count2 = 1
            
            if layer < self.load.n_cap:
                for i in range(len(pp_x)*len(pp_y)):
                    fid1.write('l1  [diff_L'+str(layer+1)+'_'+str(count1)+']'+str(left)+':'+str(right)+'\n')
                    count1 += 1
            elif layer >= self.load.n_cap:
                for i in range(len(y_valley[layer-self.load.n_cap])): 
                    for j in range(len(x_valley[layer-self.load.n_cap])): 
                        fid1.write('l1  [diff_L'+str(layer+1)+'_'+str(count2)+']'+str(left)+':'+str(right)+'\n')            
                        count2 += 1
            else:
                print 'Error: ins_files_diff!!!'
                
            fid1.close() 
            
        print 'done'

#==============================================================================
    def ins_files_head(self,output_file,left=40,right=48,obs_type='_1'):
#==============================================================================
        """
        Write PEST instruction file for head observations 
        
        Parameters
        ----------
        output_file :           str
            Outputfilename
        left :              int
            PEST marker, read from this "left" position to "right" 
        right :             int
            PEST marker, read to this "left" position from "left"
        obs_type :          str
            default(1), option to put head observations into seperate groups 
            doing PEST inversion
        """
    
        print 'ins_files_head:           ',
        
        fid1 = open(output_file+'.ins','w')
        fid1.write('pif #'+'\n')
        if obs_type == '_1':
            for i in range(self.load.n_well):
                fid1.write('l1  [well_'+str(i+1)+obs_type+']'+str(left)+':'+str(right)+'\n')
        
        elif obs_type == '_2':
            for i in range(self.load.n_pred_well):
                fid1.write('l1  [well_'+str(i+1)+obs_type+']'+str(left)+':'+str(right)+'\n')
        
        fid1.close()         
        
        print 'done'
         
#==============================================================================
    def ins_files_tem(self,output_file,pp_x,pp_y,n_tem_data,left=13,right=25):         
#==============================================================================
        """
        Write PEST instruction file for AarhusInv datafils 
        
        Parameters
        ----------
        output_file :           str
            Outputfilename
        pp_x :              List(int)  
            python list - node positions of pilot points in column direction
        pp_y :              List(int)
            python list - node positions of pilot points in row direction        
        n_tem_data  :       List(int)
            python list. Index holdes number of gates looped over total 
            number of tem-file.
        left :              int
            PEST marker, read from this "left" position to "right" 
        right :             int
            PEST marker, read to this "left" position from "left"    

        """
        print 'ins_files_tem:             ',
        
        # ---writing TEM-instruction fil---
        count = 0
        self.load
        for i in pp_y:
            for j in pp_x:
                fid1=open(output_file+str(i)+'_'+str(j)+'.ins','w')
                fid1.write('pif #\n')
                fid1.write('l11  [dbdt_'+str(i)+'_'+str(j)+'_1]'+str(left)+':'+str(right)+'\n')
                for h in range(1,n_tem_data[count],1): # number of measurement. See tem-file
                    fid1.write('l1  [dbdt_'+str(i)+'_'+str(j)+'_'+str(h+1)+']'+str(left)+':'+str(right)+'\n')
                count += 1
                fid1.close()
        
        print 'done' 

#==============================================================================
    def ins_file_rchsum(self):               
#==============================================================================
        """
        Write PEST instruction file for rchsum.exe
        Rchsum is a fortran script written by Steen Christensen for calculating
        the aveage travel time and rechage area for particle tracking (MODPATH)
        """
        print 'ins_file_rchsum:             ',
        
        fid1 = open('rchsum_dat.ins','w')
        fid1.write('pif #'+'\n')
        fid1.write('l3  [ave_age]'+str(34)+':'+str(51)+'\n')
        fid1.write('l1  [area]'+str(17)+':'+str(27)+'\n')
        fid1.close()
        
#==============================================================================
    def ins_files_river(self,output_file,left=40,right=55):       
#==============================================================================
        """
        Write PEST instruction file for rchsum.exe
        
        output_file :           str
            Outputfilename
        left :              int
            PEST marker, read from this "left" position to "right" 
        right :             int
            PEST marker, read to this "left" position from "left"            
        """
        
        print 'ins_files_river:             ',        
        
        fid1 = open(output_file+'.ins','w')
        fid1.write('pif #'+'\n')
        fid1.write('l1  [flow]'+str(left)+':'+str(right))
        fid1.close()
        
        print 'done'
        
#==============================================================================
    def int2realarray(self,fname,line = [3,4,5]):
#==============================================================================
        """
        Write real-array file with mean hk vlaues in each cell.
        Reads n realizatin (defind in HYTEB ini-file) and multiply PEST 
        estimated hk values with corresponding zones.
        The estimated hk values is read from PEST parameter file (.par)
        
        Parameters
        ----------
        fname :             str
            Outputfilename
        line :              List(int)
            Line number to read est hk values from. 
        """
        # --- load cat from hydrozones.ini---
        hydro_zones = self.load.hydro_zones(self.load.zone_par_file,self.load.n_zone_hk,4)
        
        # --- load estimated hk parameters ---
        par_vec = []
        for i in range(1,self.load.blocksis_n_real+1,1):
            tmp  = self.load.PEST_par(self.load.path_results+'\\'+'case.par.'+str(i),self.load,line=line)
            par_vec.append(tmp)
            
        # --- allocate numpy array ---
        realarray    = np.zeros((self.load.blocksis_n_real,self.load.ny,self.load.nx,self.load.n_layer),dtype=float)
        
        # --- multiply estiated parameter values to n zone realizations  ---     
        print 'realizations: '
        for n_real in range(1,self.load.blocksis_n_real,1):
            print n_real,
            intarray = self.load.zns2mat(self.load.blocksis_path+'\\'+self.load.blocksis_name+'.'+str(n_real),self.load)
            
            for layer in range(1):
                for i in range(self.load.ny):
                    for j in range(self.load.nx):
                        for k in range(len(hydro_zones[line[0]:line[-1]+1,3])):
                            if intarray[i,j,layer] == int(hydro_zones[k,3]):
                                realarray[n_real-1,i,j,layer] = par_vec[n_real-1][k]*intarray[i,j,layer] 
        
        
        # --- calc mean hk fields for n realizations ---
        x = np.mean(realarray,axis=0)
        
        # ---write realarray file ---
        fid1 = open(fname+'.ref','w')
        for i in range(self.load.ny):
            for j in range(self.load.nx):
                fid1.write('%16.7e' % (x[i,j,layer]))
            fid1.write('\n')
         
#==============================================================================
    def lpf_tpl(self,left=0,right=39):
#==============================================================================
        """
        Writes PEST template file for MODFLOW-lpf file, when groundwater model 
        is parametrized with zones.
        
        Parameters
        ----------
        left :              int
            PEST marker, read from this "left" position to "right" 
        right :             int
            PEST marker, read to this "left" position from "left"
        """        
        
        print 'modflow lpf template file; ',
        
        if self.load.parametrization == -2:
            fid2 = open(self.load.modflow_file_name+'_calib.lpf','r')
            text = fid2.readlines()
            fid2.close()
            
            fid1 = open(self.load.modflow_file_name+'_calib_lpf.tpl','w')
            fid1.write('ptf $'+'\n')
            for i in range(len(text)):
                if i in np.arange(7,200,self.load.n_layer+1) and i < (self.load.n_layer+1)*self.load.n_zone_hk+7:
                    tmp1 = text[i].split()
                    tmp2 = ('%1.'+str(len(tmp1[0]))+'s'+
                            '%'+str(len(tmp1[1])+1)+'.'+str(len(tmp1[1]))+'s'+
                            '%2.1s'+
                            '%'+str(len(tmp1[0])+1)+'.'+str(len(tmp1[0]))+'s'+
                            '%'+str(right-left)+'.1s'+
                            '%'+str(len(tmp1[len(tmp1)-1])+1)+'.'+str(len(tmp1[len(tmp1)-1]))+'s'+
                            '\n')
                    fid1.write(tmp2 % (tmp1[0],tmp1[1],'$',tmp1[0],'$',tmp1[len(tmp1)-1]))
                else:
                    fid1.write(text[i])
            
            fid1.close()
        
        print 'done'


#==============================================================================
    def model2folder(self,model_numb):
#==============================================================================        
        """
        Copy (true) geological reference structures generated with TPRoGS from 
        folder to work directory.
        
        Parameters
        ----------
        model_numb :        int
            geological reference structure number (realization number)
            
        """
        import shutil
        
        cwd = os.getcwd()
        # --- copy recharge file ---
        source  = self.load.model_path+'/'+self.load.recharge+'.'+self.load.dtype_recharge+'.'+str(model_numb)
        dest    = cwd+'/'+self.load.recharge+'.'+self.load.dtype_recharge#+'.'+str(model_numb)
        shutil.copyfile(source,dest)
        
        # --- copy hydraulic conductivity and porosity fields ---
        for i in range(self.load.n_layer):
            # --- hk-fields ---            
            source  = self.load.model_path+'/'+self.load.model_name+str(i+1)+'.'+self.load.dtype_model+'.'+str(model_numb)
            dest    = cwd+'/'+self.load.model_name+str(i+1)+'.'+self.load.dtype_model#+'.'+str(model_numb)
            shutil.copyfile(source,dest)
            # --- por-fields ---
            source  = self.load.model_path+'/'+self.load.por_name+str(i+1)+'.'+self.load.dtype_por+'.'+str(model_numb)
            dest    = cwd+'/'+self.load.por_name+str(i+1)+'.'+self.load.dtype_por#+'.'+str(model_numb)
            shutil.copyfile(source,dest)


#==============================================================================
    def modfile_tpl(self,pp_x,pp_y):      
#==============================================================================
        """
        Write modfile template file for PEST
        
        Parameters
        ----------
        pp_x :              List(int)  
            python list - node positions of pilot points in column direction
        pp_y :              List(int)
            python list - node positions of pilot points in row direction
        
        """
        print 'modfile_tpl:              ',  
        
               
        for j in pp_x:
            for k in pp_y:
                self.write.mod2tpl(j,k,self.load)
        
        print ' done'    
            
#==============================================================================
    def obs_point_in_valley(self,pp_x,pp_y):
#==============================================================================
        """
        Reads the positions of the pilot points and the MODFLOW ibound array.
        Return positions of pilot point in active cells for each layer. 
        
        Parameters
        ----------
        pp_x :              List(int)  
            python list - node positions of pilot points in column direction
        pp_y :              List(int)
            python list - node positions of pilot points in row direction
        
        Returns
        -------
        x_valley :          List(int)
            python list - node positions of pilot points in column direction
            inside the buried valley   
        y_valley :          List(int)
             python list - node positions of pilot points in row direction
             inside the buried valley
        n_pp_valley :       int
            number of active pilot point (inside buried valley)
        """
        
        print 'obs_point_in_valley:       ',        
        
        x_valley       = []
        y_valley       = []
        valley_dis_x   = np.zeros((self.load.n_layer-self.load.n_cap,2),dtype=float)        
        valley_dis_y   = np.zeros((self.load.n_layer-self.load.n_cap,2),dtype=float)       
      
      # --- location of buried valley, checking all layer ---               
        count = 0
        for layer in range(self.load.n_layer+1):
            if layer > (self.load.n_cap):
                fid2 = open(self.load.model_name+str(layer)+'.'+self.load.dtype_model,'r')
                hk = np.array(fid2.read().split(),dtype=float).reshape(self.load.ny,self.load.nx)
                fid2.close()
            
                fid2 = open(self.load.ibound+str(layer)+'.'+self.load.dtype_ibound,'r')
                ibound = np.array(fid2.read().split(),dtype=float).reshape(self.load.ny,self.load.nx)
                ibound = np.absolute(ibound)            
                fid2.close()

                hk = hk*ibound
                x = []
                y = []
                for i in range(self.load.ny):            # north-south "nodes"
                    for j in range(self.load.nx):        # east-west "nodes"
                                if hk[i,j] != 0:
                                    x.append(j)
                                    y.append(i)
                                    
                valley_dis_x[count,:]  = [min(x), max(x)]
                valley_dis_y[count,:]  = [min(y), max(y)]
            
    
                count += 1
                del hk,x,y
                
        # --------------------------------------------------
        #       --- measurement points inside valley ---
        # --------------------------------------------------    
          
        count  = 0
        count1 = 0
        count3  = 0
        for layer in range(self.load.n_layer):
            y =[]        
            
            if layer > self.load.n_cap-1:                
#                for i in range(self.load.n_pp_y):  
                for i in range(len(pp_y)):                                        
                    if  float(valley_dis_y[count,0]) <= float(pp_y[i]) <= float(valley_dis_y[count,1]):
                        y.append((pp_y[i]))
                        count1 += 1
                        #--- flip y_coor updown to fit MODFLOW/PEST input format ---
                        y = sorted(y,reverse = True)
                y_valley.append(y)
                count += 1        


            
            if layer > self.load.n_cap-1:
                
                count2 = 0
                x =[]
#                for j in range(self.load.n_pp_x):
                for j in range(len(pp_x)):
                    if  float(valley_dis_x[count3,0]) <= float(pp_x[j]) <= float(valley_dis_x[count3,1]):                               
                        x.append((pp_x[j]))
                        count2 += 1
                count3 += 1
                x_valley.append(x)
                    
        n_pp_valley = 0                
        for i in range(self.load.n_layer-self.load.n_cap):
            n_pp_valley += len(x_valley[i])*len(y_valley[i])
        
        print 'done'    
        return x_valley, y_valley,n_pp_valley



#==============================================================================
    def PEST_inversion(self,pp_x=[],pp_y=[],x_valley=[],y_valley=[],n_pp_valley=[],
                       n_tem_data=0,res_err=None,IREGADJ=1,SVDMODE=1,
                       PESTMODE='regul',grid=1,pp_file_name=[]):
#==============================================================================
        """
        Module to make pest control file (pcf)
        Creates a PEST control file for inversion --> case.pst

        Parameters
        ----------
        pp_x :              List(int)  
            python list - node positions of pilot points in column direction
        pp_y :              List(int)
            python list - node positions of pilot points in row direction
        x_valley :          List(int)
            python list - node positions of pilot points in column direction
            inside the buried valley
        y_valley :          List(int)
             python list - node positions of pilot points in row direction
             inside the buried valley
        n_pp_valley :       int
            number of active pilot point (inside buried valley)
        n_tem_data  :       List(int)
            python list. Index holdes number of gates looped over total 
            number of tem-file.
        res_err :           array_like(float)
            python list(n,z,m), 
            n = number of TEM soundings
            z = number of layers
            m = 2 (res,STD)
        IREGADJ :           int
            See PEST doc for more info
        SVDMODE :           int
            PEST singular value decomposition mode, See PEST doc for more info
            1 -> on
            0 -> off      
        PESTMODE :          str
            regul       -> for regularization mode for underdetermined problems
            estimation  -> for estimation mode for overdetermined problems
        grid :              int
            grid number used when running model, see ini file for discretization          
        pp_file_name :      Optional(List(str))
            Python list with sufix to fac2real-ini files
        
         
        """
        print 'PEST_joint                 ',
        
        # --- check inversion type ---
        for i in range(1):
            try:
                if self.load.inv_type == 'hydrologic':
                    print 'HYDROLOGIC INVERSION',
                if self.load.inv_type == 'joint':
                    print 'JOINT INVERSION',
                if self.load.inv_type == 'sequential':    
                    print 'SEQUENTIAL INVERSION',
                if self.load.inv_type != 'joint' and self.load.inv_type != 'sequential' and self.load.inv_type != 'hydrologic':
                    print 'wrong input in inversion ini file for inversion type'
                    break
            except:
                print 'something wrong in PEST_pst_control_file'
            
        # --- pre-fix for the the data type    
        gp_dtype            = 'dbdt'   
                                                 
        # --- small_scale_error (struc_noise) is multiplied to the STD when calc head weight for PEST
        small_scale_error   = self.load.struct_noise 
                                   
        #  --- number of tem/pilot point locations
        n_pos               = len(pp_x)*len(pp_y)       

        # --- allocate array for k-field and ibound field --- 
        hk = 9e9*np.ones((self.load.n_layer,self.load.ny,self.load.nx),dtype=float)
        
        ib = np.zeros((self.load.n_layer+1,self.load.ny,self.load.nx),dtype=float)
        ib[0,:,:] = np.ones((self.load.ny,self.load.nx),dtype=float)
        
        for layer in range(self.load.n_layer):            
            fid1 = open(self.load.model_name+str(layer+1)+'.'+self.load.dtype_model,'r')
            fid3 = open(self.load.ibound+str(layer+1)+'.'+self.load.dtype_ibound,'r')
            
            k_field         = np.array(fid1.read().split(),dtype=float).reshape(self.load.ny,self.load.nx)
            i_field         = np.array(fid3.read().split(),dtype=float).reshape(self.load.ny,self.load.nx)
            i_field         = np.absolute(i_field)
            hk[layer,:,:]   = k_field*np.absolute(i_field)
            ib[layer+1,:,:] = np.absolute(i_field)
            fid1.close()
            fid3.close()  
            
        # --------------------------------------------------------------------------
        #   --- calc the numb of parameter held tied under PEST INVERSION ---
        # --------------------------------------------------------------------------
        n_tied_par = 0
        
        # --- joint case, tied electical res parameters outside the "active" part of modflow (1) ...
        # --- However a extra layer is added to estimated the res of the bottom of the model ---
        if self.load.inv_type == 'joint':
            count1 = 1
            count2 = 0
            for i in pp_y:        
                for j in pp_x:
                    if self.load.tied_res == 'false':
                        n_tied_par = 0                        
                    elif self.load.tied_res == 'true':
                        for layer in range(int(self.load.n_layer+1)):
                            if ib[layer,i,j] == 0.0: 
                                count1 += 1                             
                            elif ib[layer,i,j] != 0.0:
                                count2 += 1
            n_tied_par += count1
            
        if self.load.inv_type == 'hydrologic':
            if self.load.parametrization == -1:
                pass
            elif self.load.parametrization == -2:
                try:
                    hz = self.load.hydro_zones(self.load.zone_par_file,self.load.n_zone_hk,4)
                except:
                    print 'Warning_1 hydro_zones.ini file not found or could not be loaded'
            
        if self.load.inv_type == 'sequential':
            if self.load.parametrization == -1:
                pass
            elif self.load.parametrization == -2:
                try:
                    hz = self.load.hydro_zones(self.load.zone_par_file,self.load.n_zone_hk,4)
                except:
                    print 'Warning_1 hydro_zones.ini file not found or could not be loaded'
        
        # --- Combind zones and pilot points; tied pilot point within zones (parameters) excludede from the inversion.ini file ---
        if self.load.parametrization == -3:
            
            n_tied_par = 0
            # --- zone(s) with active pilot point parametrization, see inversion.ini ----
            zones   = map(int,self.load.pp2zone.split(','))
            
            # --- zone(s) with tied parameters (homogen hk-fields) ---
            zones_tied   = map(int,self.load.pp2zone_tied.split(','))
            
            # --- load integer array with zonation ---
            int_mat = self.load.int2mat(self.load.integer_array_name,self.load)
            
            # --- calc number of tied parameters ---
            count1 =0
            for layer in range(self.load.n_layer):            
                for i in pp_y:
                    for j in pp_x:
                            
                            if int_mat[i,j,layer] in zones_tied:
                                count1 += 1
            n_tied_par += count1                                                        

        
        # --------------------------------------------------------------------------
        #              --- parameter for the pest controle file ---
        # --------------------------------------------------------------------------
                
        # --- number of parameters ---
        if self.load.inv_type == 'joint': 
            NPAR    = (self.load.n_layer+1)*len(pp_x)*len(pp_y)+self.load.n_cap*n_pos+n_pp_valley       
        if self.load.inv_type == 'sequential':
            NPAR    = self.load.n_cap*n_pos+n_pp_valley  
        if self.load.inv_type == 'hydrologic':        
            tmp = 0
            if self.load.parametrization == -3:
                NPAR    = self.load.n_cap*n_pos+n_pp_valley  
            if self.load.parametrization == -2:
                NPAR    = self.load.n_zone_hk
            if self.load.parametrization == -1:
                NPAR    = self.load.n_cap*n_pos+n_pp_valley

        if self.load.rch_est == -1:
            if self.load.parametrization == -1:
                NPAR += 2   
            if self.load.parametrization == -2:
                NPAR += self.load.n_zone_rch
            if self.load.parametrization == -3:
                NPAR += 2   
                
        elif self.load.rch_est == -2: 
            NPAR += self.load.n_zone_rch
            
                
        if self.load.petro_est == -1:
            NPAR += 2
        elif self.load.petro_est == -2:
            NPAR += 2*self.load.n_layer
        elif self.load.petro_est == -3:
            NPAR += self.load.n_cap*n_pos+n_pp_valley
            
        # --- number of measurement observations ---                    
        if self.load.inv_type == 'joint': 
            NOBS    = self.load.n_well+1+sum(n_tem_data)+n_pp_valley+n_pos*self.load.n_cap                       # Number of observations:       k+river+rhoa+diff
        if self.load.inv_type == 'sequential':
            if self.load.parametrization == -1:
                NOBS    = self.load.n_well+self.load.n_river+n_pp_valley+n_pos*self.load.n_cap                       # Number of observations:       k+river+diff
        if self.load.inv_type == 'hydrologic':
            if self.load.parametrization != -1: 
                NOBS    = self.load.n_well+self.load.n_river                
        if self.load.inv_type == 'hydrologic':
            if self.load.parametrization == -1:
                NOBS    = self.load.n_well+self.load.n_river                       # Number of observations:       k+river+diff
            if self.load.parametrization != -1: 
                NOBS    = self.load.n_well+self.load.n_river
        # --- add prediction to calibration dataset ---
        if self.load.pred2calibdata < 0:
            if self.load.pred2calibdata == -1 or self.load.pred2calibdata == -3:
                NOBS += 2
            if self.load.pred2calibdata == -1 or self.load.pred2calibdata == -2:
                NOBS += self.load.n_pred_well
                
        # --- number of parameter groups ---        
        if self.load.inv_type == 'joint':             
            NPARGP  = 2
        if self.load.inv_type == 'sequential':
            NPARGP  = 1
        if self.load.inv_type == 'hydrologic': 
            if grid == 1:
                NPARGP = 1
            else:
                NPARGP = 0
            
        if self.load.parametrization == -1:    
            if self.load.rch_est == -1:
                NPARGP += 2                            
        if self.load.parametrization == -2:    
            if self.load.rch_est == -2:
                NPARGP += 1
        if self.load.parametrization == -3:    
            if self.load.rch_est == -1:
                NPARGP += 2
            elif self.load.rch_est == -2:  
                NPARGP += self.load.n_zone_rch
                
        if self.load.petro_est < 0:
            NPARGP += 2     
            
        # --- number of observations groups ---
        if self.load.parametrization == -1: 
            if self.load.inv_type == 'joint':
                if self.load.regulasobs == None:             
                    NOBSGP  = 3 + self.load.n_layer
                else:
                    NOBSGP  = 4
            elif self.load.inv_type == 'hydrologic':             
                    NOBSGP  = 2 
            elif self.load.inv_type == 'sequential':
                if self.load.regulasobs == True:
                    NOBSGP = 3
                else:
                    NOBSGP  = 2
            # --- add prediction to calibration dataset ---
            if self.load.pred2calibdata < 0:
                if self.load.pred2calibdata == -1 or self.load.pred2calibdata == -3:
                    NOBSGP += 2
                if self.load.pred2calibdata == -1 or self.load.pred2calibdata == -2:
                    NOBSGP += 1
           
                                                                    # Number of observation groups: k+river+rho+regul       
        if self.load.parametrization == -2:                                                             # Number of observation groups: k+river+regul
            NOBSGP = 2
        if self.load.parametrization == -3:                                                             # Number of observation groups: k+river+regul
            NOBSGP = 2
        
        
        # --- number of template files ---
        if self.load.parametrization == -1:
            if self.load.inv_type == 'joint':
                NTPLFLE = n_pos+self.load.n_layer 
            if self.load.inv_type == 'sequential':                                                      
                NTPLFLE = self.load.n_layer
            if self.load.inv_type == 'hydrologic':
                NTPLFLE = self.load.n_layer
        
        if self.load.parametrization == -2:
            if self.load.inv_type == 'hydrologic':                                                      
                NTPLFLE = 1 
        if self.load.parametrization == -3:
                NTPLFLE = self.load.n_layer 
                
        if self.load.rch_est != 0:
            NTPLFLE += 1
        if self.load.petro_est == -1:
            NTPLFLE += 1
        elif self.load.petro_est == -2: 
            NTPLFLE += self.load.n_layer        
        elif self.load.petro_est == -3:    
            NTPLFLE += self.load.n_layer  
            
        #==============================================================================
        # --- number of instrustion files ---
        #==============================================================================
        if self.load.parametrization == -1:
            if self.load.inv_type == 'joint':
                NINSFLE = n_pos+self.load.n_layer+2                                                     # Number of instruction files:          
            if self.load.inv_type == 'sequential':                                                      
                NINSFLE = self.load.n_layer+2
            if self.load.inv_type == 'hydrologic':                                                      
                NINSFLE = 2    
                                                 # Number of instruction files:          
        if self.load.parametrization == -2:
            if self.load.inv_type == 'hydrologic':                                                      
                NINSFLE = 2
        if self.load.parametrization ==-3:
            if self.load.inv_type == 'hydrologic':                                                      
                NINSFLE = 2
        # --- add prediction to calibration dataset ---
        if self.load.pred2calibdata < 0:
            if self.load.pred2calibdata == -1 or self.load.pred2calibdata == -3:
                NINSFLE += 1
            if self.load.pred2calibdata == -1 or self.load.pred2calibdata == -2:
                NINSFLE += 1
                
        # --- number of command lines ---
        if self.load.inv_type == 'joint':
            NUMCOM  = 2+n_pos                                                                       # Number of command lines used to run model
        if self.load.inv_type == 'sequential':                                                      
            NUMCOM  = 1                                       
        if self.load.inv_type == 'hydrologic':                                # Number of command lines used to run model
            NUMCOM  = 1  
        
        # --- command line number to run groundwater model  ---
        if self.load.inv_type == 'joint':
            hk_num  = 1+1+n_pos                                                                     # command line to run groundwater model
        if self.load.inv_type == 'sequential':                                                      
            hk_num  = 1                                                                     # command line to run groundwater model
        if self.load.inv_type == 'hydrologic':                                                      
            hk_num  = 1     
        
        
        RLAMBDA1= 10                                                                            # Initial Marquardt Lambda
        RLAMFAC = -3                                                                            # Dictates Marquardt Lambda adjustment process
        NOPTMAX = self.load.PEST_niter                                                          # max number of iterations (-1, jacobi is calculated)
        MAXSING = NPAR-n_tied_par                                                               #  NUMBER of parameters for SVD
        FRACPHIM= 0.1
        
        
        if self.load.inv_type == 'joint':
            PHILIM  = 3
        if self.load.inv_type == 'sequential':
            PHILIM  = 2
            if self.load.regulasobs == True:
                PHILIM  += 1                
        if self.load.inv_type == 'hydrologic':
            PHILIM  = 2
            
        #-----------------------------------------------------
        # start of pcf                      <--------------
        #-----------------------------------------------------
        fid2 = open('case.pst','w')        
        fid2.write('pcf\n')
        fid2.write('* control data\n')
        
        if self.load.inv_type != 'hydrologic' or self.load.regulasobs == True:
            if PESTMODE == 'regul':
                fid2.write('restart  regularization \n')
            elif PESTMODE == 'est':
                fid2.write('restart  estimation \n')
                
        elif self.load.inv_type == 'hydrologic' or self.load.regulasobs == None:
            fid2.write('restart  estimation \n')
                    
        fid2.write(str(NPAR)+'      '+str(NOBS)+'       '+str(NPARGP)+'       0       '+str(NOBSGP)+' \n')
        fid2.write(str(NTPLFLE)+'     '+str(NINSFLE)+' single  point  '+str(NUMCOM)+' 0  0 '+'\n')
        fid2.write(str(RLAMBDA1)+'  '+str(RLAMFAC)+'  0.3  0.01  10 999 ')
        if self.load.PEST_lamforgive == 'nolamforgive':        
            fid2.write(self.load.PEST_lamforgive+' ')
        if self.load.PEST_lamforgive == 'lamforgive':        
            fid2.write(self.load.PEST_lamforgive+' ')            
        else:
            print 'something wrong with lamforgive in inversion_setup.ini'
            
        if self.load.PEST_derforgive == 'noderforgive':        
            fid2.write(self.load.PEST_derforgive+'\n')
        if self.load.PEST_derforgive == 'derforgive':        
            fid2.write(self.load.PEST_derforgive+'\n')            
        else:
            print 'something wrong with derforgive in inversion_setup.ini'            
        fid2.write('10.0  10.0  0.001 \n')
        fid2.write('0.1  noaui \n')
        fid2.write(str(NOPTMAX)+'  0.005  4  3  0.01  4  0  1\n')
        fid2.write('0  0  0  0  \n')
        
        if self.load.parametrization != -2:
            # --- SVD / singular value decomposition ---
            fid2.write('* singular value decomposition'+'\n')
            fid2.write(str(SVDMODE)+'\n')
            fid2.write(str(MAXSING)+'   5e-7'+'\n')
            fid2.write('0'+'\n')

        #-------------------------------------------------------
        # --- parameter groups ---
        #-------------------------------------------------------
        fid2.write('* parameter groups\n')
        if self.load.parametrization != -2:
            if self.load.rch_est == -2:                                       # <-------------------------------
                for i in range(1,self.load.n_zone_rch+1,1):
                    fid2.write(' rch_'+str(i)+'      relative    0.01  0.0	switch  2.0 parabolic '+'\n')
                
            elif self.load.rch_est == -1:
                fid2.write(' rch_a      relative    0.01  0.0	switch  2.0 parabolic '+'\n')
                fid2.write(' rch_b      relative    0.01  0.0	switch  2.0 parabolic '+'\n')
        
        if self.load.rch_est != 0 and self.load.parametrization == -2:
            fid2.write(' rch      relative    0.01  0.0	switch  2.0 parabolic '+'\n')
            
        if self.load.petro_est < 0:
            fid2.write(' petro_a      relative    0.01  0.0	switch  2.0 parabolic '+'\n')
            fid2.write(' petro_b      relative    0.01  0.0	switch  2.0 parabolic '+'\n')            
        
        if self.load.inv_type == 'joint':
            fid2.write(' h          relative    0.01  0.0	switch  2.0 parabolic '+'\n')
            fid2.write(' res        relative    0.01  0.0	switch  2.0 parabolic '+'\n')
            
        if self.load.inv_type == 'sequential':                                                      
            fid2.write(' h          relative    0.01  0.0	switch  2.0 parabolic '+'\n')
            
        if self.load.inv_type == 'hydrologic':
            if grid == 1:
                fid2.write(' h          relative    0.01  0.0	switch  2.0 parabolic '+'\n')
          
                    
        #-------------------------------------------------------
        # --- parameter data / STARTING MODEL --- 
        #-------------------------------------------------------
        fid2.write('* parameter data \n')
        if self.load.parametrization != -2:
            # --- using log-log linear realtionship for transformating hk to rch ---
            if self.load.rch_est == -1:
                fid2.write('%1.15s %10.9s %10.9s %20.6e %15.3e %15.3e %10.5s %10.3f %10.3f %10.0f'%('recharge_a','log','factor',4.22e-08,1e-9,1e-6,'rch_a',1.0,0.0,hk_num)+'\n')
                fid2.write('%1.15s %10.9s %10.9s %20.6e %15.3e %15.3e %10.5s %10.3f %10.3f %10.0f'%('recharge_b','log','factor',1.605,0.5,2.5,'rch_b',1.0,0.0,hk_num)+'\n')
            # --- n recharge zone parameters defined in hydro_zones.ini --- 
            elif self.load.rch_est == -2:
                rch = self.load.hydro_zones('dummy.ini',self.load.n_zone_rch,4)
                for i in range(1,self.load.n_zone_rch+1,1):
                    fid2.write('%1.15s %10.9s %10.9s %20.6e %15.3e %15.3e %10.5s %10.3f %10.3f %10.0f'%('recharge_'+str(i),'log','factor',float(rch[-self.load.n_zone_rch-1+i,2]),1e-10,1e-7,'rch_'+str(i),1.0,0.0,hk_num)+'\n')
           
            if self.load.petro_est == -1:
                fid2.write('%1.15s %10.9s %10.9s %20.6e %15.3e %15.3e %10.7s %10.3f %10.3f %10.0f'%('p_a','log','factor',1e-12,1e-20,1e-5,'petro_a',1.0,0.0,hk_num)+'\n')
                fid2.write('%1.15s %10.9s %10.9s %20.6e %15.3e %15.3e %10.7s %10.3f %10.3f %10.0f'%('p_b','log','factor',4,1e-3,100,'petro_b',1.0,0.0,hk_num)+'\n')
            elif self.load.petro_est == -2:
                for layer in range(self.load.n_layer):
                    fid2.write('%1.15s %10.9s %10.9s %20.6e %15.3e %15.3e %10.7s %10.3f %10.3f %10.0f'%('p_a_'+str(layer+1),'log','factor',1e-12,1e-20,1e-5,'petro_a',1.0,0.0,hk_num)+'\n')
                    fid2.write('%1.15s %10.9s %10.9s %20.6e %15.3e %15.3e %10.7s %10.3f %10.3f %10.0f'%('p_b_'+str(layer+1),'log','factor',4,1e-3,100,'petro_b',1.0,0.0,hk_num)+'\n')

        
        #-------------------------------------------------------
        # --- hydaulic conductivities ---        
        #-------------------------------------------------------
        p_type = 'h'
        if self.load.parametrization == -3:
            tick = np.zeros(len(zones_tied))
        
        if self.load.parametrization == -1 or self.load.parametrization == -3:
            for layer in range(self.load.n_layer):
                count3 = 0
                # --- pilot point discretization ---
                if grid == 1:
                    pp_file_name = [self.load.pp_file_name]
                    ptype = 'h'
                    parlbnd = [1e-15]
                    parubnd = [1e-1]
                    
                if grid == 2:
                    pp_file_name = pp_file_name
                    ptype =['petro_a','petro_b']
                    parlbnd = [1e-15,1]
                    parubnd = [1e-9,7]
                
                 
                for file_name in pp_file_name:
                    count1 = 0
                    count2 = 0
                                        
                    # --- pilot point name ---
                    tmp     = file_name.split('_')
                    pp_name = tmp[0]
                    
                    # - parameter type ---
                    p_type =ptype[count3]
                    
                    # --- Prameter bounds---
                    PARLBND = parlbnd[count3]
                    PARUBND = parubnd[count3]
                    
                    count3 += 1
                    
                    if layer < self.load.n_cap:
                        # --- read pilot point file for initial value of hk parameters ---
                        fid1 = open(file_name+str(layer+1)+'.pts')
                        text = fid1.readlines()
                        fid1.close()
                        ppoints = len(text)
                        for h in range(ppoints):
                            line = text[h].split()
                            j    = int((float(line[2])-self.load.dy/2.0)/self.load.dy)
                            i    = int((float(line[1])-self.load.dx/2.0)/self.load.dx)
                            zone = (int(line[3]))                   
                            cond = (float(line[4]))

                            
                            count1 += 1 
                            
                            if cond == 0.0:
                                print 'WRONG hydraulic con value!',layer, i, j
                                
                            elif cond > 0.0:
                                # --- combined zone and pilot point parametrization ---
                                if self.load.parametrization == -3:
                                    
                                    for cat in zones:                                
                                        if zone == cat:                                   
                                            fid2.write('%1.15s %10.9s %10.9s %20.6e %15.3e %15.3e %10.7s %10.3f %10.3f %10.0f'%(pp_name+str(layer+1)+'_ppt'+str(count1),'log','factor',cond,PARLBND,PARUBND,p_type,1.0,0.0,hk_num)+'\n')
                                                              
                                    for i in range(len(zones_tied)):                                
                                        if zone == zones_tied[i]:
                                            if tick[i] == 0:
                                                fid2.write('%1.15s %10.9s %10.9s %20.6e %15.3e %15.3e %10.7s %10.3f %10.3f %10.0f'%(pp_name+str(layer+1)+'_ppt'+str(count1),'log','factor',cond,PARLBND,PARUBND,p_type,1.0,0.0,hk_num)+'\n')
                                                tick[i] += 1
                                            elif tick[i] > 0:
                                                fid2.write('%1.15s %10.9s %10.9s %20.6e %15.3e %15.3e %10.7s %10.3f %10.3f %10.0f'%(pp_name+str(layer+1)+'_ppt'+str(count1),'tied','factor',cond,PARLBND,PARUBND,p_type,1.0,0.0,hk_num)+'\n')
                                                
                                # --- pilot point parametrization ---            
                                else:
                                    fid2.write('%1.15s %10.9s %10.9s %20.6e %15.3e %15.3e %10.7s %10.3f %10.3f %10.0f'%(pp_name+str(layer+1)+'_ppt'+str(count1),'log','factor',cond,PARLBND,PARUBND,p_type,1.0,0.0,hk_num)+'\n')
                                           
                    elif layer >= self.load.n_cap:
                        fid1 = open(file_name+str(layer+1)+'.pts')
                        text = fid1.readlines()
                        fid1.close()
                        ppoints = len(text)
                        for h in range(ppoints):
                            line = text[h].split()
                            j = int((float(line[2])-self.load.dy/2.0)/self.load.dy)
                            i = int((float(line[1])-self.load.dx/2.0)/self.load.dx)
                            cond = (float(line[4]))
                            count2 += 1
                            fid2.write('%1.15s %10.9s %10.9s %20.6e %15.3e %15.3e %10.7s %10.3f %10.3f %10.0f'%(pp_name+str(layer+1)+'_ppt'+str(count2),'log','factor',cond,PARLBND,1.0e2,p_type,1.0,0.0,hk_num)+'\n')
                    else:
                        print 'layer out of bound!!!'
        
        # -----------------------------------------------------------------            
        if self.load.parametrization == -2: 
                              
            for i in range(self.load.n_zone_hk):
                fid2.write('%1.15s %10.9s %10.9s %20.6e %15.3e %15.3e %10.7s %10.3f %10.3f %10.0f'%(str(hz[i,0])+'_'+str(hz[i,1]),'log','factor',float(hz[i,2]),1e-12,1e-1,p_type,1.0,0.0,hk_num)+'\n')
            if self.load.rch_est == -2:
                for i in range(self.load.n_zone_rch):
                    fid2.write('%1.15s %10.9s %10.9s %20.6e %15.3e %15.3e %10.7s %10.3f %10.3f %10.0f'%(str(hz[2*self.load.n_zone_hk+i,0])+'_'+str(hz[2*self.load.n_zone_hk+i,1]),'log','factor',float(hz[2*self.load.n_zone_hk+i,2]),1e-11,1e-7,'rch',1.0,0.0,hk_num)+'\n')
                      
                
        #-------------------------------------------------------                    
        # --- electrical resistivities ---
        #-------------------------------------------------------
        if self.load.inv_type == 'joint':
            count1 = 1
            count2 = 0
            for i in pp_y:        
                for j in pp_x:
                    count1 += 1
                    # --- read Aarhus Inv mod-file for resistivity initial parameter values ---
                    fid1 = open(self.load.mod_file+str(i)+'_'+str(j)+'.mod','r')
                    text = fid1.readlines()
                    if self.load.tied_res == 'false':
                        for layer in range(int(self.load.n_layer+1)):
                            temp = text[5+layer].split()
                            res = float(temp[0])                   
                            fid2.write('%1.15s %10.9s %10.9s %20.6f %15.3e %15.3e %10.3s %10.3f %10.3f %10.0i'%('r_'+str(i)+'_'+str(j)+'_'+str(layer+1),'log','factor',res,1.0e-1,2.0e4,'res',1.0,0.0,count1)+'\n')
                            
                    elif self.load.tied_res == 'true':
                        for layer in range(int(self.load.n_layer+1)):
                            if ib[layer,i,j] == 0.0: 
                                # makes res_parameter tied in PEST inversion 
                                temp = text[5+layer].split()
                                res = float(temp[0])                                                  
                                fid2.write('%1.15s %10.9s %10.9s %20.6f %15.3e %15.3e %10.3s %10.3f %10.3f %10.0i'%('r_'+str(i)+'_'+str(j)+'_'+str(layer+1),'tied','factor',res,1.0e-1,2.0e4,'res',1.0,0.0,count1)+'\n')
                                count2 += 1 
                            elif ib[layer,i,j] != 0.0:
                                temp = text[5+layer].split()
                                res = float(temp[0])
                                fid2.write('%1.15s %10.9s %10.9s %20.6f %15.3e %15.3e %10.3s %10.3f %10.3f %10.0i'%('r_'+str(i)+'_'+str(j)+'_'+str(layer+1),'log','factor',res,1.0e-1,2.0e4,'res',1.0,0.0,count1)+'\n')
                    fid1.close()
                       
            #-------------------------------------------------------
            # --- tied parameters for joint inversion --
            #-------------------------------------------------------                                       

            # --- electrical resistivity parameters---
            if self.load.tied_res == 'true':
                
                for i in pp_y:
                    for j in pp_x:
                        count1 = 0
                        for layer in range(self.load.n_layer+1):
                            
                            if ib[layer,i,j] == 0:
                                
                                for tied_layer in range(layer+1,self.load.n_layer+2,1):
                                    fid2.write('r_'+str(i)+'_'+str(j)+'_'+str(tied_layer)+'\t')
                                    fid2.write('r_'+str(i)+'_'+str(j)+'_'+str(layer)+'\n')
                                break
    
                            else:
                                continue                                
                            count1 += 1
            
        #-------------------------------------------------------
        # --- tied parmeters for hydraulic/sequential inversion with zonation ---
        #-------------------------------------------------------
        
        # --- hydraulic conductivity parameters ---
        if self.load.parametrization == -3:
            tick = np.zeros(len(zones_tied))
            for k in range(len(zones_tied)):
                for layer in range(self.load.n_layer):
                    fid1 = open(self.load.pp_file_name+str(layer+1)+'.pts')
                    text = fid1.readlines()
                    fid1.close()
                    ppoints = len(text)
                    
                    count1 = 0
                    for h in range(ppoints):
                        line = text[h].split()
                        j    = int((float(line[2])-self.load.dy/2.0)/self.load.dy)
                        i    = int((float(line[1])-self.load.dx/2.0)/self.load.dx)
                        zone = (int(line[3]))                   
                        cond = (float(line[4]))                   
                        count1 += 1    
                        
                                                    
                        if zone == zones_tied[k]:                                 
                            if tick[k] == 0:
#                                fid2.write('%1.15s %10.9s'%('hk'+str(layer+1)+'_ppt'+str(count1),'cat'+str(zone))+'\t')
                                t1 = layer+1
                                t2 = count1          
                                tick[k] += 1
                                
                            elif tick[k] > 0:
                                fid2.write('%1.15s '%('hk'+str(layer+1)+'_ppt'+str(count1))+'\t') 
                                fid2.write('%1.15s '%('hk'+str(t1)+'_ppt'+str(t2),)+'\n')
                                                
        #-------------------------------------------------------
        # --- observation groups ----
        #-------------------------------------------------------
        fid2.write('* observation groups'+'\n')
#        print self.load.regulasobs
        if self.load.inv_type != 'hydrologic' and self.load.regulasobs == None:   
                for i in range(self.load.n_layer):
                    fid2.write('regul_'+str(i+1)+'\n')
                
        
        if self.load.regulasobs == True:
            fid2.write('reg'+'\n')
        fid2.write('heads'+'\n')
        fid2.write('flux'+'\n')
        if self.load.inv_type == 'joint':
            fid2.write('db'+'\n')
        
        # --- add prediction to calibration dataset ---
        if self.load.pred2calibdata < 0:
            if self.load.pred2calibdata == -1 or self.load.pred2calibdata == -3:
                fid2.write('avg_age'+'\n')
                fid2.write('area'+'\n')
            if self.load.pred2calibdata == -1 or self.load.pred2calibdata == -2:
                fid2.write('pred'+'\n')
            
        
        
        #-------------------------------------------------------
        #               --- observation data ---
        #-------------------------------------------------------
        
        fid2.write('* observation data'+'\n')

        # --- head data and weight ---
        if self.load.std_well == 0.0:
            ww          = 1.0/(np.sqrt(self.load.n_well)*0.1*small_scale_error)
        else:
            ww          = 1.0/(np.sqrt(self.load.n_well)*self.load.std_well*small_scale_error)
        
        obsfile = self.load.head_obs_points.split('.')
        fid1 = open(obsfile[0]+'.obs','r')
        count = 0        
        while count < self.load.n_well:
            tmp = fid1.readline().split()  
            fid2.write('%2.15s %20.9s %20.3e %9.7s' % (tmp[0]+'_1',tmp[3],ww,'heads')+'\n')                
            count += 1
        fid1.close()
        
        # --- river flux data ---     
        fid1 = open('bud2smp_river'+'.obs','r')
        count = 0        
        while count < self.load.n_river:
            tmp = fid1.readline().split()
            if self.load.noise_river == 0.0:
                wr          = 1.0/(np.sqrt(self.load.n_river)*0.1*abs(float(tmp[3])))
            else:
                wr          = 1.0/(np.sqrt(self.load.n_river)*self.load.noise_river*abs(float(tmp[3])))
            fid2.write('%2.30s %30.30s %15.3e %8.7s' % (tmp[0],tmp[3],wr,'flux')+'\n')                
            count += 1
        
        fid1.close()
        

        if self.load.inv_type == 'joint':
            # --- Geophysical data ---
            tal2 = 0 
            for k in pp_y:
                for j in pp_x:
                    fid3 = open(self.load.tem_name+str(k)+'_'+str(j)+'.tem','r')
                    fid4 = open('fwr_'+str(k)+'_'+str(j)+'.ins','r')
                    tem_lines = len(fid4.readlines())
                    fid4.close()
                    
                    for i in range(int(tem_lines+self.load.n_data_header)):
                            if i < self.load.n_data_header:
                                text = fid3.readline()                   
                            if i > self.load.n_data_header:                
                                tex = fid3.readline().split()
                                wr = 1/(float(tex[1])*float(tex[2]))                            # weight on data
                                fid2.write('%2.15s %15.12s %20.4e %9.6s'% ('dbdt_'+str(k)+'_'+str(j)+'_'+str(int(i-self.load.n_data_header)),str(tex[1]),wr/np.sqrt(len(pp_y)*len(pp_x)),'db')+'\n')
                                tal2 += 1
                    fid3.close()  
                    
        # ---------------------------------------------------------------------
        # --- the regulratization constain between hk and res as a part of the observation data [log(k_mf)-log(k_Tem)]  ---                    
        # ---------------------------------------------------------------------
        if self.load.inv_type == 'sequential' and self.load.regulasobs == True:
            # --- only applied for sequential ---
            for layer in range(self.load.n_layer):
                count1 = 0
                count2 = 0
                
                if layer < self.load.n_cap:
                    for i in pp_y:
                        for j in pp_x:
                            count1 += 1
                            if self.load.inv_type == 'sequential' or self.load.inv_type == 'joint' and self.load.regulasobs == True:   
                                if res_err == None:
                                    tmp = 1
                                else:
                                    [ix,iy,iz] = np.shape(res_err)
                                    
                                    
                                    # --- unc calculated with the probagation of error, see Barlow eq 4.10 ---
#                                    unc = self.load.beta_2*res_err[count1+count2-1][layer][1]/(res_err[count1+count2-1][layer][0]*mat.log(10))
                                    # --- unc calculated with the probagation of error, see Barlow eq 4.8a ---
                                    unc = self.load.beta_2*res_err[count1+count2-1][layer][1]
                                    # --- weight on prior information, normalized with number of prior observations ---
                                    n_reg = np.sqrt(len(pp_y)*len(pp_x))                                    
                                    w   = 1.0/(n_reg*unc)
                                    
                                fid2.write('diff_L'+str(layer+1)+'_'+str(count1)+'        0.0000000000             '+str(w)+'        reg'+'\n')
                                
                elif layer >= self.load.n_cap:          
                    for h in range(len(y_valley[layer-self.load.n_cap])):
                        for k in range(len(x_valley[layer-self.load.n_cap])):
                            count2  += 1
                            if res_err != None:
                                
                                # --- unc calculated with the probagation of error, see Barlow eq 4.10 ---
                                unc = self.load.beta_2*res_err[count1+count2-1][layer][1]/(res_err[count1+count2-1][layer][0]*mat.log(10))
                                # --- weight on prior information, normalized with number of prior observations ---
                                n_reg = np.sqrt(len(pp_y)*len(pp_x))                                
                                w   = 1.0/(n_reg*unc)
                            
                            elif res_err == None:
                                w = 1

                            if self.load.inv_type == 'sequential' and self.load.regulasobs == True: 
                                fid2.write('diff_L'+str(layer+1)+'_'+str(count2)+'        0.0000000000      '+str(w)+'        reg'+'\n')
                            
                        
                else:
                    print 'layer out of bound!!!'
            
        #-------------------------------------------------------
        #           --- regularization observation --- 
        #-------------------------------------------------------
        
        if self.load.inv_type != 'hydrologic':
            # only applied for sequential and joint inversion
            for layer in range(self.load.n_layer):
                count1 = 0
                count2 = 0
                if layer < self.load.n_cap:
                    for i in pp_y:
                        for j in pp_x:
                            count1 += 1
                            if self.load.inv_type == 'joint' and self.load.regulasobs == None:
                                if res_err == None:
                                    w = 1
                                else:
                                    # --- unc calculated with the probagation of error, see Barlow eq 4.10 ---
                                    unc = self.load.beta_2*res_err[count1+count2-1][layer][1]/(res_err[count1+count2-1][layer][0]*mat.log(10))
                                    # --- weight on prior information, normalized with number of prior observations ---
                                    n_reg = np.sqrt(len(pp_y)*len(pp_x))
                                    w   = 1.0/(n_reg*unc)
                                fid2.write('diff_L'+str(layer+1)+'_'+str(count1)+'        0.0000000000      '+str(w)+'        regul_'+str(layer+1)+'\n')
                            
                            if self.load.inv_type == 'sequential' and self.load.regulasobs == None:   
                                if res_err == None:
                                    w = 1
                                else:
                                    # --- unc calculated with the probagation of error, see Barlow eq 4.10 ---
                                    unc = self.load.beta_2*res_err[count1+count2-1][layer][1]/(res_err[count1+count2-1][layer][0]*mat.log(10))
                                    # --- weight on prior information, normalized with number of prior observations ---
                                    n_reg = np.sqrt(len(pp_y)*len(pp_x))
                                    w   = 1.0/(n_reg*unc)
                                fid2.write('diff_L'+str(layer+1)+'_'+str(count1)+'        0.0000000000      '+str(w)+'        regul_'+str(layer+1)+'\n')
                                
                elif layer >= self.load.n_cap:          
                    for h in range(len(y_valley[layer-self.load.n_cap])):
                        for k in range(len(x_valley[layer-self.load.n_cap])):
                            count2  += 1
                            if res_err != None:
                                tmp = res_err[count1+count2-1][layer]
#                                tmp = 1+res_err[count1+count2-1][layer]
                            elif res_err == None:
                                tmp = 1
                                
                            if self.load.inv_type == 'joint':
                                fid2.write('diff_L'+str(layer+1)+'_'+str(count2)+'        0.0000000000      '+str(tmp)+'        regul_'+str(layer+1)+'\n')
                            if self.load.inv_type == 'sequential' and self.load.regulasobs == None: 
                                fid2.write('diff_L'+str(layer+1)+'_'+str(count2)+'        0.0000000000      '+str(tmp)+'        regul_'+str(layer+1)+'\n')
                            
                        
                else:
                    print 'layer out of bound!!!'
        
        # --- prediction as observations ---
        if self.load.pred2calibdata < 0:
            predfile = self.load.head_pred_points.split('.')
            
            
            
            # --- add head predictions ---
            if self.load.pred2calibdata == -1 or self.load.pred2calibdata == -2:
                fid1 = open(predfile[0]+'.obs','r')
                count = 0        
                ww = 1.0/(np.sqrt(self.load.n_pred_well)*float(tmp[3])*0.1)
                while count < self.load.n_pred_well:
                    tmp = fid1.readline().split()  
                    fid2.write('%2.15s %20.9s %20.3e %9.7s' % (tmp[0]+'_2',tmp[3],ww,'pred')+'\n')                
                    count += 1 
                fid1.close()
                
            # --- add particle tracing predictions ---    
            if self.load.pred2calibdata == -1 or self.load.pred2calibdata == -3:
                fid1 = open('rchsum_true.dat','r')
                data = np.array(fid1.read().split())
                fid1.close()
                ave_age  = data[14]
                area     = data[17]
                ww_area = 1.0/(float(area)*0.2)/np.sqrt(56000)
                ww_age  = 1.0/(float(ave_age)*0.2)/np.sqrt(56000)
                fid2.write('%2.15s %20.9s %20.3e %9.7s' % ('area',area,ww_area,'area')+'\n')  
                fid2.write('%2.15s %20.9s %20.3e %9.7s' % ('ave_age',ave_age,ww_age,'avg_age')+'\n')  
#                fid2.write('%2.15s %20.9s %20.3e %9.7s' % ('area',area,ww,'area')+'\n')  
                
                
                                 
                
            
        #-------------------------------------------------------
        #              --- model command line ---
        #-------------------------------------------------------
        fid2.write('* model command line'+'\n')
        if self.load.inv_type == 'joint':
            fid2.write('allmodels.bat'+'\n')
            for k in pp_y:
                for j in pp_x:
                    fid2.write('gp_model.bat'+' '+str(k)+' '+str(j)+'\n')
        fid2.write('gw_model.bat'+'\n')
        
        # --- model input/output ---
        fid2.write('* model input/output\n')
        
        # --- recharge template file ---
        if self.load.rch_est != 0:
            tmp = self.load.rch_est_calc.split('.')
            if self.load.parametrization == -1 or self.load.parametrization == -3:
                fid2.write(tmp[0]+'.tpl'+'\t'+tmp[0]+'.ini'+'\n')
            if self.load.parametrization == -2:
                fid2.write('recharge_est.tpl'+'\t'+self.load.modflow_file_name+'_calib.rch'+'\n')
                
        # --- petro template file ---
        if self.load.petro_est == -1:
            fid2.write('diff_hk.tpl'+'\t'+'diff_hk.ini'+'\n')
        if self.load.petro_est == -2:
            for layer in range(self.load.n_layer):
                fid2.write('diff_hk_L'+str(layer+1)+'.tpl'+'\t'+'diff_hk_L'+str(layer+1)+'.ini'+'\n')
        
        # --- hk template file ---
        if self.load.parametrization == -1:
            
            for i in range(self.load.n_layer):
                if grid == 1:
                    fid2.write(self.load.model_name_calib+str(i+1)+'.tpl'+'\t'+self.load.model_name_calib+str(i+1)+'.pts'+'\n')
                elif grid == 2:
                    print self.load.n_layer
                    for name in pp_file_name:
                        fid2.write(name+str(i+1)+'.tpl'+'\t'+name+str(i+1)+'.pts'+'\n')
            if self.load.inv_type == 'joint':
                for k in pp_y: 
                    for j in pp_x:
                        fid2.write(self.load.mod_file+str(k)+'_'+str(j)+'.tpl'+'  '+self.load.mod_file+str(k)+'_'+str(j)+'.mod'+'\n')
        
        if self.load.parametrization == -2:
                fid2.write(self.load.modflow_file_name+'_calib_lpf.tpl'+'\t'+self.load.modflow_file_name+'_calib.lpf'+'\n')            
        
        if self.load.parametrization == -3:
            for i in range(self.load.n_layer):
                fid2.write(self.load.model_name_calib+str(i+1)+'.tpl'+'\t'+self.load.model_name_calib+str(i+1)+'.pts'+'\n')
        
        # --- instruction files ---        
        tmp = self.load.head_obs_points.split('.')
        fid2.write('head_obs_model_smp.ins  '+tmp[0]+'_model.smp'+'\n')        
        fid2.write('bud2smp_model.ins  bud2smp_model.smp'+'\n')
        
        # --- add prediction to calibration dataset ---
        if self.load.pred2calibdata < 0:
            tmp = self.load.head_pred_points.split('.')
            if self.load.pred2calibdata == -1 or self.load.pred2calibdata == -3:
                fid2.write('rchsum_dat.ins rchsum.dat'+'\n')
            if self.load.pred2calibdata == -1 or self.load.pred2calibdata == -2:
                fid2.write('head_pred_model_smp.ins  '+tmp[0]+'_model.smp'+'\n')
        
        if self.load.inv_type == 'joint':
            for k in pp_y: #range(numfiles):
                for j in pp_x: #range(num_obs):
                    if self.load.geophys_software in 'Aarhusinv64':
                        fid2.write('fwr_'+str(k)+'_'+str(j)+'.ins'+'  '+self.load.mod_file+str(k)+'_'+str(j)+'00001.fwr'+'\n')
                    if self.load.geophys_software in 'em1dinv64':
                        fid2.write('fwr_'+str(k)+'_'+str(j)+'.ins'+'  '+self.load.mod_file+str(k)+'_'+str(j)+'001.fwr'+'\n')
                        
        if self.load.inv_type != 'hydrologic':
            for layer in range(self.load.n_layer):
                fid2.write('hk_diff_L'+str(layer+1)+'.ins     hk_diff_L'+str(layer+1)+'.dat'+'\n')
        
        if self.load.inv_type == 'hydrologic' and self.load.regul_hybrid== True:
            fid2.write('* prior information'+'\n')
            
        # --- regularisation ---
        if self.load.parametrization != -2 or self.load.regul_hybrid == True:
            if PESTMODE == 'regul':
                
                fid2.write('* regularisation'+'\n')    
                fid2.write(str(PHILIM)+'   '+str(PHILIM+PHILIM*0.1)+'   '+str(FRACPHIM)+'\n')
                fid2.write('1.0   1.0e-10    1.0e10 '+'\n')
                fid2.write('1.3   1.0e-2     '+str(IREGADJ)+'\n')
        
        fid2.close()
        print 'done'
        
#==============================================================================
    def obs2folder(self,model_num,pp_x=[],pp_y=[]):
#==============================================================================
        """
        Copy hydrological and geophysical measurement observation to folder with results.
        
        Parameters
        ----------
        model_numb :        int
            geological reference structure number (realization number)
        pp_x :              List(int)  
            python list - node positions of pilot points in column direction
        pp_y :              List(int)
            python list - node positions of pilot points in row direction
        """
        import shutil 
        
        # --- hydrological observation ---
            # -- heads ---
        tmp = self.load.head_obs_points.split('.')
        src = tmp[0]+'.obs'
        dst = self.load.path_results+'/'+src
        shutil.copyfile(src,dst+'.'+str(model_num)) 
        del dst
        
            # -- river discharge ---
        tmp = self.load.head_obs_points.split('.')
        src = 'bud2smp_river.obs'
        dst = self.load.path_results+'/'+src
        shutil.copyfile(src,dst+'.'+str(model_num))
        del dst
        
        # --- geophysical observation ---
        if pp_x != [] and pp_y != []:
            for i in pp_y:
                for j in pp_x:
                    src = self.load.tem_name+str(i)+'_'+str(j)+'.tem'
                    dst = self.load.path_results+'/'+src
                    shutil.copyfile(src,dst+'.'+str(model_num))
                    del dst
                    
    
#==============================================================================
    def pilot_point_file(self,pp_x,pp_y,x_valley,y_valley,grid=1):
#==============================================================================
        """       
        Writes pilot point file 
        
        PEST pilot point file
         
        pp_x :              List(int)  
            python list - node positions of pilot points in column direction
        pp_y :              List(int)
            python list - node positions of pilot points in row direction
        x_valley :          List(int)
            python list - node positions of pilot points in column direction
            inside the buried valley
        y_valley :          List(int)
             python list - node positions of pilot points in row direction
             inside the buried valley
        grid :              int
            grid number used when running model, see ini file for discretization
        """ 
        
        print 'pilot_point_file(pts)      ',
        x_pos           = pp_x
        y_pos           = pp_y
        
        n_pp = len(pp_x)*len(pp_y)
        # --- name of pilot point files ---
        if grid == 1:
            pp_file_name = self.load.pp_file_name
        
        if grid == 2:
            pp_file_name = ['a_L','b_L']
            par_ini      = [self.load.beta_1*np.ones(n_pp),self.load.beta_2*np.ones(n_pp)]
#        test = par_ini[1]
#        print test[0],test
#        return

        #  --- load inter array ---        
        if self.load.parametrization == -1:
            # --- load zonation from integer array files ---
            zone_model      = self.load.int2mat(self.load.integer_array_name,self.load) 
        if self.load.parametrization == -2:            
            # --- load zonation from integer array files ---
            zone_model      = self.load.int2mat(self.load.integer_array_name,self.load)
        if self.load.parametrization == -3:            
            # --- load zonation from integer array files ---
            zone_model      = self.load.int2mat(self.load.integer_array_name,self.load)
   
        # --- initial parameter values ---    
        if self.load.start_model == -1:
            k = self.load.model(self.load.model_name,self.load.dtype_model,self.load.nx,self.load.ny,self.load.n_layer)
        elif self.load.start_model == -2:
            res = self.load.mod2grid(self.load.mod_file,pp_x,pp_y,self.load)
            k = np.zeros((self.load.ny,self.load.nx,self.load.n_layer),dtype=float)
            for layer in range(self.load.n_layer):
                for i in range(self.load.ny):
                    for j in range(self.load.nx):
                        k[i,j,layer] = self.load.beta_1*res[i,j,layer]**self.load.beta_2  
        
        elif self.load.start_model == -3:
            # --- allocate array ---
            k = np.zeros((self.load.ny,self.load.nx,self.load.n_layer),dtype=float)

            # --- load zones parameter into array ---
            hz = self.load.hydro_zones(self.load.zone_par_file,self.load.n_zone_hk,4)
            hk = map(float,hz[0:3,-2])
            zone = map(int,hz[0:3,-1])
        
            for layer in range(self.load.n_layer):
                for i in range(self.load.ny):
                    for j in range(self.load.nx):
                        for h in range(len(hk)):
                            if zone_model[i,j,layer] == zone[h]:
                                k[i,j,layer] = hk[h]
#                                

            
        elif self.load.start_model > 0:
            k = self.load.start_model*np.ones((self.load.ny,self.load.nx,self.load.n_layer),dtype=float)
        else:
            "WARNING, pilot point file- starting model value"
            return
        
        # ***********************************************
        # --- recalc positions form notes into meters ---
        # ***********************************************
        x_cap = []
        y_cap = []        
        for j in pp_y:
            y = float(j)*self.load.dy+self.load.dy/2.0
            y_cap = y_cap + [y]
            
        y_cap = np.flipud(y_cap)                # flip y-coor to fit pilot point format
            
        for i in pp_x:
            x = float(i)*self.load.dx+self.load.dx/2.0
            x_cap = x_cap + [x]
        
        if self.load.n_layer >= self.load.n_cap:
            x_val = []
            y_val = []
            for layer in range(self.load.n_layer-self.load.n_cap):
                x_v = []
                count = 0
                for i in x_valley[layer]:
                    x       = float(i)*self.load.dx+self.load.dx/2.0
                    x_v     = x_v + [x]
                    count += 1
                x_val.append(x_v)
                
            for layer in range(self.load.n_layer-self.load.n_cap):
                y_v = []
                count = 0
                for i in y_valley[layer]:
                    y       = float(i)*self.load.dy+self.load.dy/2.0
                    y_v     = y_v + [y]
                    count += 1
                y_val.append(y_v)
                
                
        # *******************************      
        # --- write pilot point file ---
        # *******************************                
        for layer in range(self.load.n_layer):
            if layer < self.load.n_cap:
                
                hk   = []
                zone = []
                for j in y_pos:
                    for i in x_pos:
                        # --- pick hydraulic conductivity value a pilot poit positions
                        hk = hk + [k[j,i,layer]]
                        
                        # --- pick zone number at pilot point positions --- 
                        zone = zone + [zone_model[j,i,layer]]

                if grid == 1:                          
                    self.write.pts_file(pp_file_name+str(layer+1),self.load.pp_name,x_cap,y_cap,zone,hk)
                elif grid ==2:
                    for h in range(len(pp_file_name)):                   
                        self.write.pts_file(pp_file_name[h]+str(layer+1),self.load.pp_name,x_cap,y_cap,zone,par_ini[h])
                        
            elif layer >= self.load.n_cap:
                hk   = []
                zone = []
                for j in y_valley[layer-self.load.n_cap]:
                    for i in x_valley[layer-self.load.n_cap]:
                        # --- pick hydraulic conductivity value a pilot poit positions
                        hk = hk + [k[j,i,layer]]
                        
                        # --- pick zone number at pilot point positions --- 
                        zone = zone + [zone_model[j,i,layer]]
                
                if grid == 1:                          
                    self.write.pts_file(pp_file_name+str(layer+1),self.load.pp_name,x_cap,y_cap,zone,hk)
                elif grid ==2:
                    for h in range(len(pp_file_name)):                   
                        self.write.pts_file(pp_file_name[h]+str(layer+1),self.load.pp_name,x_cap,y_cap,zone,par_ini[h])

                       

        print 'done'
        
#==============================================================================
    def pilot_point_tpl(self,left=38,grid=1,pp_file_name=[]):        
#==============================================================================
        """
        
        write pilot point template file
        
        Parameters
        ----------
        left :              int
            PEST marker, read from this "left" position to "right"
        grid :              int
            grid number used when running model, see ini file for discretization
        pp_file_name :      Optional(List(str))
            Python list with sufix to fac2real-ini files
            
        """
        
        if grid == 1:
            pp_file_name = self.load.pp_file_name
        if grid == 2:
            pp_file_name = pp_file_name
        
        for layer in range(1,self.load.n_layer+1,1):        
            self.write.pts2tpl(layer,pp_file_name,left)
            
#==============================================================================
    def ppkreg2pst(self):
#==============================================================================        
        """
        Writes ini-files for ppkreg (PEST utility) and run ppkreg
        
        """
        
        import subprocess 
        
      
        print 'ppkreg2pst:  ',
        
        # write batch file for running ppkreg
        fid1 = open('run_ppkreg.bat','w')
        fid1.write('echo off \n')
        fid1.write('ppkreg.exe < ppkreg.in > null.dat')
        fid1.close()

        for i in range(1,self.load.n_layer+1,1):
            if i == 1:
                fid1 = open('ppkreg.in','w')        
                fid1.write(self.load.PEST_name+'.pst'+'\n')        
                fid1.write('regularisation_L'+str(i)+'.dat'+'\n')
                fid1.write(str(1)+'\n')        
                fid1.write('hk'+str(i)+'_'+'\n')
                fid1.write('g'+'\n')                
                fid1.write(str(1)+'\n')             
                fid1.write('regul_'+str(i)+'\n')
                fid1.write('ppr'+str(i)+'_'+'\n')
                fid1.write('y'+'\n')             
                fid1.write('case_reg_'+str(i)+'.pst'+'\n')      
                fid1.write(str(37)+'\n')
                fid1.write(str(1)+'\n')
                fid1.write(str(1e-10)+'\n')
                fid1.write(str(1e10)+'\n')        
                fid1.close()
                subprocess.call(['run_ppkreg.bat'])
            else:
                
                fid1 = open('ppkreg.in','w')        
                fid1.write('case_reg_'+str(i-1)+'.pst'+'\n')        
                fid1.write('regularisation_L'+str(i)+'.dat'+'\n')
                fid1.write(str(1)+'\n')        
                fid1.write('hk'+str(i)+'_'+'\n')
                fid1.write('g'+'\n')
                fid1.write(str(1)+'\n')                 
                fid1.write('regul_'+str(i)+'\n')
                fid1.write('ppr'+str(i)+'_'+'\n')
                fid1.write('y'+'\n') 
                fid1.write('case_reg_'+str(i)+'.pst'+'\n')      
                fid1.write(str(37)+'\n')
                fid1.write(str(1)+'\n')
                fid1.write(str(1e-10)+'\n')
                fid1.write(str(1e10)+'\n')        
                fid1.close()
                subprocess.call(['run_ppkreg.bat'])
                
#        # --- delete first pest controle file before writing new with tikonov reg
#        try:                
#            os.remove(self.load.PEST_name+'.pst')
#        except:
#            if IOError:
#                print 'pest controle does not excist!!!'
                
        # --- write new pest controle file with the "right" PHILIM       
        fid1 = open('case_reg_'+str(i)+'.pst','r')          
        lines = fid1.readlines()
        fid1.close()

        fid1 = open(self.load.PEST_name+'.pst','w')
        find = '* regularisation'       
        for i in range(len(lines)):                        
            if find in lines[i]:
                tmp         = lines[i+1].split()
                if self.load.regul_hybrid == True and self.load.inv_type == 'sequential':
                    tmp[0]      = str(3.0)
                    tmp[1]      = str(3.3)
                else:
                    tmp[0]      = str(2.0)
                    tmp[1]      = str(2.2)
                tmp.extend('\n')
                text = " ".join(str(x) for x in tmp)
                lines[i+1] = text
            fid1.write(lines[i])
        fid1.close()
        
        print 'done'            
#==============================================================================
    def ppkreg2pst_ppzone(self):
#==============================================================================        
        """
        Writes ini-files for ppkreg (PEST utility) and run ppkreg
        Combining pilot points and zones
        
        """
        
        import os
        import subprocess 
        
      
        print 'ppkreg2pst:  ',
        
        # write batch file for running ppkreg
        fid1 = open('run_ppkreg.bat','w')
        fid1.write('echo off \n')
        fid1.write('ppkreg.exe < ppkreg.in > null.dat')
        fid1.close()
        
        if self.load.parametrization == -1:
            n_zone = 1
        if self.load.parametrization == -3:
            n_zone = len(self.load.pp2zone.split(','))+len(self.load.pp2zone_tied.split(','))
        
        n_zone = 1  
        
        for i in range(1,self.load.n_layer+1,1):
            count1 = 1
            if i == 1:
                fid1 = open('ppkreg.in','w')        
                fid1.write(self.load.PEST_name+'.pst'+'\n')        
                fid1.write('regularisation_L'+str(i)+'.dat'+'\n')
                fid1.write(str(n_zone)+'\n')
                for j in range(1,n_zone+1,1):
                    fid1.write('hk'+str(i)+'_'+'\n')
                    fid1.write('g'+'\n')                
                    fid1.write(str(1)+'\n')             
                    fid1.write('regul'+str(i)+'_'+str(j)+'\n')
                    fid1.write('r'+str(i)+'_'+str(j)+'\n')
                    fid1.write('y'+'\n')
                    fid1.write('y'+'\n')
#                    fid1.write('y'+'\n')
          
                fid1.write('case_reg_'+str(i)+'.pst'+'\n')      
                fid1.write(str(37)+'\n')
                fid1.write(str(1)+'\n')
                fid1.write(str(1e-10)+'\n')
                fid1.write(str(1e10)+'\n')        
                fid1.close()
                subprocess.call(['run_ppkreg.bat'])
                
            else:
                fid1 = open('ppkreg.in','w')        
                fid1.write('case_reg_'+str(i-1)+'.pst'+'\n')        
                fid1.write('regularisation_L'+str(i)+'.dat'+'\n')
                fid1.write(str(n_zone)+'\n')
                for j in range(1,n_zone+1,1):
                    fid1.write('hk'+str(i)+'_'+'\n')
                    fid1.write('g'+'\n')                
                    fid1.write(str(1)+'\n')             
                    fid1.write('regul'+str(i)+'_'+str(j)+'\n')
                    fid1.write('r'+str(i)+'_'+str(j)+'\n')
                    fid1.write('y'+'\n') 
                    fid1.write('y'+'\n')
#                    fid1.write('y'+'\n')

                fid1.write('case_reg_'+str(i)+'.pst'+'\n')      
                fid1.write(str(37)+'\n')
                fid1.write(str(1)+'\n')
                fid1.write(str(1e-10)+'\n')
                fid1.write(str(1e10)+'\n')        
                fid1.close()
                subprocess.call(['run_ppkreg.bat'])
                count1 += 1
        # --- delete first pest controle file before writing new with tikonov reg
        try:                
            os.remove(self.load.PEST_name+'.pst')
        except:
            if IOError:
                print 'pest controle does not excist!!!'
                
        # write new pest controle file with the "right" PHILIM       
        fid1 = open('case_reg_'+str(i)+'.pst','r')          
        lines = fid1.readlines()
        fid1.close()

        fid1 = open(self.load.PEST_name+'.pst','w')
        find = '* regularisation'       
        for i in range(len(lines)):                        
            if find in lines[i]:
                tmp         = lines[i+1].split()
                tmp[0]      = str(2.0)
                tmp[1]      = str(2.2)
                tmp.extend('\n')
                text = " ".join(str(x) for x in tmp)
                lines[i+1] = text
            fid1.write(lines[i])
        fid1.close()
        
        print 'done'


#==============================================================================
    def select_obs_point(self,node2egde=5,grid=1):         
#==============================================================================
        """
        Utility to select regular distrutibutated pilot point positions and
        geophysical observations points
    

        Parameters
        ----------
        node2egde :         int
            buffer distance to model-boundaries
        
        grid :              int
            select pilot point grid discretization from model_dis.ini
        
        Returns
        -------
        pp_x :              list (int)
            node position with pilot points in x-direction
        
        pp_y :              list (int)
            node position with pilot points in y-direction
            
        """
        
        # --- select_obs_point ---
        pp_x           = []
        pp_y           = []
        
        print 'select_obs_point:          ',
        # --- load first set of pilot point from model_dis.ini ---
        if grid == 1:
            n_x = self.load.n_pp_x
            n_y = self.load.n_pp_y
        # --- load second set of pilot point from model_dis.ini ---
        elif grid == 2:
            n_x = self.load.n_pp_x_2
            n_y = self.load.n_pp_y_2
        
        # --------------------------------------------------        
        # --- position of measurement points --- 
        # --------------------------------------------------
        
        
        dis_x = int(self.load.nx/n_x/2.0)
        pp_x = np.round(np.linspace(dis_x+node2egde,self.load.nx-dis_x-node2egde,n_x))     
        pp_x = map(int,pp_x)
        
        dis_y = int(self.load.ny/n_y/2.0)
        pp_y = np.round(np.linspace(dis_y+node2egde,self.load.ny-dis_y-node2egde,n_y))
        pp_y = map(int,pp_y)
        
        print 'done'
        return pp_x, pp_y 


#==============================================================================
    def realization2folder(self,model_numb):
#==============================================================================        
        """
        Copy zone realization file (from blocksis) to work directory
        
        Parameters
        ----------
        model_numb :        int
            realization number to copy from "blocksis_path"/realization folder to work folder
            
        
        
        """
        import shutil
        
        cwd = os.getcwd()
               
        # --- copy zone file  from seperate BLOCKSIS folder---            
        source  = self.load.blocksis_path+'/'+self.load.blocksis_name+'.'+str(model_numb)
        
        # --- copy into modflow inter array file name specified in hydro_setup.ini --- 
        dest    = cwd+'/'+self.load.modflow_zone_file
        
        try:
            shutil.copyfile(source,dest)
        except:
            print 'Error in realization2folder....'
            print cwd

        

#==============================================================================
    def recharge_tpl(self,inifile='recharge_est',left=0,right=39,split='REM'):
#==============================================================================
        """
        Writes PEST template file
        
        inifile :           str
            instruction-file for calc recharge based on hk of upper most layer
        left :              int
            PEST marker, read from this "left" position to "right"
        right :             int
            PEST marker, read to this "left" position from "left"
        split :             str
            string delimiter
        
        """

        print 'recharge_tpl:                ',
        
        

        tmp = self.load.rch_est_calc.split('.')
        fid1 = open(tmp[0]+'.tpl','w')
        fid1.write('ptf $ \n')   
        
        # --- pilot point parametrization ---
        if self.load.parametrization == -1 or self.load.parametrization == -3:
            # --- load ini-file for recharge_est.py ---
            fid2 = open(tmp[0]+'.ini','r')
            text = fid2.readlines()
            fid2.close()
#            print text
            # --- write tpl file ---              
            if self.load.rch_est == -1:
                tex = text[0].split(split)   
                tmp ='%'+str(left)+'.1s %0.11s %'+str(right-left-3-3)+'.1s %5.'+str(len(tex[1])+5)+'s'
                fid1.write(tmp % ('$','recharge_a','$ ',split+' '+tex[1]))
                tex = text[1].split(split)        
                fid1.write(tmp % ('$','recharge_b','$ ',split+' '+tex[1]))
                tmp ='%'+str(left)+'.1s %0.11s %'+str(right-left-3-3)+'.1s %5.'+str(len(tex[1])+5)+'s'
                for i in range(2,6,1):
                    fid1.write(text[i])
                
                
            elif self.load.rch_est == -2:
                tex = ' Recharge zone estimates'
                tmp ='%'+str(left)+'.1s %0.11s %'+str(right-left-3-3)+'.1s %5.'+str(len(tex)+5)+'s'
                for i in range(1,self.load.n_zone_rch+1,1):
#                    print       i
                    fid1.write(tmp % ('$','recharge_'+str(i),'$ ',split+' '+tex+'\n'))
           
                for i in range(self.load.n_zone_rch,len(text),1):
                    fid1.write(text[i])
#                print i,    text[i]
        
        # --- zone parametrization ---        
        if self.load.parametrization == -2:
            fid2 = open(self.load.modflow_file_name+'_calib'+'.rch','r')
            text = fid2.readlines()
            fid2.close()
            for i in range(len(text)):
                if i < 3:                
                    fid1.write(text[i])
                if i >= 3 and i < 3*self.load.n_zone_rch:
                    if i%2 != 0:                    
                        tmp1 = text[i].split()
                        tmp2 = ('%1.'+str(len(tmp1[0]))+'s'+
                                '%'+str(len(tmp1[1])+1)+'.'+str(len(tmp1[1]))+'s'+
                                '%2.1s'+
                                '%'+str(len(tmp1[0])+1)+'.'+str(len(tmp1[0]))+'s'+
                                '%'+str(right-left)+'.1s'+
                                '%'+str(len(tmp1[len(tmp1)-1])+1)+'.'+str(len(tmp1[len(tmp1)-1]))+'s'+
                                '\n')
                        fid1.write(tmp2 % (tmp1[0],tmp1[1],'$',tmp1[0],'$',tmp1[len(tmp1)-1]))
                    else:
                        fid1.write(text[i])
                if i >= 3*self.load.n_zone_rch:
                    fid1.write(text[i])
                    
        fid1.close()
        
        print 'done'

#==============================================================================
    def recharge_est(self):
#==============================================================================
        """
        
        """
        
        fid1 = open('recharge_est_test.py','w')
        
        fid1.write('# -*- coding: utf-8 -*- \n')
        fid1.write('#\n')
        fid1.write('# recharge_est.py \n')
        fid1.write('# ---------------------------------------------------\n')
        fid1.write('#   24-01-2014\n')
        fid1.write('#   Nikolaj Kruse Christensen\n')
        fid1.write('# ---------------------------------------------------\n')
        fid1.write('#   Description:\n')
        fid1.write('#\n')
        fid1.write('#   INPUT:\n')
        fid1.write('# hk_L1.ref\n')
        fid1.write('#   OUTPUT:\n')
        fid1.write('# recharge_est.ref\n')
        fid1.write('# ---------------------------------------------------\n')
        fid1.write('import numpy as np\n')
        fid1.write('import os\n')
        fid1.write('\n')
        fid1.write('# --- import info from ini file ---\n')
        fid1.write("fid1 = open('recharge_est.ini','r')\n")
        fid1.write('text = fid1.readlines()\n')
        fid1.write('tal = np.shape(text)\n')
        fid1.write('tal = tal[0]\n')
        fid1.write('fid1.close()\n')
        fid1.write('\n')
        fid1.write('ini =[]\n')
        fid1.write('for i in range(tal):\n')
        fid1.write('    tex = text[i].split()\n')
        fid1.write('    ini.append(tex[0])\n')
        fid1.write('    \n')
        fid1.write('a           = float(ini[0])\n')
        fid1.write('b           = float(ini[1])\n')
        fid1.write('input_file  = ini[2]\n')
        fid1.write('output_file = ini[3]\n')
        fid1.write('nrow        = int(ini[4])\n')
        fid1.write('ncol        = int(ini[5])\n')
        fid1.write('\n')
        fid1.write('# --- load k-field ---\n')
        fid1.write("fid1 = open(input_file,'r')\n")
        fid1.write('hk = np.array(fid1.read().split(),dtype=float)\n')
        fid1.write('\n')
        fid1.write('hk_log = np.log10(hk)\n')
        fid1.write('\n')
        fid1.write('# --- calc recharge with parameter a and b in loglog space ---\n')
        fid1.write('recharge_log = np.log10(a)+np.log10(b)*hk_log\n')
        fid1.write('recharge_est = 10**recharge_log\n')
        fid1.write('\n')
        fid1.write('\n')
        fid1.write('# --- print rechareg file --- \n')
        fid1.write("fid1 = open(output_file,w')\n")
        fid1.write('count = 0\n')
        fid1.write('for i in range(nrow):\n')
        fid1.write('    for j in range(ncol):\n')
        fid1.write("        fid1.write(str(recharge_est[count])+'\t')\n")
        fid1.write('        count += 1\n')
        fid1.write("    fid1.write('\n')\n")
        fid1.write('fid1.close()\n')
        

    

#==============================================================================
    def replace_NOPTMAX(self,niter):
#==============================================================================
        """
        Replace number of iteration in PEST control file
        Find line containing NOPTMAX and Replace NOPTMAX with new value. 
        
        Parameters
        ----------
        niter :             int
            number of iteration
        """
        
        # --- make copy of PEST control file ---
        dos_tex = 'parrep.exe case.par case.pst case_opt_tmp.pst'
        subprocess.call(dos_tex)
        
        # ---------------------------------------------------
        #                   --- INPUT ---
        # ---------------------------------------------------

        tmp_file      = 'case_opt_tmp.pst'
        output_file   = 'case_opt.pst'
        
        # ---------------------------------------------------
        #               --- IN/OUTPUT FILES
        # ---------------------------------------------------
        fid1 = open(tmp_file,'r')
        text = fid1.readlines()
        fid1.close()
        
        
        # ---------------------------------------------------
        #       --- Find Line and replace NOPTMAX --- 
        # ---------------------------------------------------
        fid1 = open(output_file,'w')
        for i in range(len(text)):
            if i == 8:
                temp = text[i].split()
                fid1.write(str(niter))
                for j in range(len(temp)-1):
                    fid1.write(' ')
                    fid1.write(temp[j+1])
                fid1.write('\n')                                  
            else:            
                fid1.write(text[i])
        fid1.close()
        
        # --- delete temp PEST controle file with optimal par values ---
        os.remove(tmp_file)
        
        # --- run pest with estimated parameter values ---
        dos_tex = 'pest.exe case_opt.pst' 
        subprocess.call(dos_tex)

        print('done')

#==============================================================================
    def run_fac2real(self,input_file):
#==============================================================================
        """
        run PEST utility to interpolate the pilot point values intop the groundwater model grid 
        """
              
        print 'run_fac2real:              ',
        for i in range(self.load.n_layer):
            myinput = open(input_file+str(i+1)+'.in','r')
            tmp     = open('null_Fac2real_L'+str(i+1)+'.dat','w')
            subprocess.call('fac2real.exe',stdin=myinput,stdout=tmp)
            myinput.close()
            tmp.close()       
            
        print 'done'
            
#==============================================================================
    def run_ppk2fac(self,input_file):
#==============================================================================
        """
        run PEST utility to calculate regularisation factors 
        """          
        
        print 'run_ppk2fac:             ',
        for i in range(self.load.n_layer):
            myinput = open(input_file+str(i+1)+'.in','r')
            tmp     = open('null.dat','w')
            subprocess.call('ppk2fac.exe',stdin=myinput,stdout=tmp)
            myinput.close()
            tmp.close()       
            
        print 'done'
            
