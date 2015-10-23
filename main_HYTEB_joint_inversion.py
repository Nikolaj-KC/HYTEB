# -*- coding: utf-8 -*-
"""
Created on Thu Jun 05 09:37:53 2015

@author: nikolaj.christensen

Description:
Script for running joint inversion.
"""



import sys

# --- directory of HYTEB---
sys.path.append('C:/Users/nikolaj/Dropbox/PH_D/Scripts/Test_bench')

# --- import HYTEB modules ---
import hydro_module
import geophys_module
import inversion_module

hydro   = hydro_module.hydro_module(hydro_ini = 'hydro_ini_joint.ini',inversion_ini='inversion_setup_joint.ini')
gp      = geophys_module.geophys_module(hydro_ini = 'hydro_ini_joint.ini')
invs    = inversion_module.inversion_module(hydro_ini = 'hydro_ini_joint.ini',inversion_ini='inversion_setup_joint.ini')

# --- "true" reference systems for sampling geophysical and hydrological data sets ---
# --- realizations have been chosen from the clustering of the reference predictions ---     
realization = [141,202,234,359,450,464,678,727,944,991]

# --- run HYTEB analysis for multiple realizations ---
for i in realization:  
    
    # --- delete all "old input/out put files and PEST slaves ---    
    invs.delete_files(del_slave='on')
    
    # --- copy reference system from folder(see ini-files) to work directory --- 
    invs.model2folder(i)
    
    #==============================================================================
    # Defind groundwater model parametrization and boundaries
    #==============================================================================
    
    # --- create MODFLOW ibound file for defining the active part of the groundwater model ---
    hydro.ibound(x_top=1500,x_bot=500,valley_north=5000,fixed_head=280)
    
    # --- select pilot point positions for the "capping" part of groundwater model ---
    [pp_x,pp_y]                         = invs.select_obs_point()
    
    # --- select pilot point positions for the "valley" part of groundwater model ---
    [x_valley,y_valley,n_pp_valley]     = invs.obs_point_in_valley(pp_x,pp_y)
    
    #==============================================================================
    # geophysical observations
    # quasi 3D groundbased TEM 40x40 
    #==============================================================================  
    
    # ---sample quasi 3D geophysicakl models ---
    gp.tem_mod_quasi_3D(pp_x,pp_y)
    
    # --- generate TEM forward responses ---
    gp.run_aarhus_inv('modfil_',pp_x,pp_y)
    
    # --- added noise to TEM forward responses ---
    gp.fwr2tem(pp_x,pp_y)
    
    # --- proccess the TEM data as the were field data --- 
    n_tem_data = gp.process_tem(pp_x,pp_y)
    
    # --- homogeneus starting model for the geophysical inversions --- 
    gp.mod2mod('modfil_','mod_',pp_x,pp_y,homogen=40.0,n_iter = -1)

    
    #==============================================================================
    # make groundwater model files for modflow-2000
    #==============================================================================

    hydro.model_specification()
    hydro.integer_array_file()
    
    # --- make MODFLOW files both for generating observation data and in calibration mode ---
    hydro.modflow_bas(sufix='',external=0)
    hydro.modflow_bas(sufix='_calib',external=75)
    hydro.modflow_dis()
    hydro.modflow_gmg()
    
    hydro.modflow_lpf()    
    hydro.modflow_lpf(sufix='_calib',run_type = 'calib')
    hydro.modflow_mpt()
    hydro.modflow_mpn(run_type='calib',well='on' ,recharge='calib',sufix='_calib_pumping')
    hydro.modflow_mpn(run_type='calib',well='off',recharge='calib',sufix='_calib')
    hydro.modflow_wel()
    hydro.modflow_oc(head=51,drawdown=52)
    hydro.modflow_rch(run_type='true')
    hydro.modflow_rch(run_type='calib')
    
    # --- writing mf nam file and writing mf river file (.riv) ---
    hydro.modflow_nam('valley_river',run_type='true',well='on',river='off')
    hydro.modflow_river()
    
    # --- writing new mf nam file --- 
    hydro.modflow_nam('valley_true',run_type='true',well='off',river='on')
    hydro.modflow_nam('valley_pump_true',run_type='true',well='on',river='on')
    
    # --- run MODFLOW-2000 with reference system for generating observations data ---
    hydro.run_mf2k('valley_pump_true',save_initial_head='on')
    hydro.run_hydro_obs()
    
    # --- add noise to hydraulic measurements ---
    hydro.smp2obs()
    
    # --- writing new mf nam file and run mf for calibration---
    hydro.modflow_nam('valley_calib_pumping',run_type='calib',well='on',river='on',initial_head='on')
    hydro.modflow_nam('valley_calib',run_type='calib',well='off',river='on',initial_head='on')
    
    # --- modpath5_0 ---
    hydro.modflow_mpt()
    hydro.modflow_mpn(run_type='calib',well='on' ,recharge='calib',sufix='_calib_pumping')
    hydro.modflow_mpn(run_type='calib',well='off',recharge='calib',sufix='_calib')
    hydro.modflow_rsp('valley_calib_pumping.mpn','particles.ini',sufix='_pumping')
    hydro.modflow_rsp('valley_calib.mpn','particle.ini',sufix='')
    
    #==============================================================================
    #  PEST inversion setup
    #==============================================================================
    
    # --- make PEST pilot point files --- 
    invs.pilot_point_file(pp_x,pp_y,x_valley,y_valley)
    
    # --- make ini files for PEST utility for interpolating estimated K and R from ---
    # --- pilot point positions into groundwater model grid ---
    invs.in_fac2real()
    invs.in_ppk2fac('ppk2fac_L')
    invs.run_ppk2fac('ppk2fac_L')
    
    # --- make PEST instruction files --- 
    invs.ins_files_diff('hk_diff_L',pp_x,pp_y,x_valley,y_valley)
    invs.ins_files_head('head_obs_model_smp')
    invs.ins_files_river('bud2smp_model')
    invs.ins_files_tem('fwr_',pp_x,pp_y,n_tem_data)
    
    # --- make PEST template files ---
    invs.modfile_tpl(pp_x,pp_y)
    invs.pilot_point_tpl() 
    invs.diff_hk_tpl()
    invs.recharge_tpl()
    invs.diff_hk_tpl()
    
    # --- make cmd batch files for calling different model responses ---
    invs.batch_all_models('valley_calib_pumping',pp_x,pp_y)
    invs.batch_gw_model('valley_calib_pumping',pp_x,pp_y)
    invs.batch_gp_model()
    
    # --- make PEST control file ---
    invs.PEST_inversion(pp_x,pp_y,x_valley,y_valley,n_pp_valley,n_tem_data=n_tem_data)    
    
    # --- generating slave for paralization --- 
    invs.copy2slaves(pp_x,pp_y,tmp_numb=i) 
    
    # --- run joint inversion with BeoPEST ---
    invs.batch_beopest(run_inversion='on')
    
    # --- update estimated parameter file and hydraulic responses ---
    invs.replace_NOPTMAX(0)
    
    # --- run hydraulic model prediction --- 
    hydro.run_prediction(i,run_type='calib')   
    
