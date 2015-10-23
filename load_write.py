# -*- coding: utf-8 -*-
"""
Created on Mon Aug 11 13:01:49 2014

@author: Nikolaj
"""

import os
import numpy as np

# *****************************************************************************
# *****************************************************************************
class load_file:
    """
    class with loading and writing modules. 
    Load ini-files values into handle "load" 
    


    """
#==============================================================================
    def eas_borehole(self,name_list,kote_bottom,kote_top,prop_clay,utm_x=0,utm_y=0,utm_z=0):    
#==============================================================================
        """
        Load eas borehole file and returns numpy array.
        
        Paramaters
        ----------
        name_list :         List(str)
            python list[n]. List with n names of "clay-fraction" bor-files.
            One bor-file for each layer in the model.  
        kote_bottom :       float
            Bottom of layer for a givin depth interval            
        kote_top :          float
            Top of layer for a givin depth interval
        prop_clay :         List(float)
            Python list(n). List with n fraction of the sediments in a depth 
            interval.
        utm_x :             int
            definde coordinate system. Southern-east corner in depth of n_layer is basis 
        utm_y :             int
            definde coordinate system. Southern-east corner in depth of n_layer is basis         
        utm_z :             int
            definde coordinate system. Southern-east corner in depth of n_layer is basis
        
        Returns
        -------
        text :              List(float)
            Python list(m,n). 
            m is the number of boreholes spaning a depth interval.
            n = 7. utmx, utmy, utmz, minimun CF, CF, maximum CF, well-name.    
        
        Notes
        -----
        Clay-fraction are using bottom layer of the south-west corner as Origo.
        
        """
        
        data = []
        for fil in name_list:
            fid2 = open(fil,'r')
            tmp = fid2.readlines()
            fid2.close()
            line    = tmp[1].split() 
            for i in range(int(line[0]),len(tmp),1):
                data.append(tmp[i].split())
                
        STD = 0.05
        count = 7
        count = 0
        text= []
        
        while count < len(data):    
            try:
                if kote_bottom <= float(data[count][2]) <= kote_top:
                    count2 = 0
                    if kote_top > 150:
                        for i in range(len(prop_clay)):
                            if int(data[count][3+i]) == 1:
                                CT = prop_clay[count2]*(kote_top-kote_bottom) 
                            count2 += 1
                        CT_min = CT-CT*STD
                        CT_max = CT+CT*STD
                        well_dgu = data[count][6]
                    if kote_top <= 150:
                        for i in range(len(prop_clay)):
                            if int(data[count][3+i]) == 1:            
                                CT = prop_clay[count2]*(kote_top-kote_bottom)
                            count2 += 1
                        
                        CT_min = CT-CT*STD
                        CT_max = CT+CT*STD
                        well_dgu = data[count][7]
                        
                    text.append((utm_x+float(data[count][0]),utm_y+float(data[count][1]),utm_z+float(data[count][2]),CT_min,CT,CT_max,'well_'+well_dgu))
                count += 1
            except:
                pass        
                count += 1
        return text
#==============================================================================
    def emofile(self,fname,line,n_layer):
#==============================================================================
        """
        Loads the AarhusInv inversion results and return parameters,STD and DOI. 
        
        Parameters
        ----------
        fname :             str
            Emo-file name.
        line :              int
            Line number in emo file, that holdes the number of iterations.
        n_layer :           int
            Number of layers in groundwater model. The geophysical model 
            consist of "n\_layer+1" layers.
            
        Returns
        -------
        inv_result :        List(float)        
            Estimated electirical resistivities 
        par_err :           List(float)
            Err from inversion data analysis
        doi :               List(float)
            Depth of investigation (DOI)
        """

        # --- load emo file ---        
        fid1 = open(fname,'r')
        text = fid1.readlines()  
        fid1.close()

        # --- go to line and find estimated parameters ---
        n_iter      = int(text[line])
        tmp         = map(float,text[46+5+2*n_iter].split())
        
        inv_result  = []
        for i in range(2*n_layer+1):
            inv_result.append(tmp[i+1])        
        
        del tmp

        # --- adding depth to results ---
        depth = 10
        for i in range(n_layer+1,2*n_layer+1,1):
            inv_result.append(depth)
            depth  += inv_result[i]
        
        # --- go to line and find inversion analysis result ---
        tmp = map(float,text[46+5+2*n_iter+3].split())
        par_err = tmp[0:n_layer+1]
        
        # --- find DOI ---
        count = 0
        doi = []
        find = 'DOI absolute'
        while count < len(text):
            try: 
                if find in text[count]:
                    tmp = text[count+1].split()
                    doi = (float(tmp[0]))#+float(tmp[1]))/2.0 
                count += 1
            except:
                print 'Error in load emofile - DOI!'
        
        # --- return inversion result ---        
        return inv_result,par_err,doi

#==============================================================================
    def getIPAddresses(self):
#==============================================================================
        """
        Returning IP adress of machine
        
        """
        from ctypes import Structure, windll, sizeof
        from ctypes import POINTER, byref
        from ctypes import c_ulong, c_uint, c_ubyte, c_char
        MAX_ADAPTER_DESCRIPTION_LENGTH = 128
        MAX_ADAPTER_NAME_LENGTH = 256
        MAX_ADAPTER_ADDRESS_LENGTH = 8
        class IP_ADDR_STRING(Structure):
            pass
        LP_IP_ADDR_STRING = POINTER(IP_ADDR_STRING)
        IP_ADDR_STRING._fields_ = [
            ("next", LP_IP_ADDR_STRING),
            ("ipAddress", c_char * 16),
            ("ipMask", c_char * 16),
            ("context", c_ulong)]
        class IP_ADAPTER_INFO (Structure):
            pass
        LP_IP_ADAPTER_INFO = POINTER(IP_ADAPTER_INFO)
        IP_ADAPTER_INFO._fields_ = [
            ("next", LP_IP_ADAPTER_INFO),
            ("comboIndex", c_ulong),
            ("adapterName", c_char * (MAX_ADAPTER_NAME_LENGTH + 4)),
            ("description", c_char * (MAX_ADAPTER_DESCRIPTION_LENGTH + 4)),
            ("addressLength", c_uint),
            ("address", c_ubyte * MAX_ADAPTER_ADDRESS_LENGTH),
            ("index", c_ulong),
            ("type", c_uint),
            ("dhcpEnabled", c_uint),
            ("currentIpAddress", LP_IP_ADDR_STRING),
            ("ipAddressList", IP_ADDR_STRING),
            ("gatewayList", IP_ADDR_STRING),
            ("dhcpServer", IP_ADDR_STRING),
            ("haveWins", c_uint),
            ("primaryWinsServer", IP_ADDR_STRING),
            ("secondaryWinsServer", IP_ADDR_STRING),
            ("leaseObtained", c_ulong),
            ("leaseExpires", c_ulong)]
        GetAdaptersInfo = windll.iphlpapi.GetAdaptersInfo
        GetAdaptersInfo.restype = c_ulong
        GetAdaptersInfo.argtypes = [LP_IP_ADAPTER_INFO, POINTER(c_ulong)]
        adapterList = (IP_ADAPTER_INFO * 10)()
        buflen = c_ulong(sizeof(adapterList))
        rc = GetAdaptersInfo(byref(adapterList[0]), byref(buflen))
        if rc == 0:
            for a in adapterList:
                adNode = a.ipAddressList
                while True:
                    ipAddr = adNode.ipAddress
                    if ipAddr:
                        yield ipAddr
                    adNode = adNode.next
                    if not adNode:
                        break
        
#==============================================================================
    def ini_geophys(self,geophys_file):
#==============================================================================
        """       
        Load ini file for geophysical settings 
        
        """
        
        fid1 = open(geophys_file,'r')
        text = fid1.readlines()
        tal = len(text)        
        fid1.close()
        ini = []
        for i in range(tal): 
            tex = text[i].split()
            ini.append(tex[0])  
        try:            
            self.geophys_software       = str(ini[0])
            self.beta_1                 = float(ini[1])
            self.beta_2                 = float(ini[2])
            self.noise_petro            = float(ini[3])
            self.res_clay               = float(ini[4])
            self.mod_header             = str(ini[5])
            self.n_mod_header           = int(ini[6])
            self.mod_file_temp          = str(ini[7])
            self.mod_file               = str(ini[8])
            self.niter_mod_temp_gp      = int(ini[9])
            self.niter_mod_gp           = int(ini[10])
            self.tem_name_raw           = str(ini[11])
            self.tem_name                   = str(ini[12])
            self.data_template              = str(ini[13])
            self.n_data_header              = int(ini[14])
            self.current                    = float(ini[15])
            self.noise_gauss                = float(ini[16])
            self.noise_em                   = float(ini[17])
            self.tem_angel                  = float(ini[18])
            self.gp_noise_model             = int(ini[19])
            self.gp_noise_model_file        = str(ini[20])
            self.gp_relative_std_as_weight  = float(ini[21])
        except:
              if IndexError:
                  print 'INDEX ERROR in geo_phys.ini'
                  pass
        
        
        
#==============================================================================    
    def ini_hydro(self,hydro_file):
#==============================================================================    
        """       

        Load ini file for hydraulic settings 
       
        Parameters
        ----------
        hydro_file :        str
            string with full name of ini file
        
        Returns
        -------
        pred2calibdata :   int
            add predictions to calibration dataset.
            -1 -> heads and particel tracing.
            -2 -> head only.
            -3 -> particle tracing only.
            
        """
        fid1 = open(hydro_file,'r')
        text = fid1.readlines()
        tal = len(text)        
        fid1.close()
        
        ini = []
        for i in range(tal): 
            tex = text[i].split()
            ini.append(tex[0])  
        try:            
            self.model_name         = str(ini[0])
            self.dtype_model        = str(ini[1])
            self.ibound             = str(ini[2])
            self.dtype_ibound       = str(ini[3])
            self.por_name           = str(ini[4])
            self.dtype_por          = str(ini[5])
            self.recharge           = str(ini[6])
            self.dtype_recharge     = str(ini[7])
            self.model_name_calib   = str(ini[8])
            self.modflow_file_name  = str(ini[9])
            self.n_well             = int(ini[10])
            self.std_well           = float(ini[11])
            self.n_river            = int(ini[12])
            self.noise_river        = float(ini[13])
            self.model_spc          = str(ini[14])
            self.head_obs_points    = str(ini[15])
            self.head_pred_points   = str(ini[16])
            self.integer_array_name = str(ini[17])
            self.settings_file      = str(ini[18])
            self.path_results       = str(ini[19])
            self.n_zone_hk          = int(ini[20])
            self.n_zone_rch         = int(ini[21])
            self.zone_par_file      = str(ini[22])
            self.modflow_zone_file  = str(ini[23])
            self.hydro_noise_model            = int(ini[24])
            self.hydro_noise_model_file       = str(ini[25])
            self.pred2calibdata     = int(ini[26])
            self.n_pred_well             = int(ini[27])
        except:
              if IndexError:
                  
                  pass
        
#==============================================================================
    def ini_inversion(self,inversion_file):
#==============================================================================
        """
        Load inversion setup settings
            
        Parameters
        ----------
        inversion_file :    str
            Full name of ini file
            
        Returns
        -------            
        rch_est :           int
            estimation of recharge. 
             0 -> no, 
            -1 -> calc rch from upper most hk layer with liner log-log relationship.
            -2 -> use zones and est rch for each zone
        
        petro_est :         int
            estimation of petro-physical realtionship. a and b in linear log-log
            realitionship; b*log(res)+log(a)
             0 -> no
            -1 -> est only one set of a and b
            -2 -> est a and b for each layer
            -3 -> est a and b in pilot point format and krig values into grid. 
            
        """
        
        import ast        
        
        fid1 = open(inversion_file,'r')
        text = fid1.readlines()
        tal = len(text)        
        fid1.close()
        

        ini = []
        for i in range(tal): 
            tex = text[i].split()
            ini.append(tex[0])  
 
        try: 
            self.pp_file_name       = str(ini[0])
            self.pp_name            = str(ini[1])        
            self.inv_type           = str(ini[2])
            self.start_model        = float(ini[3])
            self.PEST_niter         = int(ini[4]) 
            self.PEST_derforgive    = str(ini[5])
            self.PEST_lamforgive    = str(ini[6]) 
            self.rch_est            = int(ini[7])        
            self.petro_est          = int(ini[8])
            self.struct_noise       = float(ini[9])        
            self.parametrization    = int(ini[10])
            self.tied_res           = str(ini[11])
            self.model_path         = str(ini[12])
            self.n_slaves           = int(ini[13])
            self.PEST_name          = str(ini[14])
            self.IP_addr            = str(ini[15])
            self.IP_port            = int(ini[16])
            self.reg_calc           = str(ini[17])
            self.rch_est_calc       = str(ini[18])
            self.pp2zone            = str(ini[19])
            self.pp2zone_tied       = str(ini[20])
            self.blocksis_name      = str(ini[21])
            self.blocksis_n_real    = int(ini[22])
            self.blocksis_path      = str(ini[23])
            tmp                     = str(ini[24])
            self.regulasobs         = ast.literal_eval(tmp)
            del tmp 
            tmp                     = str(ini[25])
            self.regul_hybrid       = ast.literal_eval(tmp)
            self.petro_est_calc     = str(ini[26])       
            pass
        except:
              if IndexError:
                  print "IndexError - Inversion ini"
                  pass
              else:
                  print 'Error in "ini_inversion.ini"'

              
#==============================================================================    
    def ini_model_dis(self,model_dis_file):
#==============================================================================    
        """
        Load model discretization
        
        Parameters
        ----------
        model_dis_file :    str
            string with full name of ini file
            
        Returns
        -------
        dx :                float
            Cell size in x-direction.            
        dy :                float
            Cell size in y-direction.            
        dz :                float
            Cell size in z-direction.
        n_layer :           int
            Number of layers of model
        n_cap :             int
            The model is genereted with TPROGS in goes. First run generated the 
            capping part of the model consisting of n layers.
            Below the capping part is a n_layer-n_cap tick layer of heayv clay 
            with a burid valley erroreded into it. 
        x :                 float
            Size of model domain in meters in east-west direction
        y :                 float
            Size of model domain in meters in north-south direction
        nx :                int
            Number of nodes in x-direction.
        ny :                int
            Number of nodes in y-direction.
        n_pp_x :            int
            Number of pilot points in x-direction. Referred to as "grid1"
        n_pp_y :            int
            Nuber of pilot points in y-direction. Referred to as "grid1"
        n_pp_x_2 :          int
            Number of pilot points in x-direction. Referred to as "grid2"
        n_pp_y_2:            int
            Nuber of pilot points in y-direction. Referred to as "grid2"
        """
        
        
        fid1 = open(model_dis_file,'r')
        text = fid1.readlines()
        tal = len(text)        
        fid1.close()
        
        ini = []
        for i in range(tal): 
            tex = text[i].split()
            ini.append(tex[0])  
                    
        self.dx                 = int(ini[0])
        self.dy                 = int(ini[1])
        self.dz                 = int(ini[2])
        self.n_layer            = int(ini[3])
        self.n_cap              = int(ini[4])
        self.x                  = int(ini[5])      
        self.y                  = int(ini[6])
        self.nx                 = int(ini[7])
        self.ny                 = int(ini[8])
        self.n_pp_x             = int(ini[9])
        self.n_pp_y             = int(ini[10])
        try:
            self.n_pp_x_2             = int(ini[11])
            self.n_pp_y_2             = int(ini[12])
        except: 
            print 'pp estimated parameters not loaded => models_dis.ini!!!' 
                       
            pass
            
        
#==============================================================================
    def int2mat(self,fname,load):        
#=============================================================================
        """
        Load integer arrays into matrice
        
        Parameters
        ----------
        fname :             str
            string with file name of intger array
        load :              obj
            Object holding all initilized values read from ini files
            
        Returns
        -------
        int_mat :           numpy.array(int)
            numpy array with zones in 3D (ny,nx,nz)
            
            
        """
        int_mat = np.zeros((load.ny,load.nx,load.n_layer),dtype=int)
        for i in range(load.n_layer):
            tmp = fname.split('.')
            fid1 = open(tmp[0]+'_'+str(i+1)+'.'+tmp[-1],'r')        
            data = np.array(fid1.read().split()).reshape(load.ny,load.nx)
            fid1.close()

            int_mat[:,:,i] = data
                        
            
        return int_mat
        
#==============================================================================
    def hydro_zones(self,fname,n_par,n_col):
#==============================================================================
        """
        Read ini-file with zone values for hydraulic conductivity and recharge
        
        Parameters
        ----------
        fname :             str
            ini file name
        n_par :             int
            number of hk parameters. The number og recharge parameters must
            match the number of hk parameter values  
        n_col :             int
            number of columns in ini-file
        
        Returns
        -------
        hydro_zones :       numpy.array(str)
            Numpy array with zone parameter values (3*n_par,n_col).
            
        """        
        
        try:            
            fid2 = open('hydro_zones.ini','r')
            hydro_zones = np.array(fid2.read().split()).reshape(3*n_par,n_col)
            fid2.close()              
                           
        except:
            print 'WARING HYDRO_ZONES.ini does not excist'
            
        return hydro_zones   
#==============================================================================
    def model(self,model_name,dtype_model,nx,ny,n_layer):        
#==============================================================================
        """
        Load 3D model into numpy array
        
        Parameters
        ----------
        model_name :        str
            Name of real-array file.
        dtype :             str
            data type of real-array files
        nx :                int
            Number of nodes in x-direction.
        ny :                int
            Number of nodes in y-direction.
        n_layer :           int
            Number of layers of model.
            
        Returns
        -------
        model :    numpy.array((ny,nx,n_layer),float)
            3D array of model
        """
        model = np.zeros((ny,nx,n_layer),dtype=float)
        
        for i in range(n_layer):
            fid1 = open(model_name+str(i+1)+'.'+dtype_model,'r')
            model[:,:,i] = np.array(fid1.read().split(),dtype=float).reshape(ny,nx)
            fid1.close()
            
        return model
        
#==============================================================================    
    def modfile(self,fname,header_lines):
#==============================================================================        
        """
        Read AarhusInv mod-files and retruns two List.
        - Header-part.
        - Model-part.
        
        Parameters
        ----------
        fname :             str
            Input name of mod-file.
        header_lines:       int
            Number of header lines in mod-file
        
        Returns
        -------
        h_text :            List(str)
             python list containing the header of the mod-file.
        mod_par :           List(float)
             Mod-file  parameters and constrains.

        """        
        h_text   = []        
        mod_par  = []        

        fid1    = open(fname,'r')               
        # --- read header file ---
        for i in range(header_lines):
            text = fid1.readline().split('\n')
            h_text.extend(text[0:-1])

        # --- read mod-data part ---
        n_layer = int(h_text[-1][0])        
        n_layer = int(h_text[-1]) 
        for i in range(3*n_layer-2):    
            mod_par.append(map(float,fid1.readline().split()))
                
        fid1.close()
        
        
                
        return h_text,mod_par   

#==============================================================================
    def mod2grid(self,fname,pp_x,pp_y,load):
#==============================================================================
        """
        Load AarhusInv mod-files (1D model) and store resistivities in grid and
        return grid
        
        Parameters
        ----------
        fname :             str
            name of AarhusInv mod-file.
        pp_x :              List(int)  
            python list - node positions of pilot points in column direction.
        pp_y :              List(int)
            python list - node positions of pilot points in row direction.
        load :             obj
            Object holding all initilized values read from ini files

        Returns
        -------
        mod2grid :          numpy.array((ny,nx,n_layer),float)
            3D array holding geophysical parameters at each measurement position
            (pp_x,pp_y).
        
        
        """        
        
        
        print 'mod2grid:            ',
        
        mod2grid = np.zeros((load.ny,load.nx,load.n_layer),dtype=float)
        
        for i in pp_y:
            for j in pp_x:
                [header,par] = self.modfile(fname+str(i)+'_'+str(j)+'.mod',load.n_mod_header)    
                for k in range(load.n_layer):
                    mod2grid[i,j,k] = par[k][0]                                 
        
        return mod2grid
        
        print 'done'
    #==============================================================================
    def noisemodel2vec(self,fname):    
    #==============================================================================
        """
        Load noise model file into list an return as numpy vector

        Parameters
        ----------
        fname :             str
            name of noise model file
            
        Returns
        -------
        noisemodel :        numpy.array((n),float)
            Gaussian model with mean of zero and spread of 1. n independend 
            generated values. 
        """
        
        fid1 =open(fname,'r')
        tmp1 = fid1.readlines()
        fid1.close()
        
        noisemodel = 999999*np.ones((len(tmp1)),dtype=float)
        for i in range(len(tmp1)):
            line = tmp1[i].split()
            noisemodel[i] = float(line[0])
            
        return noisemodel
            
    #==============================================================================
    def PEST_par(self,fname,load,line=[]):                 
    #==============================================================================
        """
        load pest parameter file into numpy vector
        
        Parameters
        ----------
        fname :             str
            Filename of PEST parameter file.
        load :             obj
            Object holding all initilized values read from ini files.
        line :              List[n](int)
            python list. Lines to read from PEST parameter fil.
        
        Returns
        -------
        par_vec :           numpy.array((n),float)
            
        """
        
        # --- read all lines ----
        fid1 = open(fname,'r')
        lines = fid1.readlines()        
        fid1.close()

        # --- allocate array ----
        par_vec = np.zeros((len(line)),dtype=float)       
        
        # --- read values from "line" ---
        count1 = 0       
        for i in line:

            tmp1 = lines[i+1].split()
            par_vec[count1] = float(tmp1[1])
            count1 += 1
        
        return par_vec
        
#==============================================================================
    def strat(self,fname,load):
#==============================================================================
         """
         Load strat file holding categorical lithology units 
         
         Parameters
         ----------
         fname :            str
             String with file name of strat file
         load :             obj
             Object holding all initilized values read from ini files
         
         Returns
         -------
         mat :      numpy.array((ny,nx,n_layer),int)
            3D array with lithologies.
         """
         
         mat = np.zeros((load.ny,load.nx,load.n_layer),dtype=int)
         # --- load strat file into numpy array ---
         fid1 = open(fname,'r')
         for i in range(load.n_layer):
             for j in range(load.ny):
                 tmp = fid1.readline().split()
                 mat[j,:,i]  = map(abs,np.array(tmp,dtype=int))
         fid1.close()
         
         return mat

#==============================================================================
    def zns2mat(self,fname,load): 
#==============================================================================        
        """
        Load modlflow zone file and store in numpy matrice
        
        Parameters
        ----------
        fname :             str
            modflow zone array file 
        load :             obj
            Object holding all initilized values read from ini files
        
        Returns
        -------
        mat :               numpy.array((ny,nx,n_layer),int)
            3D array with MODFLOW zone values 
            
            
        """
        
        # --- allocate numpy array for zone model ---
        mat = np.zeros((load.ny,load.nx,load.n_layer),dtype=int)
        
        # --- load modflow zone file ---
        fid1 = open(fname,'r')
        lines = fid1.readlines()
        
        data = []
        for i in range(len(lines)):
            tmp = lines[i].split()
            if len(tmp) < 5:
                pass
            else:
                for j in range(len(tmp)):
                    data.append(tmp[j])
        # --- store modflow zone file in array ---
        count = 0
        for i in range(load.n_layer):
            for j in range(load.ny):
                for k in range(load.nx):                
                    mat[j,k,i] = data[count]
                    count += 1
        return mat
        
# *****************************************************************************
# *****************************************************************************
class write_file:        
    """
    """
#==============================================================================
    def batch_beopest(self,addr,load):
#==============================================================================
        """
        Writes BATCH-script for running beoPEST on n-slaves (machines)
        
        Parameters
        ----------
        addr :              str
            IP adress of "slave"-machine
        load :              obj
            Object holding all initilized values read from ini files              

        """
        fid1 = open('run_beopest.bat','w')
        fid1.write('\n')
        fid1.write('REM Start n slaves'+'\n')
        fid1.write('\n')
        fid1.write('for /L %%i in (1,1,'+str(load.n_slaves)+') do ('+'\n')
        fid1.write('	cd slave%%i'+'\n')
        fid1.write('	start beopest64.exe '+load.PEST_name+'.pst /h '+str(addr)+':'+str(load.IP_port)+'>> beopest_output.dat'+'\n') 
        fid1.write('	cd..'+'\n')
        fid1.write(')'+'\n')
        fid1.write('\n')
        fid1.write('beopest64.exe '+load.PEST_name+'.pst /h :'+str(load.IP_port)+'\n')
        fid1.write('\n')
        fid1.write('exit'+'\n')
        fid1.close()
        
#==============================================================================
    def bud2smp_in(self,output_file,modflow_file_name,gs_file,n_layer,ibound_file): 
#==============================================================================
        """
        Write in-file for PEST utility bud2smp. bud2smp are reading MODFLOW 
        budget(often binary) file and writes the results to a asci file.
        
        Parameters
        ----------        
        output_file : str
            name of bud2smp ini-file
        modflow_file_name:  str
            name of MODFLOW file. 
        gs_file :   str
            grid specification file
        n_layer :   int
            number of layers in model
        i_bound :   str
            name of MODFLOW boundary array
        """
        
        fid1 = open(output_file+'.in','w')
        fid1.write(gs_file+'\n')
        fid1.write(str(n_layer)+'\n')
        fid1.write(modflow_file_name+'.cbb'+'\n')
        fid1.write(str(9)+'\n')
        fid1.write(str(3)+'\n')
        fid1.write('river'+'\n')                
        fid1.write('1/1/2000\n')
        fid1.write('00:00:00\n')
        fid1.write('s\n')
        for i in range(n_layer):
            fid1.write(ibound_file+str(i+1)+'.inf'+'\n')
        fid1.write('flow'+'\n')
        fid1.write('non-flow'+'\n')
        tmp = output_file.split('.')
        fid1.write(tmp[0]+'.smp'+'\n')
        fid1.write(str(1)+'\n')
        fid1.write('f'+'\n')
        fid1.write(tmp[0]+'_rec.dat'+'\n')
        fid1.close()
        
#==============================================================================
    def bor_file(self,fname,kote_bottom,kote_top,text):
#==============================================================================
        """
        Write bor file for caly fraction modeling
        
        Parameters
        ----------
        fname : str
            output file name for glm file
        kote_bottom : int
            bottom interval of bor-file
        kote_top :  int
            top interval of bor-file
        text :  numpy.array            
            homogen matric list with text to be written in file
        """

        fid1 = open(fname+'_'+str(kote_bottom)+'_'+str(kote_top)+'.bor','w')
        n_boreholes = len(text)
        fid1.write('Test_borehole_file'+'\n')
        fid1.write(str(n_boreholes)+'\t'+str(kote_top)+'\t'+str(kote_bottom)+'\n')
        for i in range(len(text)):
            fid1.write('%15.2f %15.2f %15.2f %15.2f %15.2f %15.2f %15.12s \n' % (text[i]))
        fid1.close()


#==============================================================================
    def glm_file(self,fname,pp_x,pp_y,model,load,utm_x=0,utm_y=0,utm_z=0):
#==============================================================================
        """
        Write glm file for caly fraction modeling
        
        Parameters
        ----------
        fname:  str
            output file name for glm file
        pp_x :              List(int)  
            python list - node positions of pilot points in column direction    
        pp_y :              List(int)
            python list - node positions of pilot points in row direction
        model : numpy.array(nrow,ncol,nlayer)
            3D array with resistivities for each cell 5 values is specifyed
        load :             obj
            Object holding all initilized values read from ini files
        utm_x :   int
            definde coordinate system. Southern-east corner in depth of n_layer is basis
        utm_y :   int
            definde coordinate system. Southern-east corner in depth of n_layer is basis
        utm_z :   int
            definde coordinate system. Southern-east corner in depth of n_layer is basis
        
        """
        print ' write glm-file:              ',        
        
        text    = 'skyTEM in HYTEB' 
    
        fid1 = open(fname+'.glm','w')
        fid1.write('glm file \n')
        fid1.write(str(len(pp_x)*len(pp_y))+'\n')
        for i in pp_y:
            for j in pp_x:
                fid1.write('%20.2f  %20.2f %20.2f %1.40s \n'% (utm_x+(j)*load.dx-load.dx/2.0,utm_y+(load.ny-i)*load.dy-load.dy/2.0,utm_z+load.n_layer*load.dz,text))
                fid1.write('%1.0i  %2.1i %5.2i \n'%(load.n_layer,0,model[i,j,0,4]))
                depth = 0
                for k in range( len(model[0,0,0,:])):
                    for layer in range(load.n_layer):
                        if layer < load.n_layer:
                            if k == 0: 
                                fid1.write('%20.2f  %20.2f %20.2f \n' %(model[i,j,layer,0],model[i,j,layer,2],model[i,j,layer,3]))
                        if layer < load.n_layer-1:
                            if k == 1:
                                fid1.write('%20.2f  %20.2f \n'%(model[i,j,layer,1],model[i,j,layer,2]))
                            if k == 2:
                                depth += model[i,j,layer,1] 
                                fid1.write('%20.2f  %20.2f \n'%(depth,model[i,j,layer,1])) 
            fid1.close()                
        print 'done'
#==============================================================================
    def modfile(self,fname,h_text,mod_par):
#==============================================================================
        """
        Writes AarhusINV model files (mod-files)
        
        Parameters
        ----------
        fname : str
            Output file name of mod-file
        h_text :    str
            Header text of mod-file. Each line is given as a "indgang" in a list format
        mod_par :   float
            Model parameters and constrains in the mod-file.     

        """        
        
        fid1        = open(fname,'w')        
        
        for i in range(len(h_text)):
            fid1.write(h_text[i]+'\n')
        
        for i in range(len(mod_par)):
            n_col = len(mod_par[i])
            for j in range(n_col):
                fid1.write(str(mod_par[i][j])+'\t')
            fid1.write('\n')
        fid1.close()

        
#==============================================================================
    def modfile_clay_fraction(self,fname,bor_name_list,Nx,Ny,load,upper=[],lower=[],utm_x=0.0,utm_y=0.0,utm_z=0.0,n_iter=50):        
#==============================================================================
        """
        Write mod file for clay fraction modeling 
        
        Parameters
        ----------
        fname : str
            Output file name of mod-file
        bor_name_list : str
            Header text of mod-file. Each line is given as a "indgang" in a list format
        Nx :    int
            Model parameters and constrains in the mod-file.
        Ny :    int
            Model parameters and constrains in the mod-file.
        utm_x :   int
            definde coordinate system. Southern-east corner in depth of n_layer is basis
        utm_y :   int
            definde coordinate system. Southern-east corner in depth of n_layer is basis
        utm_z :   int
            definde coordinate system. Southern-east corner in depth of n_layer is basis
        n_iter :    int
              number of iterations        
        """ 
        
        if upper ==[] or lower == []:
            upper.append(70.0)
            lower.append(30.0)
        prior   = -1.00E+000 
        latx    = 2.00E+000
        laty    = 2.00E+000
        ver     = 0.5
        
        fid1 = open(fname+'.mod','w')
        fid1.write('%1.30s %74.60s \n'                       % ('header text','!header text'))
        fid1.write('%4.1i  %5.3i %21.30s %77.60s \n'         % (3,700,'40exp(400)','!3->kriging, Search Radius, Variogram'))
        fid1.write('%4.0i  %5.0i  %111.50s \n'               % (2,1,'!Model Type (1->linear cut-off,2->compl.err.func.3->err.func.), Kriging factor (not used)'))        
        fid1.write('%4.0i  %5.0i %112.50s \n'                % (4,20,'!CalcGrid, ClayThickness over/under grids (only if CalcGrid=1 (under) or 2 (over)) '))
        fid1.write('%4.0i  %95.60s \n'                       % (load.n_layer,'!number of datafiles (.bor)'))
        for i in bor_name_list:
            fid1.write('%15.40s %75.50s \n'                  % (i,'!borehole info layer #'))
        fid1.write('%1.60s %93.50s \n'                      % (fname+'.glm','!File with geoph. layered models'))    
        fid1.write('%7.0i %71.60s \n'                        % (n_iter,'!NIte'))
        fid1.write('%7.0i %2.0i %84.60s \n'                  % (Nx,Ny,'!Ny,Nx for all layers'))
        fid1.write('%7.1i %6.0i %80.60s \n'                  % (utm_x,utm_x+load.x,'!x1,x2 for all layers'))
        fid1.write('%7.1i %6.0i %80.60s \n'                  % (utm_y,utm_y+load.y,'!y1,y2 for all layers'))
        for i in range(len(bor_name_list)):
            tmp1= bor_name_list[i].split('.')
            tmp2= tmp1[0].split('_')
            bot = float(tmp2[-2])
            top = float(tmp2[-1])
            fid1.write('%3.1f %3.1f %85.60s \n'                  % (top,bot,'!top,bottom for layer #'))
            for j in range(2*Ny*Nx):
                if j%2 == 0:
                    fid1.write('%6.2e %2.2e  %2.2e %2.2e %2.1f %62.60s \n'       %(lower[0],prior,latx,laty,ver,'!value, prior, latx, laty, ver'))   
                else:
                    fid1.write('%6.2e %2.2e  %2.2e %2.2e %2.1f %62.60s \n'       %(upper[0],prior,latx,laty,ver,'!value, prior, latx, laty, ver'))   
        fid1.close()
        
#==============================================================================
    def modflow_zone(self,fname,mat,load):
#==============================================================================
        """
        writes modflow zone file
        
        Parameters
        ----------
        fname :             str
            String with moflow zone file name
        mat :               array_like, int
            matrice holding a 3D cathegorical model (x,y,z)
        load :              obj
            object holding all initilized values read form ini files
        
        Returns
        -------
        fname.zns
        
        """        
        # --- write modflow zone file ---
        fid1 = open(fname+'.zns','w')
        fid1.write(str(load.n_layer)+'\n')
        for layer in range(load.n_layer):
            fid1.write('ZLAY_'+str(layer+1)+'\n')
            fid1.write('INTERNAL 1 (FREE) -1'+'\n')
            for i in range(load.ny):
                for j in range(load.nx):
                    fid1.write(str(mat[i,j,layer])+' ')
                fid1.write('\n')
        fid1.close()
    
#==============================================================================
    def mod_tem_q3D(self,nx_pos,ny_pos,dim,load,res_ver_con):
#==============================================================================
        """   
        Makes mod files for AarhusINV. \n 
        The script loads in the groundwater model(GWM) grid and makes a \n
        N - layer model file based on the GWM - grid. \n
        The Foot-print of the EM signal is calculated by Pythagoras. \n 
        The hydraluic conductivities is transformed to resistivities and the 
        average res within the foot print in ease depth is writen to a mod-file.\n
         
        Parameters
        ----------
        nx\_pos :             array\_like(int)
            python list type, containt "east-west"- node positions
        ny\_pos :             array\_like(int)
            python list type, containt "north-south"- node positions
        load :              obj
            python object containing parameters values read from ini files        
        res\_ver\_con :       array\_like, float
            python list type, vertical constarins on resistivity         
        dim :              str, (q3D,1D)
            defulat set to quasi 3D. Can be changed to 1D.
        
        Returns
        -------
        *.mod
        
        Notes
        -----
        Last Change:
        16-07-2013
        Nikolaj Kruse Christensen
         
        
        """

        import math as mat
        import random as ran
        
        # -----------------------------------------------------------

        # --- resisivity constrains ----
        res_STD         = -1*np.ones(load.n_layer)
        
        # -- vertical constrains ---
        res_ver_STD = res_ver_con

            
        
        # --- thickenss ---
        thickness_STD   =  1e-3*np.ones(load.n_layer+1)
        thk_ver_STD     =  9e9*np.ones(load.n_layer+1)
        
        # --- 1D model dis in z-direction ---
        thickness       = load.dz*np.ones(load.n_layer+1)
    
        # -----------------------------------------------------------
        # --- load k-field and ibound in to matrix ---
        # -----------------------------------------------------------
    
        k_field = np.zeros( (load.ny,load.nx,load.n_layer) )
        for i in range(load.n_layer):
            fid1            = open(load.model_name+str(i+1)+'.'+load.dtype_model,'r')
            fid2            = open(load.ibound+str(i+1)+'.'+load.dtype_ibound,'r')
            k_data          = np.array(fid1.read().split(),dtype=float).reshape(load.ny,load.nx)
            i_field         = np.array(fid2.read().split(),dtype=float).reshape(load.ny,load.nx)
            i_field         = np.absolute(i_field)
            k_field[:,:,i]  = k_data*i_field
            fid1.close()
            fid2.close()
    
        # -----------------------------------------------------------
        # --- calculat radius of foot print of EM signal in depth ---
        # -----------------------------------------------------------
    
        depth = np.linspace(load.dz-load.dz/2.0,load.n_layer*load.dz-load.dz/2.0,load.n_layer)
        fp = depth*mat.tan(float(load.tem_angel)/180.0*mat.pi)
        
        # -----------------------------------------------------------
        #           --- sample resistivities ----
        # The Avg_res for each layer is caluclated in log10(res)
        # look in layer around position(i,j) smaller than foot-print
        # of EM signal.
        # 
        # -----------------------------------------------------------
    
    
        cells = fp/load.dx
        
        t = np.floor(cells)
        
        res_layer = []
        for layer in range(load.n_layer):
            num_cells = int(t[layer])

            # --- sample the 1D model as 1D, oherwise as quasi 3D ---            
            if dim == '1D':
                num_cells = 0
                
            count = 0
            con = 0
            if num_cells > 0:
                for k in range(-num_cells,num_cells+1,1):
                    for h in range(-num_cells,num_cells+1,1):
                        r = ((k*load.dx)**2+(h*load.dx)**2)**(1/2)
                        if r < fp[layer]:
                            hk = k_field[ny_pos+k,nx_pos+h,layer]
                            if hk == 0.0:
                                std     =  mat.log10(1+load.noise_petro)
                                e       =  ran.gauss(0,std)
                                con += np.log10(load.res_clay)+e/load.beta_2
                            else:
                                std     =  mat.log10(1+load.noise_petro)
                                e       =  ran.gauss(0,std)
                                con     += (mat.log10(hk)-mat.log10(load.beta_1)-e)/load.beta_2
                            count += 1.0
            elif num_cells == 0.0:
                r = ((load.dx)**2+(load.dx**2))**(1/2)
                if r < fp[layer]:
                    hk = k_field[ny_pos,nx_pos,layer]

                    if hk == 0.0:
                        std     =  mat.log10(1+load.noise_petro)
                        e       =  ran.gauss(0,std)
                        con += np.log10(load.res_clay)+e/load.beta_2
                    else:
                        std     = mat.log10(1+load.noise_petro)
                        e       = ran.gauss(0,std)
                        con     += (mat.log10(hk)-mat.log10(load.beta_1)-e)/load.beta_2
                    count += 1.0
            else:
                print 'Error in Writefile.py --> write_mod_tem_q3D'
           
            res_layer += [(con/count)]
         
        # -----------------------------------------------------------
        #               --- writing mod file ----
        # -----------------------------------------------------------
        fid1 = open(load.mod_header,'r')
        fid2 = open(load.mod_file_temp+str(ny_pos)+'_'+str(nx_pos)+'.mod','w')
    
        for i in range(3):
            tex = fid1.readline()
            fid2.write(tex) 
    
        fid2.write(str(load.niter_mod_temp_gp)+'\n')
        fid2.write(str(load.n_layer+1)+'\n')
        
        for i in range(load.n_layer):
            fid2.write('%4.7f %15.3e %15.3e\n' % (10**res_layer[i], res_STD[i],res_ver_STD[i]))
        fid2.write('%4.7f %15.3e\n' % (load.res_clay, res_STD[0]))
    
        for i in range(load.n_layer):
            fid2.write('%1.1f %20.3e %20.3e\n' % (thickness[i],thickness_STD[i],thk_ver_STD[i]))
            
        for i in range(load.n_layer):
            fid2.write('%1.1f %13.1i\n' % (depth[i]+load.dz/2.0, -1))
    
        fid2.close()
        fid1.close()
        
#==============================================================================
    def mod_voxel(self,fname,pp_x,pp_y,load,grid_file='dummy_DEM.grd',altitude=True,run_type='dummy'):
#==============================================================================
        """
        write AarhusInv mod-file for voxel inversion
        
        Parameters
        ----------
        fname : str
            output file name of AarhuINV mod file for voxel inversion
        pp_x :              List(int)  
            python list - node positions of pilot points in column direction
        pp_y :              List(int)
            python list - node positions of pilot points in row direction
        load :             obj
            Object holding all initilized values read from ini files
        grid_file : str,default('dummy_DEM.grd')
            DEM(elevation)-grid
        altitude :  logical,default(True)
            True  --> include flight altitude in inversion
            None  --> do not include flight altitude in inversion
        run_type :  str, default('dummy')
            dummy --> use dummy/template TEM file
            
        
        """
        # --- initial resistivity model ---
        initial_res = 100*np.ones((len(pp_y),len(pp_x),load.n_layer),dtype=float)
        if altitude == True:
            ModelType = 5
            MT_norms  = 5
        if altitude == False:
            ModelType = 1        
            MT_norms  = 3
            
        # --- write mod file ---
        fid1 = open(fname+'.mod','w')
        fid1.write('%1.2i %3.2i %3.1f %3.60s  '%(load.dx,load.dy,0.0,grid_file))
        for i in range(1,MT_norms+1,1):        
            fid1.write('%3.1i ' % (i))                
        fid1.write('%3.2f %3.2f %3.2f %3.2f \n '% (load.dy/2.0,load.x-load.dx/2.0,load.dy/2.0,load.y-load.dy/2.0))
        fid1.write('%1.1i %2.1i %3.60s \n' % (2*len(pp_x)*len(pp_y),2,'!number of files, constarit width , sharpness settings'))         
        count1 = 0
        for i in pp_y:        
            for j in pp_x:
                count1 +=1
                for h in range(2):
                    if h == 0:
                        if run_type != 'dummy':
                            fid1.write('%1.1i %7.1i %25.60s \n' %(count1,ModelType,load.tem_name+str(i)+'_'+str(j)+'_low.tem'))
                        else:
                            fid1.write('%1.1i %7.1i %25.60s \n' %(count1,ModelType,'tem_example'+'_low_moment.tem'))
                    elif h== 1:
                        if run_type != 'dummy':                        
                            fid1.write('%1.1i %7.1i %25.60s \n' %(count1,ModelType,load.tem_name+str(i)+'_'+str(j)+'_high.tem'))
                        else:
                            fid1.write('%1.1i %7.1i %25.60s \n' %(count1,ModelType,'tem_example'+'_high_moment.tem'))
        fid1.write('%1.1i %25.60s \n' %(load.niter_mod_gp,'! number of iterations' ))
        count1 = 0
        for i in pp_y:        
            count2 = 0
            for j in pp_x:
                fid1.write('%1.1i %10.1f %10.1f %10.1f %25.60s \n' %(load.n_layer,j*load.dx-load.dx/2.0,i*load.dy-load.dy/2.0, 0,'   !number of layers, x, y, z'))
                if altitude == True:
                    fid1.write('%1.1i %10.2f %10.2f %25.60s \n' %(ModelType,0.04,0.3,''))
               
               # --- resisisvity section; res, STD_prior, vertical_con, horizontal_con ---    
                ver_con = np.linspace(1,0.2,load.n_layer-1) 
                hor_con = np.linspace(0.1,0.3,load.n_layer) 
                
                for layer in range(load.n_layer):
                    fid1.write('%1.2i    ' %(initial_res[count1,count2,0]))
                    fid1.write('%10.1f   ' %(-1))
                    if layer < load.n_layer-1:
                        fid1.write('%10.1f   ' %(ver_con[layer]))
                        if j != pp_x[-1]:
                            fid1.write('%10.1f \n' %(hor_con[layer]))
                        if j == pp_x[-1]:
                            fid1.write('\n')                        
                    else:
                        if j != pp_x[-1]:
                            fid1.write('%23.1f \n' %(hor_con[layer]))
                        if j == pp_x[-1]:
                            fid1.write('\n')
                    
                                         
                # --- thickness section; res, STD_prior, vertical_con, horizontal_con ---    
                for layer in range(load.n_layer-1):    
                    fid1.write('%1.2i    ' %(load.dz,))
                    fid1.write('%11.3e   ' %(0.001))
                    if layer < load.n_layer-2:
                        fid1.write('%10.1f   ' %(-1))
                        if j != pp_x[-1]:
                            fid1.write('%10.1f \n' %(-1))
                        else: 
                            fid1.write('\n')
                    else:
                        if j != pp_x[-1]:
                            fid1.write('%23.1f \n' %(-1))
                        else:
                            fid1.write('\n')
                # --- depth section; res, STD_prior, vertical_con, horizontal_con ---                    
                depth = 0.0
                for layer in range(load.n_layer-1):
                    depth += load.dz
                    fid1.write('%1.2i    ' %(depth))
                    fid1.write('%11.3e   ' %(-1))
                    if j != pp_x[-1]:
                        fid1.write('%23.1f \n' %(-1))
                    if j == pp_x[-1]:
                        fid1.write('\n')
                    
                        
                    
                count2 += 1
            count1 += 1
    
                
        fid1.close()
           
#==============================================================================
    def mod2smp_in(self,output_file,modflow_file_name,gs_file,bc_file):
#==============================================================================
        """
        Writes in-file for PEST utility MOD2SMP
        
        Parameters
        ----------
        output_file :   str
            name of MOD2SMP in-file
        modflow_file_name:  str
            name of MODFLOW file. 
        gs_file :   str
            grid specification file
        bc_file:    str
            bore coordinate file
        """
        
        fid1 = open(output_file,'w')
        fid1.write(gs_file+'\n')
        fid1.write(bc_file+'\n')
        fid1.write(bc_file+'\n')
        fid1.write(modflow_file_name+'.hds'+'\n')        
        fid1.write('f\n')
        fid1.write('1\n')
        fid1.write('1e30\n')
        fid1.write('s\n')
        fid1.write('1/1/2000\n')
        fid1.write('00:00:00\n')
        tmp = bc_file.split('.')
        fid1.write(tmp[0]+'.smp'+'\n')
        fid1.close()
        
#==============================================================================
    def mod2obs_in(self,outputfile,modflow_file_name,gs_file,bc_file,n_layer):         
#==============================================================================
        """
        Writes in-file for PEST utility MOD2OBS
        
        Parameters
        ----------
        output_file :   str
            name of MOD2OBS in-file
        modflow_file_name:  str
            name of MODFLOW file. 
        gs_file :   str
            grid specification file
        bc_file:    str
            bore coordinate file
        n_layer :   int
            number of layers in groundwater model
        """
        
        #  --- write head ini file for mod2obs ---
        fid1    = open(outputfile+'.in','w')        
        fid1.write(gs_file+'\n')
        fid1.write(bc_file+'\n')
        fid1.write(bc_file+'\n')
        tmp = bc_file.split('.')
        fid1.write(tmp[0]+'.smp'+'\n')
        fid1.write(modflow_file_name+'.hds'+'\n')
        fid1.write('f'+'\n')
        fid1.write('1e30'+'\n')
        fid1.write('s'+'\n')
        fid1.write('1/1/2000'+'\n')
        fid1.write('00:00:00'+'\n')
        fid1.write(str(n_layer)+'\n')
        fid1.write('3.5e-5'+'\n')
        fid1.write(tmp[0]+'_model.smp'+'\n')
        fid1.close()        
                        
#==============================================================================
    def mod2tpl(self,nx_pos,ny_pos,load):
#==============================================================================
        """
        Writes PEST template files for AarhusINV model files (mod-file)
        Reads a mod-file for each geophysical measurement position and
        writes a corrosponding template file
        
        Parameters
        ----------
        nx_pos :
            node positions of pilot points in column direction
        ny_pos :
            node positions of pilot points in row direction
        load :             obj
            Object holding all initilized values read from ini files
            
        Returns
        -------
        "mod_name"+position.tpl
         
        """

        
        #  -----------------------------------------------------------
        #               --- writing mod tpl file ----
        # -----------------------------------------------------------
        fid1 = open(load.mod_file+str(ny_pos)+'_'+str(nx_pos)+'.mod','r')
        fid2 = open(load.mod_file+str(ny_pos)+'_'+str(nx_pos)+'.tpl','w')

        
        fid2.write('ptf #\n')
        
        count1 = 0
        count2 = 1
        while 1:
            tex = fid1.readline()
            if count1 < load.n_mod_header:
                fid2.write(tex)
            if count1 >=load.n_mod_header and count1 < load.n_mod_header+load.n_layer:
                tmp = tex.split()
                fid2.write('#r_'+str(ny_pos)+'_'+str(nx_pos)+'_'+str(count2)+'    #  '+tmp[1]+'  '+tmp[2]+'\n')
                count2 += 1
            if count1 == load.n_mod_header+load.n_layer:
                tmp = tex.split()
                fid2.write('#r_'+str(ny_pos)+'_'+str(nx_pos)+'_'+str(count2)+'    #  '+tmp[1]+'\n')
            if count1 >= load.n_mod_header+load.n_layer+1:
                fid2.write(tex)
                
            count1 +=1
            
            if not tex:
                    break
            
        fid1.close()
        fid2.close()

#==============================================================================
    def pts_file(self,pp_file_name,pp_name,x_pos,y_pos,flag_1,hk): 
#==============================================================================
        """
        Write PEST pilot point file
         
        Parameters
        ----------
        pp_file_name :  str
            pilot point file name
        pp_name :   str
        x_pos :     List[]
        y_pos :     List[]
        flag_1 :    List[]
        hk :        List[]
        

        Notes
        -----        
        12-08-2014
        Nikolaj Kruse Christensen
        nkc07@phys.au.dk
        
        """ 
        
        fid1 = open(pp_file_name+'.pts','w')
        count = 0
        for j in range(len(y_pos)):
            for i in range(len(x_pos)):
                if i < 9: 
                    fid1.write('%1.10s %11.6s %10.6s %4.1i %18.6e\n' % (pp_name+str(count+1), x_pos[i], y_pos[j],flag_1[count],hk[count]))
                if 9 <= i < 99: 
                    fid1.write('%1.10s %10.6s %10.6s %4.1i %18.6e\n' % (pp_name+str(count+1), x_pos[i], y_pos[j],flag_1[count],hk[count]))
                if 99 <= i < 999: 
                    fid1.write('%1.10s %9.6s %10.6s %4.1i %18.6e\n' % (pp_name+str(count+1), x_pos[i], y_pos[j],flag_1[count],hk[count]))
                count += 1
        
        fid1.close()

#==============================================================================
    def pts2tpl(self,layer,pp_file_name,left): 
#==============================================================================
        """      
        Write PEST pilot point template file
         
        Parameters
        ----------
        layer : int
            number of layers in groundwater model
        pp_file_name :  str
            pilot point file name        
        left :  int
            left indent for start reading
  
        Notes
        -----
        12-08-2014
        Nikolaj Kruse Christensen
        nkc07@phys.au.dk
        
        """
        #  -----------------------------------------------------------
        #               --- writing tpl file ----
        # -----------------------------------------------------------
        fid1 = open(pp_file_name+str(layer)+'.pts','r')
        fid2 = open(pp_file_name+str(layer)+'.tpl','w')
        fid2.write('ptf #\n') 
        name = pp_file_name.split('_')
        count = 1        
        while 1:
            tex = fid1.readline()
            if not tex:
                break
            else:                
                fid2.write(tex[0:left]+'  #'+name[0]+str(layer)+'_ppt'+str(count)+'                                      # \n')            
            count += 1 
               
        fid2.close()
        fid1.close()        
       

        
#==============================================================================
    def sortbyept(self,endpoint_file,nx,ny,n_par,well_x_node,well_y_node,layer):        
#==============================================================================
        """
        Write input file for small program sortbyept.exe
        sortbyept read a MODPATH endpoint file and sorting particles going to 
        a pumping well or going to a lake.
        
        Parameters
        ----------
        endpoint_file : str
            name of MODPATH endpoint file
        nx :    int
            number of columns in groundwater model grid
        ny :    int
            number of rows in groundwater model grid
        n_par : int
            number of particle of grid cell. 
        well_x_node :    int
            column number for location of pumping well
        well_y_node :    int
            row number for location of pumping well
        layer : int
            layer where the pumping well is screened
        """
        
        fid1 = open('sortbyept.inp','w')
        fid1.write(endpoint_file+'      MODPATH endpoint-file \n')
        fid1.write(str(ny)+' '+str(nx)+' '+str(n_par)+'      Number of rows and columns in grid, and number of particles per grid cell'+'\n')
        fid1.write(str(well_x_node)+' '+str(well_y_node)+' '+str(layer)+'\n')  
        fid1.write('in-well.dat '+'\n')
        fid1.write('in-sea.dat'+'\n')
        fid1.close()
        
        

#==============================================================================
    def vox(self,fname,ModelMeshName,ForwardMeshName,ConstraintMeshName,load,pp_x,pp_y,altitude=True,itarations=0,run_type='forward'):
#==============================================================================
        """       
        write AarhusINV voxel inputfile file

        Parameters
        ----------
        fname :
        ModelMeshName :     str
            List with vtk mesh for resitivity and altitude
        ForwardMeshName :   str
            name of forward mesh
        ConstraintMeshName :    str
            name of constrian mesh
        load :             obj
            Object holding all initilized values read from ini files
        pp_x :              List(int)  
            python list - node positions of pilot points in column direction 
        pp_y :              List(int)
            python list - node positions of pilot points in row direction
        itarations :    int
            number of iterations
        run_type :  str
            forward   --> pointing to "dummy.tem" file for skyTEM setup.
            inversion --> pointing to forward skyTEM data file (.tem).
            
        Returns
        -------
        fname.vox

        Notes
        -----
        Last change 2015-06-19
        

        """
        
        print ' AarhusInv vox file       ',        
        
        # --- %1 arguments ---
        if altitude == True:
            Models          = 1
            ModelMeshes     = 2
            Interpolations  = 2
            Norms           = 6
            Ranges          = 1
            ForwardMeshes   = 1
            ModelType       = 5
            NormID          = 6
        
        elif altitude == False:
            Models          = 1
            ModelMeshes     = 1
            Interpolations  = 1
            Norms           = 4
            Ranges          = 1
            ForwardMeshes   = 1
            ModelType       = 1
            NormID          = 4

        
        fid1 = open(fname+'.vox','w')
        fid1.write('% HYTEB voxel inversion \n')
        fid1.write('%1  #Models  #ModelMeshes  #Interpolations  #Norms    #Ranges   #ForwardMeshes  #Itarations \n')
        fid1.write('%11.1i %13.1i %16.1i %7.1i %10.1i %16.1i %12.3i \n' % (Models,ModelMeshes,Interpolations,Norms,Ranges,ForwardMeshes,itarations))
        fid1.write('%2      ModelMeshID ModelID #Regions  Connectivity  ModelType      #Parameters     ParTypes   LogLin   RangeIDs  InterpolationIDs  RegionSettings   ModelMeshName \n')
        fid1.write('                  1       1        1             1          '+str(ModelType)    +'                1            1        0          1                 1             0 0      '+ModelMeshName+' \n')
        if altitude == True:
            fid1.write('                  2       1        1             0          '+str(ModelType)+'                1            3        0          1                 2             0 0      '+ModelMeshName+' \n')
        
        fid1.write('%3    ForwardMeshID   Dimensionality                ModelType ConnectedModelID   #DataFiles                       ForwardMeshName \n')
        fid1.write('                  1                1                        '+str(ModelType)+'                1            '+str(2*len(pp_x)*len(pp_y))+'                          ' +ForwardMeshName+' \n')
        fid1.write('%4 ConstraintMeshID                                                      ConstraintMeshName \n')
        fid1.write('                  1                                                              '+ConstraintMeshName+' \n')
        if altitude == True:
            fid1.write('                  2                                                              '+ConstraintMeshName+' \n')
        
                          
        fid1.write('%5           NormID           NormType           NormSettings                      NormName \n')
        fid1.write('                  1                  1                                                   L2  \n')
        fid1.write('                  2                  1                                                   L2 \n')
        fid1.write('                  3                  1                                                   L2 \n')
        fid1.write('                  4                  1                                                   L2 \n')
        if altitude == True:
            fid1.write('                  5                  1                                                   L2 \n')
            fid1.write('                  6                  1                                                   L2 \n')                 
        fid1.write('%6  InterpolationID  InterpolationType  InterpolationSettings             InterpolationName \n')
        fid1.write('                  1                  2         1 4  100.0 2.0              inverse_distance \n')             
        if altitude == True:
            fid1.write('                  2                  1                                    Nearest_neighbour \n')     
                
        fid1.write('%7          RangeID           MinValue               MaxValue  \n')               
        fid1.write('                  1            default                default \n')                
        fid1.write('%8   ForwardMeshID       ForwardMeshSubID    DataType      NormID  LogLin       DataFileName \n')
        count1 = 0
        for i in pp_y:
            for j in pp_x:   
                count1 +=1
                for h in range(2):                   
                    if h == 0:
                        if run_type == 'inversion':
                            fid1.write('%17.1i %17.1i %17.1i %11.1i %7.1i %25.60s \n' %(1,count1,1,NormID,0,load.tem_name+str(i)+'_'+str(j)+'_low.tem'))
                        elif run_type == 'forward':
                            fid1.write('%17.1i %17.1i %17.1i %11.1i %7.1i %25.60s \n' %(1,count1,1,NormID,0,'tem_example_low_moment.tem'))
                    elif h== 1:
                        if run_type == 'inversion':
                            fid1.write('%17.1i %17.1i %17.1i %11.1i %7.1i %25.60s \n' %(1,count1,1,NormID,0,load.tem_name+str(i)+'_'+str(j)+'_high.tem'))
                        elif run_type == 'forward':
                            fid1.write('%17.1i %17.1i %17.1i %11.1i %7.1i %25.60s \n' %(1,count1,1,NormID,0,'tem_example_high_moment.tem'))
        fid1.close()
        
        print 'done'
        

    #==============================================================================
    def integer_array_file(self,fname,zone_mat):
    #==============================================================================
        """
        Writes MODFLOW integer array for def zones of groundwater model
        
        Parameters
        ----------
        fname : str
            output name of integer array files
        zone_mat :  numpy.array(int)
            matrice def the zones of the groundwater model
        """            
        
        # --- size of one array matrice ---
        [n_row,n_col,n_layer] = np.shape(zone_mat)
        
        # --- write inter array for each layer
        for layer in range(n_layer):
            tmp = fname.split('.')            
            outputname = tmp[0]+'_'+str(layer+1)+'.'+tmp[-1]
            fid1 = open(outputname,'w')

            for i in range(n_row):
                for j in range(n_col):
                    fid1.write(str(zone_mat[i,j,layer])+' ')
                fid1.write('\n')
            fid1.close()
