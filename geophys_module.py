# -*- coding: utf-8 -*-
#"""
#Created on Mon Aug 11 14:35:49 2014
#
#@author: Nikolaj
#"""

import os
import numpy as np
import subprocess
import shutil

import load_write

# *****************************************************************************
# *****************************************************************************
class geophys_module:
    """
    class geophys module bla bla bla
    
    Parameters
    ----------
    model_dis : bla bla string
            bla bla
    
    """
    def __init__(self,model_dis ='model_dis.ini', hydro_ini = 'hydro_ini.ini',geophys_tem='geophys_tem.ini'):
                
        self.write  = load_write.write_file()
        self.load   = load_write.load_file()
        
        self.load.ini_model_dis(model_dis)
        self.load.ini_hydro(hydro_ini)
        self.load.ini_geophys(geophys_tem)
        
        
#==============================================================================
    def clay_fraction(self,fname,pp_x,pp_y,grid,nx,ny) :        
#==============================================================================
        """
        The petrophysical connection between hydrological parameters and 
        geophysical parameters varies spatially within survey areas and to an
        even larger degree between survey areas. This means that a global, 
        fixed petrophysical connection is non-exisisting. 
        In this modeling approach, the connection between hydrological and 
        geophysical parameters is managed by a translator function with 
        spatially variable parameters. 
        
        Parameters 
        ----------
        fname :             str
            name of clay-fraction files
        pp_x :              int
            python list - node positions of pilot points in column direction
        pp_y :              int 
            python list - node positions of pilot points in row direction
        grid :              numpy.array     
            3D voxel grid of model domain. 
        nx :                int
            number of grid-nodes in clay-fraction transformation grid, See [1]
        ny :                int
            number of grid-nodes in clay-fraction transformation grid, See [1]
        

        Returns
        -------
        fname.bor   -   Borefile        (AarhusINV)
        fname.glm   -   Geophys file    (AarhusINV)
        fname.mod   -   Mod-file        (AarhusINV)


        References
        ----------
        [1] Foged_2014
        [2] Maker_2015
        """
                       
        print 'clay-fraction modeling:'
        
        # definde coordinate system. Southern-east corner in depth of n_layer is basis
        utm_x   = 0.0
        utm_y   = 0.0
        utm_z   = 0.0

        # --------------------------------------------------------------------
        # --- load eas borehole files and write clay-fraction borhole file ---
        # --------------------------------------------------------------------              
        bor_name_list = []    
        for i in range(self.load.n_layer,0,-1):
            
            kote_bottom = (i-1)*self.load.dz
            kote_top    = i*self.load.dz
            
            # --- store borehole files in list ---           
            bor_name_list.extend([fname+'_'+str(kote_bottom)+'_'+str(kote_top)+'.bor'])        
            
            if kote_bottom > 150:
                # --- ammount of clay in different units in cap layer ---- 
                sand        = 0.1
                silt        = 0.5
                till        = 0.7
                prop_clay = [sand,silt,till]
            
            if kote_top <= 150:
                # --- ammount of clay in different units in valley ---- 
                gravel      = 0.2
                sand        = 0.4
                silt        = 0.5
                till        = 0.7
                prop_clay = [gravel,sand,silt,till]
    
            # --- read eas borehole file and store information about each kote interval ---
            name_list = ['borehole_cap.eas','borehole_valley.eas']
            text = self.load.eas_borehole(name_list,kote_bottom,kote_top,prop_clay)
            
            # --- write borefile for clay-fraction modeling ---        
            self.write.bor_file(fname,kote_bottom,kote_top,text)
        
        # -----------------------------------------------------------    
        # --- write geophys file (glm) for clay fraction modeling ---
        # -----------------------------------------------------------
        self.write.glm_file(fname,pp_x,pp_y,grid,self.load)
        
        # -----------------------------------------------------------
        # --- write mod for clay fraction modeling ---
        # -----------------------------------------------------------
        self.write.modfile_clay_fraction(fname,bor_name_list,nx,ny,self.load,utm_x=utm_x,utm_y=utm_y)
        

#==============================================================================
    def cvmesh2cvmesh(self,inputfile,outputfile,hor_con,factor=1,index=3,
                      nbytes=3,index_sharp=-1):       
#==============================================================================
        """
        Read voxel inversion constrain mesh and change horziontal constrains 
        along x and reduce the constrain by the factor
        new hor_con_x = factor*hor_con
        
        Parameters
        ----------
        inputfile :         str
        
        hor_con :           numpy.array
        
        factor :            float
        
        index :             int
            find index,
            -1 --> dummy value, "do nothing",
            2 --> vertical constrains,
            3 --> horziontal constrains,
             
        nbytes :            int
            read the nbytes+1 entry measured in bytes from the beginning of the present line
            if entry "4" is equal to "1 (1e-4)" it is a "back door" to get north-sourt" constrains 
            
        index_sharp:        int
            change constrain index 3 when running SHARP inversion ---
            -1               --> do nothing
            index_sharp > -1 --> constrain index is changed to index_sharp
            
        Returns
        -------
        outputfile :        str
        
    
        """
        
        fid1    = open(inputfile,'r')
        lines   = fid1.readlines()
        fid1.close()
        
        fid2    = open(outputfile,'w')        
        
        # --- find line in vtk file ---        
        find     = 'Res_#'
        count1   = 0

        
        while count1 < len(lines):
            # --- find and write res into vtk file
            if find in lines[count1]:
                fid2.write(lines[count1])
                tmp     = lines[count1].split() 
                n_par   = int(tmp[-2])                
                for i in range(1,n_par+1,1):
                    tmp     = lines[count1+i].split()
                    
                    # --- find v/h constrain ---
                    if int(tmp[1]) == index:
                          
                        # --- get constrain value ---
                        tmp2 = tmp[0]
                        
                        # --- find north-south constrains ---
                        if int(tmp2[nbytes]) == 1:

                            hc  = np.round(float(tmp2[0:2]+tmp2[4:9])*factor,decimals=2)  
#                            tmp[0] = str(hc+0.0001)
                            tmp[0] = str(hc+1e-4)
                            
                            # --- change N-S constrians index if for sharp inversion for ---
                            if index_sharp != -1:
                                tmp[1] = index_sharp
                    tmp = map(float,tmp)
                    fid2.write('%5.3e %10.0i %10.0i %5.0i %5.0i %5.0i %5.0i %5.0i %5.0i %15.7e \n' % (tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[5],tmp[6],tmp[7],tmp[8],tmp[9]))
                
                # --- break after npar iterations ---
                break
            
            # --- write first part of the vtk file ---
            else:
                fid2.write(lines[count1])
            count1 += 1
            
        fid2.close()
        
        print 'done'

        
        
#==============================================================================
    def emo2mod(self,inputfile,outputfile,pp_x,pp_y,STDbelowDOI=1000):
#==============================================================================
        """
        Read emo file and mod file and merge the inversion reslut (estimated parameters) into a new modfile.
        Check that the settings in the em1dinv.con file is set to print emo file
        
        
        Parameters 
        ----------
        inputfile :         str
                            mod and emo file name. The name of the mod and emo file names are identical for writing new mod-file with inversion result
        
        pp_x :              int
                            python list - node positions of pilot points in column direction
        
        pp_y :              int
                            python list - node positions of pilot points in row direction
        
        STDbelowDOI :       float
                            value to change the calculated parameter STD from AarhusINV. Value to give par_STD below less DOI; less faith [4]

        
        Returns
        -------
        outputfile :        str
                            new mod-file name
        par_est :           array_like(float)
            python list(n,z,m), 
            n = number of TEM soundings
            z = number of layers
            m = 2 (res,STD)     
        
        References
        ----------
        [3] AarhusINV - manual: http://www.hgg.geo.au.dk/HGGsoftware/em1dinv/em1dinv_manual.pdf \n
        [4] Anders Vest Christiansen - DOI
        
        """
        
        print 'emo2mod:                   ',
        
        # --- decleration of python list type ---  
        par_err = []
        par_mat = []
        doi     = []
        
        
        for i in pp_y:
            for j in pp_x:
                # --- load mod file ---      
                [h_text,mpar] = self.load.modfile(inputfile+str(i)+'_'+str(j)+'.mod',self.load.n_mod_header)        
                
                # --- load par and unc from emo file ---
                [par,err,doi_tmp] = self.load.emofile(inputfile+str(i)+'_'+str(j)+'.emo',9,self.load.n_layer)
                par_err.append(err)
                par_mat.append(par)
                doi.append(doi_tmp)


                # --- write mod-file ---
                mod = []
                for h in range(len(mpar)):
                    tmp = []
                    for k in range(len(mpar[h])):
                        if k == 0:
                            tmp += [par[h]]
                        if k > 0:
                            tmp += [mpar[h][k]]
                    mod.append(tmp)
                            
                self.write.modfile(outputfile+str(i)+'_'+str(j)+'.mod',h_text,mod)
        
        # --- make list structure with parameter estimation results---
        [ix,iy]=np.shape(par_err)
        par_est = []
        for i in range(ix):
            tmp1 = []
            for j in range(iy):
                tmp2 = []
                
                if j*self.load.dz >= doi[i]:
                    tmp2.extend([par_mat[i][j],STDbelowDOI])
                else:
                    if par_err[i][j] < 0:    
                        tmp2.extend([par_mat[i][j],STDbelowDOI]) 
                    else:     
                        tmp2.extend([par_mat[i][j],par_err[i][j]])#*par_mat[i][j]])
                tmp1.append(tmp2)    
                    
            par_est.append(tmp1)

        print 'done'
        
        return par_est
        
#==============================================================================
    def fwr2tem(self,pp_x,pp_y):
#==============================================================================       
        """
        A one-line summary that does not use variable names or the
        function name.
        Several sentences providing an extended description. Refer to
        variables using back-ticks, e.g. `var`.
        
        Parameters
        ----------
        pp_x :              int
                            python list - node positions of pilot points in column direction
        
        pp_y :              int
                            python list - node positions of pilot points in row direction

        
        Returns
        -------
        tem_name_raw+pp_y[i]_pp_x[j]+.tem

        
        References
        ----------
        [5] Munkholm & Auken 1996
        [6] Auken et al. 2008

        """
        
        print 'fwr2tem:         ',    
        
        
        
        
        # --- open on forward file to get gate center times ---
        try:
            fname = self.load.mod_file_temp+str(pp_y[0])+'_'+str(pp_x[0])
            
            if self.load.geophys_software == 'em1dinv64':
                fid1 = open(fname+'001.fwr','r')
            else:
                fid1 = open(fname+'00001.fwr','r')
                
               
        except IOError:
            print fname +' does not excist!!!'  
        
        time =[]
        count = 0
        while True:
            try:
                data = fid1.readline().split()
                if count > self.load.n_data_header-1:
                    time.append(data[0])
                count += 1
            except:
                break       
        fid1.close() 
        
        # --- Background noise contribution --- 
        V_noise = []
        for i in range(len(time)):
            V_noise.append(self.load.noise_em*(float(time[i])/1e-3)**(-1/2))
    
        # ----------------------------------------------------------------------
        # -- Write tem files with forward respons as data-file (tem) ---
        # ----------------------------------------------------------------------
        for k in pp_y: # range(1,self.n_pp_y+1,1):                                                                          
            for j in pp_x: #range(1,self.n_pp_x+1,1):
                fname = self.load.mod_file_temp+str(k)+'_'+str(j)
                try:
                    if self.load.geophys_software == 'em1dinv64':
                        fid3 = open(fname+'001.fwr','r')
                    else:
                        fid3 = open(fname+'00001.fwr','r')
                except IOError:
                    print fname +' does not excist!!!'  
                    break          
                        
                fid4 = open(self.load.tem_name_raw+str(k)+'_'+str(j)+'.tem','w')  
                
                # --- Rewrites fwr to tem_raw file and print as tem files ---
                # --- Adding Gausian and background noise to data ---                              
                if self.load.gp_noise_model == -1:
                    G = np.random.normal(0,1,len(time))
                if self.load.gp_noise_model == 0:
                    G= self.load.noisemodel2vec(self.load.gp_noise_model_file)
                                      
                count = 0
                count1 = 0
                print self.load.gp_relative_std_as_weight
                while count < len(time)+self.load.n_data_header:
#                    try:
                        if count <= self.load.n_data_header-1:
                            tex=fid3.readline()                
                            fid4.write('%s'%(str(tex)))
                            
                        if count > self.load.n_data_header-1:
                            data=fid3.readline().split()
                            # --- change the moment of the signal by adding amplitude of the current to the AarhusInv signal!!!                    
                            V = float(data[1])*self.load.current             
                            
                            # --- assigning STD as same magnitude as the noise added to data ---
                            if self.load.gp_relative_std_as_weight == 0.0:
                                sigma       = self.load.noise_gauss
                                std =((sigma)**2.0+(float(V_noise[count-self.load.n_data_header])/V)**2.0)**(1.0/2.0)
                            
                            # --- assigning STD as bigger value  that noise on data ---    
                            elif self.load.gp_relative_std_as_weight > 0.0: 
                                sigma       = self.load.gp_relative_std_as_weight
                                std     =((sigma)**2.0+(float(V_noise[count-10])/V)**2.0)**(1.0/2.0)
                                                            
                            V_resp =V + V*std*G[count1]
                            V_norm = V_resp/self.load.current
                            fid4.write('%8.5e %10.4E %3.2E %1.1i %1.1i \n' %(float(data[0]),V_norm,std,1,1))
                            count1 += 1
                        count += 1 
                        
#                    except:
#                        # --- break while loop, for totol number of TEM locations --- 
#                        if count >= len(pp_x)*len(pp_y):
#                            print 'fwr2tem  loop break', count
#                            break
                            
                fid3.close()
                fid4.close()
                
        print 'done' 

#==============================================================================
    def mod2mod(self,inputfile,outputfile,pp_x,pp_y,n_iter=None,homogen=0):
#==============================================================================
        """
        Read AarhusInv mod-files and writes new files with changed number of
        iterations, and model parameters.
        
        Parameters 
        ----------
        inputfile :         str
                    
        pp_x :              int
            python list - node positions of pilot points in column direction

        pp_y :              int 
            python list - node positions of pilot points in row direction
        
        n_iter :            int
            if n_iter== None the number of iteration in the in the new mod-file will not be changed, otherwise number of iterations = n_iter        
                            
        homogen :           float
            homogen = 0 --> mod-file will hold the "true 1D models",
            homogen > 0 --> homogen is the initial parameter values for in mod-file

  
        Returns
        -------
        outputfile :        str
                            new mod-file with changed n_iter and parametervalues. For change in other part of the mod-file change mod-file template
                                    
        """
        print 'mod2mod:                 ',

        for k in pp_y: 
            for j in pp_x:
                # --- read mod-file ---
                fid1 = open(inputfile+str(k)+'_'+str(j)+'.mod','r');
                temp = fid1.read().splitlines()
                fid1.close()
                     
                # --- declare python list type for holding  mod-header lines ---              
                text = []
                # --- loop over the number ofheader lines in mod-file ---
                for i in range(self.load.n_mod_header):
                    # --- change tem-file to tem-file with synthetic data with noise ---
                    if i == 2:
                        text.append('1 1 '+self.load.tem_name+str(k)+'_'+str(j)+'.tem')
                    # --- change number of iterations ---
                    elif i == 3:
                        if n_iter == None:
                            text.append(str(self.load.niter_mod_gp))
                        elif n_iter >= 0 or n_iter <= 0:
                            text.append(str(n_iter))   
                    # --- store mod-header text in list type ---
                    else:                        
                        text.append(temp[i])
                
                # --- declare python list type for holding mod-parameter lines ---       
                mod = []
                count = 1
                while count < self.load.n_layer*3+2:
                    try:
                        line = []
                        tmp = map(float,temp[i+count].split())
                        if homogen > 0:
                            if count < self.load.n_layer+2:
                                tmp[0] = homogen
                        line.extend(tmp)
                        mod.append(line)
                        count +=1
                    except EOFError :
                        break
                
                # --- write AarhusINV mod-file ---
                self.write.modfile(outputfile+str(k)+'_'+str(j)+'.mod',text,mod)                

        print 'done'
#==============================================================================
    def noisemodel2file(self,mean=0,std=1,number=40):
#==============================================================================
        """
        Writes a noise-ini file with a noise model for fwr2tem.
        
        
        Parameters 
        ----------
        mean :              float
            mean of gaussian distribution       
        std :               float
            spread of gaussian distribution (1*STD)        
        number :            int
            length of vector, number of draws from gaussian distribution
        

        """
        
        print 'noise2model              ',
        
        # --- noise model as numpy vector ---
        G = np.random.normal(mean,std,number)
        
        # --- write G vector to file ---
        fid1 = open(self.load.gp_noise_model_file,'w')
        for i in range(number):
            fid1.write('%9.11f \n'%(G[i]))
        fid1.close()
        
    
        print 'done'
        
#==============================================================================
    def process_tem(self,pp_x,pp_y,err=0.8,data_min=0.0):
#==============================================================================
        """  
        New tem data files are created.   
        Tem is autoprocessed.

        Parameters 
        ----------
        pp_x :              int
            python list. Node positions of pilot points in column direction.
        pp_y :              int 
            python list. Node positions of pilot points in row direction.
        err :               float
            Spread of gaussian distribution, 1 STD.
        data_min :          int
            length of vector. Number of draws from gaussian distribution.

        Returns
        -------
        n_tem_data  :       List(int)
            python list. Index holdes number of gates looped over total 
            number of tem-file.
        outputfile  :       dummy.tem
 
        Notes
        -----
        The codes cuts data if the following conditions are meet.
        - data_i <= data_i+1
        - data_i <= 0
        - err <= error_i
        
        See master Thesis by Jacobsen 
        """  

               


        print 'process_tem:               '        

        
        self.n_tem_data = []
        count3 = 0
        for i in pp_y:
                for j in pp_x:
                        fid1 = open(self.load.tem_name_raw+str(i)+'_'+str(j)+'.tem','r')
                                                
                        data     = []
                        error    = []
                        count1 = 0
                        while True:                     
                                try:
                                        if count1 < self.load.n_data_header:
                                                text = fid1.readline()
                                                
                                        if count1 > self.load.n_data_header:
                                                text = fid1.readline().split()
                                                data.extend(  [float(text[1])] )
                                                error.extend( [float(text[2])] )                                                        
                                except:
                                        break
                                count1 += 1                       
                        fid1.close()
        
                        fid1 = open(self.load.tem_name_raw+str(i)+'_'+str(j)+'.tem','r')
                        fid2 = open(self.load.tem_name+str(i)+'_'+str(j)+'.tem','w')
                        
                        
                        count2 = 0
                        for k in range(1,count1-1,1):
                                if k-1 < self.load.n_data_header:
                                        text = fid1.readline()
                                        fid2.write('%s'%(str(text)))
                                if k-1 > self.load.n_data_header:
                                        if err >= error[count2] and data_min <= data[count2] and data[count2] >= data[count2+1]:                                
                                                text = fid1.readline().split()

                                                fid2.write('%s   %s    %s    %s   %s \n'%( text[0],text[1],text[2],text[3],text[4]))
                                                count2 += 1
                                                count3 += 1
                                        else:
                                                text = fid1.readline().split()
                                                fid2.write('%s   %s    %s    %s   %s \n'%( '%'+text[0],text[1],text[2],text[3],text[4]))
                                
                        fid1.close()
                        fid2.close()
                        self.n_tem_data.append(count2)
        
        print 'done'
        return self.n_tem_data  
        
        
#==============================================================================
    def voo2vox(self,fname,outputname,pp_x,pp_y,norms,n_iter,sharp):
#==============================================================================
        """
        Read old AarhusINV voo file. 
        Voo file is a last/optimal run file similar to vox.
               
        
        Parameters
        ----------
        fname :             str
            voo filename from AarhusINV run        
        outname :           str
            vox filename        
        pp\_x :              int
            python list - node positions of pilot points in column direction
        pp\_y :              int 
            python list - node positions of pilot points in row direction       
        n\_iter :            int
            number of iterations in vox file
        
        Returns
        -------
        fname.vox
        
        References
        ----------
        [8] Fidencia, Auken - 2015 (draft)
        
        """
        
        # --- item to find in voo-file ---
        find1 = '%1'
        find2 = '%2'
        find3 = '%3'
        find4 = '%4'
        find5 = '%5'
        find6 = '%6'
        find7 = '%7'
        find8 = '%8'
        
        # --- load voo file ---
        fid1 = open(fname+'.voo','r')
        text = fid1.readlines()
        fid1.close()
        
        # --- open new vox file ---
        fid2 = open(outputname+'.vox','w')
        
        skip = []
        for i in range(len(text)):
            
            if find1 in text[i]:
                skip.extend([i+1])
                fid2.write(text[i])
                
                tmp1    = text[i+1].split()
                tmp1[0] =' '+tmp1[0]

                # --- change number of norms ---
                tmp1[3] =str(norms)
                
                # --- change number of iterations ---
                tmp1[6] =str(n_iter)
                tmp1[0] =' '+tmp1[0]
                
                tex1    = '   '.join(tmp1)
                fid2.write(tex1+'\n')
            

            elif find2 in text[i]:
                skip.extend([i+1])
                fid2.write(text[i])
                fid2.write(text[i+1])
                

            elif find3 in text[i]:
                skip.extend([i+1])
                fid2.write(text[i])
                fid2.write(text[i+1])
                
                tmp1    = text[i+1].split()
                datafiles = int(tmp1[4])
                
            elif find4 in text[i]:
                skip.extend([i+1])
                fid2.write(text[i])
                
                tmp1    = text[i+1].split()
                tmp1[0] =' '+tmp1[0]
                
                # --- change constrain Mesh name ---
                tmp1[1] = outputname+'_CVMESH.vtk'
                
                tex1 = '   '.join(tmp1)
                fid2.write(tex1+'\n')
           
            elif find5 in text[i]:
                for j in range(norms):
                    skip.extend([i+j+1])
                
                fid2.write(text[i])
                if sharp == 0:
                    for j in range(1,norms+1,1):
                        fid2.write(text[i+j])
                        
                elif sharp == 1:
                    fid2.write('                  1                      1                                            L2  \n')
                    fid2.write('                  2                     31            1      0.25 RescaledMinimumSupport \n')
                    fid2.write('                  3                     31            1      0.25 RescaledMinimumSupport \n')
                    fid2.write('                  4                     31          0.5      1.25 RescaledMinimumSupport \n')
                    fid2.write('                  5                      1                                            L2 \n')
                
                   
           
            elif find6 in text[i]:
                skip.extend([i+1])
                fid2.write(text[i])
                fid2.write(text[i+1])
          
            elif find7 in text[i]:
                skip.extend([i+1])
                fid2.write(text[i])
                fid2.write(text[i+1])
                
            elif find8 in text[i]:
                skip.extend([i+1])
                fid2.write(text[i])
                
                for j in range(int(datafiles/2)):
                    skip.extend([i+1+j*2])
                    skip.extend([i+2+j*2])
                    
                    tmp1    = text[i+j*2+1].split()
                    tmp1[0] =' '+tmp1[0]
                    tmp2    = text[i+j*2+2].split() 
                    tmp2[0] =' '+tmp1[0]                
                    
                    for h in range(11):
                        if h != 5:
                            tmp1[h] = int(tmp1[h])
                            tmp2[h] = int(tmp2[h])
                            # --- change normid --
                            if h == 3:
                                tmp1[h] = 5
                                tmp2[h] = 5
                                    
                            
                            
                        
                    delim = '%26.1i %22.1i %22.1i %22.1i %22.1i %29.60s %7.1i %7.1i %7.1i %7.1i %7.1i \n'
                    fid2.write(delim %(tmp1[0],tmp1[1],tmp1[2],tmp1[3],tmp1[4],tmp1[5],tmp1[6],tmp1[7],tmp1[8],tmp1[9],tmp1[10]))
                    fid2.write(delim %(tmp2[0],tmp2[1],tmp2[2],tmp2[3],tmp2[4],tmp2[5],tmp2[6],tmp2[7],tmp2[8],tmp2[9],tmp2[10]))
                # --- break loop after writing datafiles to vox-file ---
                break 

            else:
                if i in skip:
                    continue
                elif 'AarhusInv, version' in text[i]:
                    break
                else:
                    fid2.write(text[i])
            
        fid2.close()
        
        
#==============================================================================
    def voxel_setup(self,pp_x,pp_y,nite1=50,sharp=0,nite2=40):         
#==============================================================================
        """        
        Main module for generating and running a voxel inversion            
        
        
        Parameters
        ----------
        pp_x :              int
            python list - node positions of pilot points in column direction

        pp_y :              int 
            python list - node positions of pilot points in row direction
                            
        nite1 :             int
            number of iterations for smooth inversion
        
        sharp :             int
            run sharp inversion after smooth inversion
            0 --> do no run sharp inversion
            1 --> run sharp inversion
        
        nite2 :             int
            number of iterations of sharp inversion
        
        Returns
        -------
        
        Notes
        -----
        Last change:
        2015-06-03
        
        -name convension is not optimal
        -generation of forward data, how does AarhusINV add noise?
        --> check *.con file 
        --> vertical average, see quasi_3D.....
        
        
        """
        
        print '###############################################################'
        print 'voxel_setup:'
        
        
#        # --- file name for voxel inversion ---
        fname = 'HYTEB_voxel'
        voxel_sufix ='_MVMesh01_Ite00_Complete'
#        
#        
#        shutil.copyfile('tem_example_high_moment.ini','tem_example_high_moment.tem')
#        shutil.copyfile('tem_example_low_moment.ini','tem_example_low_moment.tem')
#
#        # --- write vox-file ---
#        modelmeshname       = fname+'.mod'
#        forwardmeshname     = fname+'.mod'
#        constraintmeshname  = fname+'.mod'
#        self.write.vox(fname,modelmeshname,forwardmeshname,constraintmeshname,
#                       self.load,pp_x,pp_y,itarations=-100,run_type='forward',altitude=False)        
#
#       
#        # --- resisisvity section; res, STD_prior, vertical_con, horizontal_con ---    
#        ver_con = np.round(np.linspace(0.85,0.1,self.load.n_layer-1),decimals=1 )  
        hor_con = np.round(np.linspace(0.1,0.3,self.load.n_layer),decimals=1 )     
#        
#        # --- write voxel mod-file ---
#        self.write.mod_voxel(fname,pp_x,pp_y,ver_con,hor_con,self.load,run_type='dummy',altitude=False)       
#              
#        # --- run aarhusInv to build model ---
#        print ' Build model:             ',
#        self.run_aarhus_inv_voxel(fname)
#        print 'done'
#               
#        # ---------------------------------------------------------------------
#        # --- run AarhusINV voxel setup, set iteration = 0 and use output model 
#        #    from res2vtk to generate forward data ---
#        # --------------------------------------------------------------------
#       
#        # --- rewrite AarhusINV vtk file with "true" model --- 
#        self.res2vtk(fname+voxel_sufix,fname)
#            
#        # --- voxel setup for generating forward data ---
#        modelmeshname       = fname+'.vtk'
#        forwardmeshname     = fname+'_FVMesh01_Ite00.vtk'
#        constraintmeshname  = fname+'.mod'
#        self.write.vox('_1_'+fname,modelmeshname,forwardmeshname,constraintmeshname,
#                       self.load,pp_x,pp_y,itarations=0,run_type='forward',altitude=False)        
#
#        print ' Generate forward data:   ',
#        self.run_aarhus_inv_voxel('_1_'+fname)
#        print 'done'
#        
#        # --- copy data from fwr2tem into subfolder ---
#        path_tem_data = os.getcwd()+'\\TEM_data\\'
#        try:
#            old_dir = os.getcwd()           
#            os.mkdir(path_tem_data)
#
#        except OSError as e:
#            try:
#                shutil.rmtree(path_tem_data)
#                os.mkdir(path_tem_data)
#            except shutil.Error as e:
#                print 'Error in copy2slave:'
#                print('Directory not copied. Error: %s' % e) 
#                print 'folder does not excist!!! '
#
#        count1 = 0
#        for i in pp_y:
#            for j in pp_x:
#                for h in range(2):
#                    count1 += 1
#                    if count1 < 10:
#                        numb = '0000'
#                    elif count1 >=10  and count1 <100:
#                        numb = '000'
#                    elif count1 >=100  and count1 <1000:
#                        numb = '00'
#                    elif count1 >=1000  and  count1 <10000:
#                        numb = '0'
#                    else:
#                        numb = ''
#                    fwr_response = '_1_'+fname+numb+str(count1)+'.fwr'    
#                    if np.mod(count1,2) != 0:
#                        
#                        shutil.copyfile(fwr_response,
#                                        path_tem_data+self.load.tem_name+str(i)+'_'+str(j)+'_low.tem')
#                        shutil.move(fwr_response,path_tem_data+fwr_response)  
#                                      
#                    elif np.mod(count1,2) == 0:
#                        shutil.copyfile(fwr_response,
#                                        path_tem_data+self.load.tem_name+str(i)+'_'+str(j)+'_high.tem')
#                        shutil.move(fwr_response,path_tem_data+fwr_response)
#                    else:
#                        print 'error in voxel_setup: copy fwr2tem'
#       
#       
#        # -----------------------------------------------
#        # --- voxel inversion ---
#        # -----------------------------------------------
#
#        # --- write voxel mod file ---        
#        self.write.mod_voxel('_2_'+fname,pp_x,pp_y,ver_con,hor_con,self.load,run_type='inversion',altitude=False,path_data=path_tem_data)
#        
#        # --- write voxel constain mesh file from forward run--
#        # --- write different horziontal constrains along x and y ---
#        inputfile   = '_1_HYTEB_voxel_CVMesh01_Ite00_Complete.vtk'
#        outputfile  = fname+'CVmesh_diff_xy.vtk'
#        self.cvmesh2cvmesh(inputfile,outputfile,hor_con,factor=0.2,index=3)
                            
#        
#        # --- write vox file --- 
#        modelmeshname       = fname+'_homo'+'.vtk'
#        forwardmeshname     = fname+'_FVMesh01_Ite00.vtk'
#        constraintmeshname  = outputfile
#        self.write.vox('_2_'+fname,modelmeshname,forwardmeshname,constraintmeshname,
#                       self.load,pp_x,pp_y,itarations=nite1,run_type='inversion',altitude=False,path_data=path_tem_data)        
#
#        # --- rewrite AarhusINV vtk file with homogeneus initial model for inversion --- 
#        self.res2vtk(fname, fname+'_homo',homogen=40)        
#
#        #--- run AarhusINV voxel inversion ---
#        print ' Voxel inversion:         ',
#        self.run_aarhus_inv_voxel('_2_'+fname)
#        print 'done'
#        
        
        # ---run AarhusInv voxel SHARP inversion ---
        if sharp == 1:
            
            # --- rewrites smooth inversion results to new inputfile for SHARP... 
            # --- inversion ---
            self.voo2vox('_2_HYTEB_voxel','_3_SHARP_'+fname,pp_x,pp_y,5,nite2,sharp)

#            # rewrites smooth inversion CVMESH file from smooth to SHARP constrains ---
#            inputfile   = '_2_HYTEB_voxel_CVMesh01_Ite0'+str(nite1)+'_Complete.vtk'
#            outputfile  = '_3_SHARP_HYTEB_voxel_CVMESH.vtk'
#            self.cvmesh2cvmesh(inputfile,outputfile,hor_con,factor=1,index=3,index_sharp=4)
        
            # ---run AarhusInv voxel SHARP inversion ---
#            self.run_aarhus_inv_voxel('_3_SHARP_'+fname)
        
        print '###############################################################'
        
#==============================================================================
    def res2vtk(self,fin,fout,homogen=-1):
#==============================================================================
        """
        Load ArhusINV voxel vtk model file.
        substituting resistivity model made from the mod-template file, with a user defind model
        
     
        Parameters
        ----------
        fin :               str
            filename of 3D grid
        
        Returns
        -------
        fout :              str
            output filename of 3D voxel vtk model, definded by the user 
                        
        """
        
        print ' res2vtk:                 ',
        

        fid1    = open(fin+'.vtk','r')
        lines   = fid1.readlines()
        fid1.close()
        
        fid2 = open(fout+'.vtk','w')
        
        # --- load 3D model ---
        if homogen == -1:
            hk  = self.load.model(self.load.model_name,self.load.dtype_model,self.load.nx,self.load.ny,self.load.n_layer)
            ib  = self.load.model(self.load.ibound,self.load.dtype_ibound,self.load.nx,self.load.ny,self.load.n_layer)
            model = np.flipud(hk*ib)
            hk_vec = np.reshape(model,(self.load.nx*self.load.ny*self.load.n_layer))
            del hk, ib
        
        # --- find line in vtk file ---        
        find = 'Res_#'
        count1 = 0
        
        while count1 < len(lines):
            # --- bug in AarhusINV, writes out the dx dy dis wrong, ---
            if count1 == 1:
                tmp1 = lines[count1].split()
                tmp1[-1] = str(self.load.dx)
                tmp1[-2] = str(self.load.dy)
                tmp1 = map(float,tmp1)
                fid2.write('%15.0i %5.0i %5.0i %15.4f %5.7e %5.7e %5.7e %5.1f %5.1f \n' % (tmp1[0],tmp1[1],tmp1[2],tmp1[3],tmp1[4],tmp1[5],tmp1[6],tmp1[7],tmp1[8]))

            # --- find and write res into vtk file
            elif find in lines[count1]:
                fid2.write(lines[count1])
                tmp     = lines[count1].split() 
                n_par   = int(tmp[-2])                
                for i in range(1,n_par+1,1):
                    tmp     = lines[count1+i].split()
                    
                    
                    if homogen == -1:
                        # --- when hk is impermable(0.0) --> res = res_clay ---
                        if hk_vec[i-1] == 0.0:
                            res = self.load.res_clay
                            tmp[0] = res
    
                        # --- ---                    
                        elif hk_vec[i-1] < 0.0:
                            res = -1*hk_vec[i-1]
                            tmp[0]  = (res/self.load.beta_1)**(1/self.load.beta_2)
                        
                        # --- transform hydraulic conducitivty to resistivity by a petrophysical relationship ---
                        else:
                            res = hk_vec[i-1]
                            tmp[0]  = (res/self.load.beta_1)**(1/self.load.beta_2)                        
                    
                    # --- homogen model ---
                    elif homogen > 0:
                        tmp[0] = homogen
                    # --- write new resistivity values into "find" section ---
                    tmp = map(float,tmp)
                    fid2.write('%15.8e %15.8e %5.2e %5.2e %5.2e %5.1i %5.1i %5.2e %5.1i %5.1i %5.1i \n' % (tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[5],tmp[6],tmp[7],tmp[8],tmp[9],tmp[10]))
                
                # --- break after npar iterations ---
                break
            
            # --- write first part of the vtk file ---
            else:
                fid2.write(lines[count1])
            count1 += 1
            
        fid2.close()
        
        print 'done'
        
#==============================================================================    
    def run_aarhus_inv(self,inputfile,pp_x,pp_y,config='em1dinv.con'):     
#==============================================================================
        """
        Run AahusINV

        Parameters
        ----------
        inputfile :         str
            mod-filename
        
        pp_x :              int
            python list - node positions of pilot points in column direction

        pp_y :              int 
            python list - node positions of pilot points in row direction
        
        config :            str            
            name of AarhusINV configuration file (*.con)
        

        """
        
        print 'run_aarhus_inv:          ',        
        
        for ny_pos in pp_y: 
            for nx_pos in pp_x: 
                program = self.load.geophys_software+'.exe'
                input1  =  inputfile+str(ny_pos)+'_'+str(nx_pos)+'.mod'
                input2  = config
                dos_tex = [program,input1,input2]
                outputtext= subprocess.call(dos_tex)
        
        print 'done'
        


#==============================================================================
    def run_aarhus_inv_voxel(self,fname):
#==============================================================================
        """
        Run voxel model
        
        Parameters
        ----------
        fname :             str
            vox-filename        
        
        Returns
        -------
        HYTEB_voxel.err
        HYTEB_voxel.dat
                
        Notes
        -----
        
        References
        ----------
        
        
        """
        program = self.load.geophys_software.split('.')
        input1  =  fname+'.vox'
        input2  =  program[0]+'.con'
        dos_tex = [program[0]+'.'+program[1],input1,input2]
        
#        log_file = 'HYTEB_voxel_log.err'
#        file_out = 'HYTEB_voxel_log.dat'
#        fid1 = open(file_out,'w+')
#        fid2 = open(log_file,'w+')
#        subprocess.call(dos_tex,stdout= fid1,stderr =fid2)
        subprocess.call(dos_tex)
#        fid1.close()
#        fid2.close()  
#           

#==============================================================================
    def tem_mod_quasi_3D(self,pp_x,pp_y,res_ver_con=[],dim='q3D'):     
#==============================================================================
        """
        Writes mod-files for AarhusInv. The 1D models are sampled as quasi 3D 
        
        Parameters
        ----------
        pp"_"x :              int
            python list  node positions of pilot points in column direction
        pp"_"y :              int 
            python list  node positions of pilot points in row direction
        res"_"ver"_"con :       array_like(float)
            python list type. Vertical constarins on resistivity values
        
        """
        
        
        print 'tem_mod_quasi_3D:' 
        
        # --- vertical constrains ---
        if res_ver_con  == []:
            res_ver_STD =  -1*np.ones(self.load.n_layer)
        else:            
            res_ver_STD = res_ver_con[0]*np.ones(self.load.n_layer)
            
        # --- number of sounding locations ---
        n = len(pp_x)*len(pp_y)
        
        # --- copy TEM template file, changed sufix from *.ini to *.tem ---
        tmp = self.load.data_template.split('.')
        try:           
            shutil.copyfile(tmp[0]+'.'+tmp[1],tmp[0]+'.tem')
        except:
            if IOError:
                print '     geophysical data template files is missing'
            
        # --- write quasi 3D mod-files for n sound locations ---
        print '     number of TEM-soundings: ',n
        count1 = 1
        for ny_pos in pp_y:
            for nx_pos in pp_x:
                print '         ',count1,'/',n
                self.write.mod_tem_q3D(nx_pos,ny_pos,dim,self.load,res_ver_con=res_ver_STD)
                count1 += 1    
        print '         done'
