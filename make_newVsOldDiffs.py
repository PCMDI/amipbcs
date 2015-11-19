#!/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu May 21 18:30:38 2015

Paul J. Durack 21st May 2015

This script generates a movie of monthly mean differences

PJD 21 May 2015     - Started
PJD 26 May 2015     - Updated with old (2012) AMIP data for comparisons
PJD  3 Jun 2015     - Finalised with diff plots
PJD  4 Jun 2015     - Moved x.close() to prevent issues
PJD 16 Jun 2015     - Updated to deal with latest files and BCs rather than rewritten obs
PJD 17 Jun 2015     - Grabbed tweaks from doutriaux1 to reduce grid reproduction
PJD 17 Jun 2015     - Updated to deal with 4 variables obs and bcs
PJD 18 Jun 2015     - Added time and memory statements
PJD 18 Jun 2015     - Updated 'new' data to 360x180_150618 path
PJD 24 Jun 2015     - Updated levels for 'sic' and corrected outFiles indentation to fix over-runs in mp4 files
PJD 15 Jul 2015     - Added UV-CDAT version attribution to be logged
PJD 16 Jul 2015     - Added delFudge variable
PJD 18 Nov 2015     - Added pyObj variable

@author: durack1
"""
import cdat_info,EzTemplate,gc,glob,os,time,re,resource,vcs
import cdms2 as cdm
import numpy as np

#%% Turn on purging of VCS objects?
delFudge = True

#%% Define functions
def initVCS(x,levs1,levs2,split):
    x.setcolormap("bl_to_darkred")
    iso1                = x.createisofill()
    #iso         = x.createboxfill()
    #iso.boxfill_type="custom"
    iso1.levels         = levs1
    iso1.ext_1          = True
    iso1.ext_2          = True
    cols                = vcs.getcolors(iso1.levels,split=split)
    iso1.fillareacolors = cols

    iso2                = x.createisofill()
    levs                = levs2
    iso2.levels         = levs
    iso2.ext_1          = True
    iso2.ext_2          = True
    cols                = vcs.getcolors(iso2.levels,split=1)
    iso2.fillareacolors = cols

    leg                 = x.createtextorientation()
    leg.halign          = "left"
    leg.height          = 8

    tmpl                = x.createtemplate()
    tmpl.blank()
    tmpl.scalefont(.9)

    tmpl.legend.textorientation = leg
    for a in ["data","legend","box1","xlabel1","xtic1","ylabel1","ytic1"]:
        setattr(getattr(tmpl,a),"priority",1)

    Ez                  = EzTemplate.Multi(rows=3,columns=1,x=x,template=tmpl)
    Ez.legend.direction ='vertical'
    Ez.margins.left     =.05
    Ez.margins.right    =.05
    Ez.margins.top      =.05
    Ez.margins.bottom   =.05

    title               = x.createtext()
    title.height        = 14
    title.halign        = "center"
    title.x             = [.5]
    title.y             = [.975]

    t1                  = Ez.get(legend="local")
    t2                  = Ez.get(legend="local")
    t3                  = Ez.get(legend="local")
    
    return iso1,iso2,title,t1,t2,t3


#%% Create input file list
newList    = sorted(glob.glob('/work/durack1/Shared/150219_AMIPForcingData/360x180/*.nc'))
sicList,sicbcList,tosList,tosbcList = [[] for _ in range(4)]
yrTest      = re.compile('[0-9]{4}.nc')
for x,filePath in enumerate(newList):
    if 'amip' in filePath.split('/')[-1].split('_')[0] and yrTest.search(filePath.split('/')[-1].split('_')[-1]):
        if filePath.split('/')[-1].split('_')[1] == 'sic':
            if 'amipobs' in filePath.split('/')[-1].split('_')[0]:
                sicList.append(filePath)
            else:
                sicbcList.append(filePath)
        if filePath.split('/')[-1].split('_')[1] == 'sst':
            if 'amipobs' in filePath.split('/')[-1].split('_')[0]:
                tosList.append(filePath)
            else:
                tosbcList.append(filePath)
del(newList,yrTest); gc.collect()

oldList    = sorted(glob.glob('/work/durack1/Shared/150219_AMIPForcingData/150616_AMIPweb/www-pcmdi.llnl.gov/360x180/*.nc'))
sicList2,sicbcList2,tosList2,tosbcList2 = [[] for _ in range(4)]
yrTest      = re.compile('[0-9]{4}.nc')
for x,filePath in enumerate(oldList):
    if 'amip' in filePath.split('/')[-1].split('_')[0] and yrTest.search(filePath.split('/')[-1].split('_')[-1]):
        if filePath.split('/')[-1].split('_')[1] == 'sic':
            if 'amipobs' in filePath.split('/')[-1].split('_')[0]:
                sicList2.append(filePath)
            else:
                sicbcList2.append(filePath)
        if filePath.split('/')[-1].split('_')[1] == 'sst':
            if 'amipobs' in filePath.split('/')[-1].split('_')[0]:
                tosList2.append(filePath)
            else:
                tosbcList2.append(filePath)
del(oldList,yrTest); gc.collect()

outPath = '/work/durack1/Shared/150219_AMIPForcingData'

#%% Purge existing files
for root, dirs, files in os.walk('./pngs', topdown=False):
    for name in files:
        #print os.path.join(root,name)
        os.remove(os.path.join(root, name))

#%% Loop through vars and files
counter = 1
for var in ['sic','sst']:
    if var == 'sst':
        varName     = 'tos'
        levs1       = list(np.arange(270,310,2.5)) ; # TOS
        levs2       = list(np.arange(-.2,.21,.025))
        levs2[8]    = 0. ; # Fix middle point
        split       = 1
    else:
        varName     = var
        levs1       = list(np.arange(-5,115,10)) ; # SIC
        levs2       = list(np.arange(-5,5.5,.5))
        split       = 0

    #%% Setup canvas options and plot
    #bg = True ; # For 1 yr uses ~2.4GB
    # [durack1@oceanonly 150219_AMIPForcingData]$ xeyes ; # Should display in foreground
    # [durack1@oceanonly 150219_AMIPForcingData]$ Xvfb :2 -screen 0 1600x1200x16
    # [durack1@oceanonly 150219_AMIPForcingData]$ bg ; # Ctrl-z called to send to background
    # [durack1@oceanonly 150219_AMIPForcingData]$ setenv DISPLAY :2
    # [durack1@oceanonly 150219_AMIPForcingData]$ xeyes ; # Should display in background
    # [durack1@oceanonly 150219_AMIPForcingData]$ jobs
    # [1]  + Running                       spyder
    # [2]  - Running                       Xvfb :2 -screen 0 1600x1200x16
    bg = False ; # For 1 yr uses ~260MB
    y2 = 2013 ; #1871; #2013
    outFiles = []

    #%% Set obs vs bcs
    for data in ['obs','bcs']:
        if data == 'bcs':
            BC = 'bc'
            varNameRead = ''.join([varName,'bcs'])
        else:
            BC = ''
            varNameRead = varName
        newList = eval(''.join([varName,BC,'List']))
        oldList = eval(''.join([varName,BC,'List2']))
        x = vcs.init()
        basic_tt = vcs.elements["texttable"].keys()
        basic_to = vcs.elements["textorientation"].keys()
        basic_tc = vcs.elements["textcombined"].keys()
        for count,y in enumerate(range(1870,y2)):
            f1  = cdm.open(newList[count]) ; # New files
            f2  = cdm.open(oldList[count]) ; # Downloadable files ~2012
            for m in range(12):
                startTime                   = time.time()
                printStr                    = 'processing: %i-%.2i' % (y,m+1)
                s1                          = f1(varNameRead,slice(m,m+1))
                s2                          = f2(varNameRead,slice(m,m+1))
                diff                        = s2-s1
                iso1,iso2,title,t1,t2,t3    = initVCS(x,levs1,levs2,split)
                title.string                = '%i-%.2i' % (y,m+1)
                x.plot(title,bg=bg)
                x.plot(s1,t1,iso1,bg=bg); #,ratio="autot"); #,vtk_backend_grid=g)
                x.plot(diff,t2,iso2,bg=bg); #,ratio="autot") ; #,vtk_backend_grid=g)
                x.plot(s2,t3,iso1,bg=bg); #,ratio="autot") ; #,vtk_backend_grid=g)
                fnm                         = '%i-%.2i.png' % (y,m+1)
                fileName                    = os.path.join(outPath,'pngs',varNameRead,fnm)
                outFiles.append(fileName)
                x.png(fileName)
                x.clear()
                for k in vcs.elements.keys():
                    pass
                    #print k, len(vcs.elements[k].keys())
                if delFudge:
                    del(vcs.elements["isofill"][iso1.name])
                    del(vcs.elements["isofill"][iso2.name])
                    del(vcs.elements["textcombined"][title.name])                
                    del(vcs.elements["template"][t1.name])
                    del(vcs.elements["template"][t2.name])
                    del(vcs.elements["template"][t3.name])
                    for nm in vcs.elements["texttable"].keys():
                        if not nm in basic_tt:
                            del(vcs.elements["texttable"][nm])
                    for nm in vcs.elements["textorientation"].keys():
                        if not nm in basic_to:
                            del(vcs.elements["textorientation"][nm])
                    for nm in vcs.elements["textcombined"].keys():
                        if not nm in basic_tc:
                            del(vcs.elements["textcombined"][nm])
                endTime                 = time.time()
                timeStr                 = 'Time: %06.3f secs;' % (endTime-startTime)
                memStr                  = 'Max mem: %05.3f GB' % (np.float32(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)/1.e6)
                counterStr              = '%05d' % counter
                pyObj                   = 'PyObj#: %07d;' % (len(gc.get_objects()))
                if counter == 1:
                    print 'UV-CDAT version: ', cdat_info.get_version()
                    print 'UV-CDAT prefix: ', cdat_info.get_prefix()
                    print 'delFudge: ',delFudge
                print counterStr,printStr,varName.ljust(6),BC,timeStr,memStr,pyObj
                #del() ; # Do a cleanup
                counter                 = counter+1
            f1.close()
            f2.close()
            gc.collect() ; # Attempt to force a memory flush
        x.close()
        outMP4File = ''.join(['AMIPBCS_newVsOld_',varNameRead,'.mp4'])
        print 'Processing: ',outMP4File
        x = vcs.init()
        x.ffmpeg(os.path.join(outPath,outMP4File),outFiles,rate=5,bitrate=2048); #,options=u'-r 2') ; # Rate is frame per second - 1/2s per month
        #x.animation.create()
        #x.animation.save('demo.mp4')
        x.close()
        outFiles = [] ; # Reset for obs vs bcs

#%%
'''
# Karl
Make movies
diff = [New(t,y,x)-Old(t,y,x)]
max = abs(diff(y,x))
rms = sum(diff(t,y,x)^2)/n(t)
'''
