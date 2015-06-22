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

@author: durack1
"""
import EzTemplate,gc,glob,os,time,re,resource,vcs
import cdms2 as cdm
import numpy as np

#%% Create input file list
newList    = sorted(glob.glob('/work/durack1/Shared/150219_AMIPForcingData/360x180_150618/*.nc'))
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
del(newList,yrTest,x,filePath); gc.collect()

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
del(oldList,yrTest,x,filePath); gc.collect()

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
        varName = 'tos'
        levs    = list(np.arange(270,310,2.5)) ; # TOS
        split   = 1
    else:
        varName = var
        levs    = list(np.arange(-5,115,10)) ; # SIC
        split   = 0

    #%% Configure VCS - and set iso levels per variable
    x           = vcs.init()
    x.setcolormap("bl_to_darkred")
    iso         = x.createisofill()
    #iso         = x.createboxfill()
    #iso.boxfill_type="custom"
    iso         = x.createisofill()
    iso.levels  = levs
    iso.ext_1   = True
    iso.ext_2   = True
    cols        = vcs.getcolors(iso.levels,split=split)
    iso.fillareacolors = cols

    iso2        = x.createisofill()
    levs        = list(np.arange(-.2,.21,.025))
    levs[8]     = 0 ; # Fix middle point
    iso2.levels = levs
    iso2.ext_1  = True
    iso2.ext_2  = True
    cols        = vcs.getcolors(iso2.levels)
    iso2.fillareacolors = cols

    leg         = x.createtextorientation()
    leg.halign  = "left"
    leg.height  = 8

    tmpl = x.createtemplate()
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

    title           = x.createtext()
    title.height    = 14
    title.halign    = "center"
    title.x         = [.5]
    title.y         = [.975]

    t1 = Ez.get(legend="local")
    t2 = Ez.get(legend="local")
    t3 = Ez.get(legend="local")

    #%% Setup canvas options and plot
    #bg = True ; # For 1 yr uses ~2.4GB
    # [durack1@oceanonly 150219_AMIPForcingData]$ xeyes ; # Should display in foreground
    # [durack1@oceanonly 150219_AMIPForcingData]$ Xvfb :2 -screen 0 1600x1200x16
    # [durack1@oceanonly 150219_AMIPForcingData]$ bg
    # [durack1@oceanonly 150219_AMIPForcingData]$ setenv DISPLAY :2
    # [durack1@oceanonly 150219_AMIPForcingData]$ xeyes ; # Should display in background
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
        for count,y in enumerate(range(1870,y2)):
            f1  = cdm.open(newList[count]) ; # New files
            f2  = cdm.open(oldList[count]) ; # Downloadable files ~2012
            #print newList[count]
            #print oldList[count]
            #continue
            for m in range(12):
                startTime = time.time()
                printStr = 'processing: %i-%.2i' % (y,m+1)
                s1 = f1(varNameRead,slice(m,m+1))
                s2 = f2(varNameRead,slice(m,m+1))
                #diff = s2+np.random.random(s1.shape)*20.-10. -s1
                diff = s2-s1
                x.plot(title,bg=bg)
                x.plot(s1,t1,iso,bg=bg); #,ratio="autot"); #,vtk_backend_grid=g)
                x.plot(diff,t2,iso2,bg=bg); #,ratio="autot") ; #,vtk_backend_grid=g)
                x.plot(s2,t3,iso,bg=bg); #,ratio="autot") ; #,vtk_backend_grid=g)
                fnm = '%i-%.2i.png' % (y,m+1)
                fileName = os.path.join(outPath,'pngs',varNameRead,fnm)
                outFiles.append(fileName)
                x.png(fileName)
                x.clear()
                endTime = time.time()
                timeStr = 'Time: %06.3f secs;' % (endTime-startTime)
                memStr  = 'Max mem: %05.3f GB' % (np.float32(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)/1.e6)
                counterStr = '%05d' % counter
                print counterStr,printStr,varName.ljust(6),BC,timeStr,memStr
                counter = counter+1
            f1.close()
            f2.close()
            gc.collect() ; # Attempt to force a memory flush
        outMP4File = ''.join(['AMIPBCS_newVsOld_',varNameRead,'.mp4'])
        print 'Processing: ',outMP4File
        x.ffmpeg(os.path.join(outPath,outMP4File),outFiles,rate=5,bitrate=2048); #,options=u'-r 2') ; # Rate is frame per second - 1/2s per month
    x.clear()
    outFiles = [] ; # Reset for obs vs bcs

#%%
'''
# Chuck tips
x=vcs.init()
iso = x.createisofill()
iso.levels = [0,1,2,3,4]
cols = vcs.getcolors(iso.levels)
iso.fillareacolors = cols
x.plot(s,iso,bg=1) # add gms in here
x.animation.create()
x.animation.save('paul.mp4')

# Karl
Make movies
diff = [New(t,y,x)-Old(t,y,x)]
max = abs(diff(y,x))
rms = sum(diff(t,y,x)^2)/n(t)
'''
