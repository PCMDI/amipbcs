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
PJD  6 Jun 2016     - Updated for new input data format
PJD 27 Jun 2016     - Updated to run against vtk bug fixes (uvcdatNightly)
PJD 28 Jun 2016     - Updated to write mp4's to png subdir, test with bg=False
PJD 28 Jun 2016     - Corrected pngs path for mp4s
PJD 30 Jun 2016     - Added donotstoredisplay argument
PJD 22 Aug 2017     - Added code tweaks following VCS changes #235
PJD 28 Aug 2017     - Further tweaks following info in https://github.com/UV-CDAT/vcs/pull/237#issuecomment-325354235
PJD 30 Oct 2017     - Updated to bring changes across from make_newVsOldDiffsMP.py
PJD 30 Oct 2017     - Corrected glob call to only return a single version
PJD 31 Oct 2017     - Updated TOS plot increments 270 to 310 (K) to -2.5 to 37.5 (degC); Updated ver info
PJD 31 Oct 2017     - split =1 needs to change to 0 to allow even colour splitting for unbalanced temp range (-2.5 to 35)
PJD 31 Oct 2017     - Explicitly initialize canvas size vcs.init(bg=True,geometry=())

@author: durack1
"""
import cdat_info,EzTemplate,gc,glob,os,time,resource,sys,vcs
import cdms2 as cdm
import numpy as np
sys.path.append('/export/durack1/git/durolib/lib/')
from durolib import mkDirNoOSErr
from string import replace

#%% Turn on purging of VCS objects?
delFudge = True
outPathVer = 'pngs_v1.1.3'
ver = 'v20171031' ; # Update for each run
verOld = 'v20170419' ; # Update for each run

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

    oldOrientation = tmpl.legend.textorientation
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

    del(Ez) ; # Purge EzTemplate object
    vcs.removeobject(vcs.elements['textorientation'][oldOrientation])

    return iso1,iso2,title,t1,t2,t3,tmpl


#%% Create input file list
newList    = sorted(glob.glob(''.join(['/work/durack1/Shared/150219_AMIPForcingData/CMIP6/input4MIPs/PCMDI/SSTsAndSeaIce/CMIP/mon/*/PCMDI-AMIP-1-1-3/*/gn/',ver,'/*.nc'])))
for x,filePath in enumerate(newList):
    if 'siconc' in filePath.split('/')[13]:
        if 'bcs' in filePath.split('/')[13]:
            sicbcList = filePath
        else:
            sicList = filePath
    if 'tos' in filePath.split('/')[13]:
        if 'bcs' in filePath.split('/')[13]:
            tosbcList = filePath
        else:
            tosList = filePath
del(filePath,newList,x); gc.collect()

oldList    = sorted(glob.glob(''.join(['/work/durack1/Shared/150219_AMIPForcingData/CMIP6/input4MIPs/PCMDI/SSTsAndSeaIce/CMIP/mon/*/PCMDI-AMIP-1-1-2/*/gn/',verOld,'/*.nc'])))
for x,filePath in enumerate(oldList):
    if 'siconc' in filePath.split('/')[13]:
        if 'bcs' in filePath.split('/')[13]:
            sicbcList2 = filePath
        else:
            sicList2 = filePath
    if 'tos' in filePath.split('/')[13]:
        if 'bcs' in filePath.split('/')[13]:
            tosbcList2 = filePath
        else:
            tosList2 = filePath
del(filePath,oldList,x); gc.collect()

outPath = '/work/durack1/Shared/150219_AMIPForcingData'

#%% Purge existing files - lock to version number
for root, dirs, files in os.walk(os.path.join('./pngs',outPathVer), topdown=False):
    for name in files:
        #print os.path.join(root,name)
        os.remove(os.path.join(root, name))

#%% Loop through vars and files
counter = 1
for var in ['sic','tos']:
    if var == 'tos':
        varName     = 'tos'
        levs1       = list(np.arange(-2.5,37.5,2.5)) ; # TOS
        levs2       = list(np.arange(-.2,.21,.025))
        levs2[8]    = 0. ; # Fix middle point
        split       = 0 ; # 0 for -2.5 to 35 scale
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
    bg = True ; # For 1 yr uses ~260MB
    #delFudge = False ; #Turn on purging of VCS objects?
    donotstoredisplay = True ; # Fix from fries2
    y2 = 2017 ; #1871; #2013
    outFiles = []

    #%% Set obs vs bcs
    for data in ['obs','bcs']:
        if data == 'bcs':
            BC = 'bc'
            varNameRead = ''.join([varName,'bcs'])
        else:
            BC = ''
            varNameRead = varName
        # Fix new naming convention
        if 'sic' in varNameRead:
            varNameNewRead = replace(varNameRead,'sic','siconc')
            inflationFactor = 1. #1e2
            unitFactor = 0.
        else:
            varNameNewRead = varNameRead
            inflationFactor = 1.
            unitFactor = 273.16 ; # Fudge errored <1.1.2 and earlier
            unitFactor = 1.
        newList = eval(''.join([varName,BC,'List']))
        oldList = eval(''.join([varName,BC,'List2']))
        x = vcs.init(bg=True,geometry=(1200,1560)) ; # Add bg and geometry
        basic_tt = vcs.elements["texttable"].keys()
        basic_to = vcs.elements["textorientation"].keys()
        basic_tc = vcs.elements["textcombined"].keys()
        # Open new input file
        monthCount = 0
        f1  = cdm.open(newList) ; # New files
        s1  = f1(varNameNewRead,time=('1870',str(y2)))
        f2  = cdm.open(oldList) ; # New files
        s2  = f1(varNameNewRead,time=('1870',str(y2)))
        for count,y in enumerate(range(1870,y2)):
            for m in range(12):
                startTime                   = time.time()
                printStr                    = 'processing: %i-%.2i' % (y,m+1)
                #s1                          = f1(varNameNewRead,slice(monthCount,monthCount+1))
                #s1                          = s1*inflationFactor ; # Correct siconc variables for unit difference
                s1s                         = s1[m:m+1,]
                s1s                         = s1s*inflationFactor ; # Correct siconc variables for unit difference
                s2s                         = s2[m:m+1,]
                s2s                         = s2s*inflationFactor ; # Correct siconc variables for unit difference
                # Test times
                #print 'new:',varNameNewRead.ljust(9),s1s.getTime().asComponentTime()
                #print 'old:',varNameRead.ljust(9),s2s.getTime().asComponentTime()
                diff                        = s2s-s1s
                iso1,iso2,title,t1,t2,t3,tmpl = initVCS(x,levs1,levs2,split)
                title.string                = '%i-%.2i' % (y,m+1)
                x.plot(title,bg=bg)
                x.plot(s1s,t1,iso1,bg=bg); #,ratio="autot"); #,vtk_backend_grid=g)
                x.plot(diff,t2,iso2,bg=bg); #,ratio="autot") ; #,vtk_backend_grid=g)
                x.plot(s2s,t3,iso1,bg=bg); #,ratio="autot") ; #,vtk_backend_grid=g)
                fnm                         = '%i-%.2i.png' % (y,m+1)
                fileName                    = os.path.join(outPath,'pngs',outPathVer,varNameRead,fnm)
                # Create directory tree - version
                if not os.path.exists(os.path.join(outPath,'pngs',outPathVer)):
                    mkDirNoOSErr(os.path.join(outPath,'pngs',outPathVer))
                outFiles.append(fileName)
                # Create directory tree - variable
                if not os.path.exists(os.path.join(outPath,'pngs',outPathVer,varNameRead)):
                    mkDirNoOSErr(os.path.join(outPath,'pngs',outPathVer,varNameRead))
                # Check file exists
                if os.path.exists(fileName):
                    #print "** File exists.. removing **"
                    os.remove(fileName)
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
                #vcs.removeobject(iso1) ; # Error thrown here and below v2.12 triggered this
                #vcs.removeobject(iso2)
                #vcs.removeobject(title)
                #vcs.removeobject(t1)
                #vcs.removeobject(t2)
                #vcs.removeobject(t3)
                vcs.removeobject(tmpl)
                endTime                 = time.time()
                timeStr                 = 'Time: %06.3f secs;' % (endTime-startTime)
                memStr                  = 'Max mem: %05.3f GB' % (np.float32(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)/1.e6)
                counterStr              = '%05d' % counter
                pyObj                   = 'PyObj#: %07d;' % (len(gc.get_objects()))
                if counter == 1:
                    print 'UV-CDAT version:'.ljust(21),cdat_info.get_version()
                    print 'UV-CDAT prefix:'.ljust(21),cdat_info.get_prefix()
                    print 'delFudge:'.ljust(21),delFudge
                    print 'Background graphics:'.ljust(21),bg
                    print 'donotstoredisplay:'.ljust(21),donotstoredisplay
                print counterStr,printStr,varName.ljust(6),BC,timeStr,memStr,pyObj
                #del() ; # Do a cleanup
                counter                 = counter+1
            gc.collect() ; # Attempt to force a memory flush
        x.backend.renWin = None ; # @danlipsa fix UV-CDAT/vcs#237
        x.close()
        f1.close()
        f2.close()
        outMP4File = os.path.join('pngs',''.join(['AMIPBCS_newVsOld_',varNameRead,'.mp4']))
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