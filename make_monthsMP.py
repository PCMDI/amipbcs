#!/bin/env python
# -*- coding: utf-8 -*-
"""
This code calls make_newVsOldDiffsMP.py

Will require edits for paths and input temporal periods
"""

import glob,os,shlex,subprocess,sys,time,vcs
sys.path.append('/export/durack1/git/durolib/lib/')
from durolib import mkDirNoOSErr
from multiprocessing import Pool

homeDir = '/work/durack1/Shared/150219_AMIPForcingData'
dataVer = 'v1.1.2'
homeDir = os.path.join(homeDir,'pngs','_'.join(['pngs',dataVer]))

if not os.path.exists(homeDir):
    mkDirNoOSErr(homeDir)

#%% Cleanup output dirs
for root, dirs, files in os.walk(homeDir, topdown=False):
    for name in files:
        print os.path.join(root,name)
        os.remove(os.path.join(root, name))

#%% Create loop
cmds = []
for i in sorted(range(1870,2016)): # Original data only extends to 2012
    cmd = "python make_newVsOldDiffsMP.py %i" % i
    cmds.append(cmd)
    cmds.sort()

print 'start:',cmds[0].split(' ')[-1]
print 'end:  ',cmds[-1].split(' ')[-1]

#%% Def function
def make_year(cmd):
    print "running",cmd
    p = subprocess.Popen(shlex.split(cmd))
    p.wait()
    return p

#%% Run
if __name__ == "__main__":
    # Specify step count
    poolStep = 20
    # Run using pool step
    p = Pool(poolStep)
    time.sleep(1)
    p.map(make_year,cmds)

#%% Make movies
for fileList in ['sic','sicbcs','tos','tosbcs']:
    # Create outFileName
    outMP4File = ''.join(['AMIPBCS_newVsOld_',fileList,'.mp4'])
    print 'Processing: ',outMP4File
    # Get files
    outFiles = sorted(glob.glob(os.path.join(homeDir,fileList,'*.png')))
    x = vcs.init()
    x.ffmpeg(os.path.join(homeDir,outMP4File),outFiles,rate=5,bitrate=2048); #,options=u'-r 2') ; # Rate is frame per second - 1/2s per month
    x.close()

#%%
print "Done"
