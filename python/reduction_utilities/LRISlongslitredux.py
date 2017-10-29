from LRIS.resample import resample
import pyfits,numpy
from LRIS.LRIStools import *
from LRIS.XSolve import *

indir = 'data'
pref = 'r110501_'
flat = 41 
arc = 36
outname = 'red_top'
slitID(indir,pref,[flat,arc,67],outname,side='top',slits=[[50,100]])
oldName = None
sf = True
for img in [67,68,69]:
    newName = '%s_%2d'%(outname,img)
    XSolve(outname,newName,indir,pref,[flat,arc,img])
    SlitCross(newName)
    WaveSolve(newName,oldName,showFit=sf)
    resample(newName,nobgsub=True,clobber=True)
    oldName = newName

