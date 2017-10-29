from LRIS.resample import resample
from LRIS.LRIStools import *
from LRIS.XSolve import *

indir = '/Users/eduardrusu/Desktop/Reduction/1206m1'
pref = 'r170322_'
flat = 31
arc = 32
outname = '1206m1_red'
slitID(indir,pref,[flat,arc,28],outname,side='bottom')
oldName = None
sf = True # Update to determine whether you want to see the fit for every exposure
for img in [28,29,33]:
    newName = '%s_%2d'%(outname,img)
    XSolve(outname,newName,indir,pref,[flat,arc,img])
    SlitCross(newName,showFit=sf)
    WaveSolve(newName,oldName,showFit=sf)
    resample(newName,nobgsub=False,clobber=True)
    oldName = newName
    sf = True # Don't show the fit anymore if True...
