#!/usr/bin/env python 
import os,sys
from math import *
from periodica import structure
from scipy import interpolate
from mysub import inputs
import numpy as np
from glob import glob


if ['soo','soo-desktop','cam1-office','carlosj-desktop','Han'].count(os.uname()[1]):vaspexe='vasp '
elif os.uname()[1].split(".")[-1]=="scsystem":vaspexe='srun /home/cam1/vasp/parallel_vasp.4.6/vasp '
else: print "I do not recognize this machine so I dont know which vasp to use...";sys.exit()

oo=inputs("spline=str relax=1 tol=0.25 sinc=5 maxiter=180 crit=1 kmesh=20 krelax=15 kstatic=43 pathb=100 mesh=0 rguess=0 rot=''")
oo.update(sys.argv[1:])
oo.report()


def generate_kmesh(out='KPOINTS'):
  nx=oo.kmesh; tempe=[]
  pos=structure(open('POSCAR').read()) 
  tt=pos.rvec.copy()/(2*np.pi)#/pos.rvec[0,0]
  p1=tt[0][0]/2.; m=-tt[1][0]/tt[1][1]; b=tt[1][1]/2.-m*tt[1][0]/2.
  #print 'hi',np.dot(tt[1],(0.5,0.5*m,0))
  densx=p1/float(nx)
  for i in range(nx+1):
    #densy=(b+m*i*densx)/int((b+m*i*densx)/densx)
    densy=densx
    incy=0
    while 1: 
      #print "%s %s"%tuple(np.dot((i*densx,incy*densy,0),np.linalg.inv(tt))[:2]) 
      tempe.append("%.8f %.8f 0 1"%(i*densx,incy*densy))
      if incy*densy>m*i*densx+b:break 
      incy+=1
  OUT=open(out,'w');OUT.write("from graphene_curve.py\n%s\nCartesian\n"%len(tempe)+"\n".join(tempe));OUT.close()
  OUT=open('../kinput','w');OUT.write("0 0 0\n%s 0 0\n%s %s 0\n0 %s 0   "%(p1,p1,m*p1+b,b))

def generate_highk():
  pos=structure(open('POSCAR').read()) 
  tt=pos.rvec.copy()/(2*np.pi)
  p1=tt[0][0]/2.; m=-tt[1][0]/tt[1][1]; b=tt[1][1]/2.-m*tt[1][0]/2.
  OUT=open('kinput','w');OUT.write("0 0 0\n%s 0 0\n%s %s 0\n0 %s 0\n0 0 0   "%(p1,p1,m*p1+b,b));OUT.close()
  os.system("kgenerator.pl %s;mv kpoints KPOINTS;sed -i 's/Reciprocal/Cartesian/' KPOINTS"%oo.pathb)

def guess(x):
    theta = (x/180.)*pi
    if dir(oo).count('spline'):
      if os.path.exists(oo.spline):
        try: temp=[[float(k) for k in s.split()] for s in open(oo.spline)];x=[s[0] for s in temp];y=[s[1] for s in temp]
        except: print "error reading spline file %s"%oo.spline;sys.exit()
    else:
      x=[float(s) for s in "0.03 1.04 6.96 9.24 13.50 25.03 33.53 46.61 58.05 69.81 77.60 81.15 83.09 83.91 85.08 86.10 87.17 88.93 89.90 91.16 92.08 92.69 93.79 95.04 95.87 97.06 98.01 99.12 100.06 101.12 102.11 102.94 103.79 104.94 106.09 110.11 112.10 112.93 115.04 118.06 121.02 123.15 126.51 128.14 128.98 130.38 130.91 132.05 133.25 133.81 136.12 138.82 154.73 156.04 164.82 165.97 174.31 179.95".split()] 
      y=[float(s) for s in " 31.06 33.77 38.65 38.14 37.04 35.55 35.06 35.57 37.47 43.28 50.17 55.89 59.44 58.07 60.22 58.68 60.69 59.09 60.52 57.10 59.11 57.24 50.07 64.55 58.02 57.88 58.82 59.09 62.45 62.67 41.12 41.35 64.33 65.78 60.51 61.91 37.54 63.21 64.15 64.63 64.43 66.16 65.21 64.57 63.82 34.75 63.37 61.83 54.62 33.30 33.11 36.62 31.97 30.94 30.41 29.51 29.40 28.89".split()]
    tck = interpolate.splrep(x,y,s=0)
    r =interpolate.splev([theta*180./pi],tck,der=0)
    return r, r*cos(theta), r*sin(theta)

def get_zeros(inp='./'):
    #input=[[float(s) for s in k.split()] for k in 
     # os.popen("plotbands.py pkpt=t %sOUTCAR"%inp+""" | jkawk.py 'if l[0]!="#":pp l[0],l[1]," ".join(["%.8f"%k for k in l[4:]])'""").readlines()]
    os.system("plotbands.py pkpt=t %sOUTCAR"%inp+"""|jkawk.py 'if l[0]!="#":pp l[0],l[1]," ".join(["%.8f"%k for k in l[4:]])' >temp11""")
    input=[[float(s) for s in k.split()] for k in open('temp11')]
    temph=["%.8f %.8f"%(s[0],s[1]) for s in input if len([k for k in s[2:] if k<=0])<=3 ]
    tempe=["%.8f %.8f"%(s[0],s[1]) for s in input if len([k for k in s[2:] if k<=0])>=5 ]
    OUT=open("fermih.out",'w');OUT.write("\n".join(temph)+'\n ');OUT.close()
    OUT=open("fermie.out",'w');OUT.write("\n".join(tempe)+'\n ');OUT.close()
    if temph and tempe: return 1
    else: return 0

def runvasp(inp):
  incar="NSW= 0 \nISMEAR = 0 \nPREC = Normal\nLWAVE=.FALSE.\nNEDOS=2000\nSIGMA=0.20\n"
  incar_relax="NSW= 41 \nISMEAR = 0 \nIBRION=2 \nISIF=2 \nPREC = Normal\nLCHARG=.FALSE.\nNEDOS=2000\nSIGMA=0.20\n"
  if os.uname()[1].split(".")[-1]=="scsystem":incar+="IALGO=48\n";incar_relax+="IALGO=48\n"
  kpoints="coord.13.02\n 0\nGamma point shift\n %s %s 1\n 0 0 0\n"
  kinput="0 0 0\n0.5 0 0\n0 0.5 0\n0 0 0\n-0.5 1 0\n"
  def getpot():os.system("cp /home/soo/vasp/POTCAR_pawlda/C_s/POTCAR .")
  if inp=='relax':
    # relaxation run...
    os.mkdir("relax"); os.chdir("relax")
    pos.pposcar("POSCAR")
    OUT=open("INCAR",'w'); OUT.write(incar_relax);OUT.close()
    OUT=open("KPOINTS",'w'); OUT.write(kpoints%(oo.krelax,oo.krelax));OUT.close(); getpot()
    print "\n\033[46;30m step %s\nperforming relaxation calculation...  \033[0m\n"%nstep
    #os.system("vasp 2>/dev/null | tee vasp.out")
    os.system("%s 2>/dev/null | tee vasp.out"%vaspexe)
    temp=get_zeros()
    os.chdir("..")
    return temp
  elif inp=='mesh':
    os.mkdir("mesh"); os.chdir("mesh")
    os.system("cp ../relax/CONTCAR POSCAR")
    generate_kmesh()
    OUT=open("INCAR",'w'); OUT.write(incar);OUT.close(); getpot()
    print "\n\033[46;30m computing mesh eigenvalues...  \033[0m\n"
    os.system("%s 2>/dev/null | tee vasp.out"%vaspexe)
    temp=get_zeros()
    os.chdir("..")
    return temp
  elif inp=='pathb':
    print "\n\033[46;30m Compute bands along path...  \033[0m\n"
    #temp=sorted(glob("rho*"))[-1]
    os.mkdir("path_bands"); os.chdir("path_bands")
    OUT=open("INCAR",'w'); OUT.write(incar);OUT.close(); getpot()
    #os.system('cp ../%s/static/CONTCAR  POSCAR;cp ../%s/static/CHGCAR .;echo ICHARG=11 >> INCAR'%(temp,temp))
    if oo.mesh:os.system('cp ../mesh/CONTCAR  POSCAR;cp ../mesh/CHGCAR .;echo ICHARG=11 >> INCAR')
    else:os.system('cp ../static/CONTCAR  POSCAR;cp ../static/CHGCAR .;echo ICHARG=11 >> INCAR')
    generate_highk()
    os.system("%s 2>/dev/null | tee vasp.out"%vaspexe)
    temp=get_zeros('../static/')
    os.chdir("..")
    return temp
  elif inp=='static':
    print "\n\033[46;30m static run...  \033[0m\n"
    os.mkdir("static"); os.chdir("static")
    OUT=open("INCAR",'w'); OUT.write(incar);OUT.close()
    OUT=open("KPOINTS",'w'); OUT.write(kpoints%(oo.kstatic,oo.kstatic));OUT.close(); getpot()
    os.system('cp ../relax/CONTCAR  POSCAR;mv ../relax/WAVECAR .')
    os.system("%s 2>/dev/null | tee vasp.out"%vaspexe)
    os.system("rm WAVECAR")
    temp=get_zeros()
    os.chdir("..")
    return temp
  else:print "sorry... vasp run option %s not understood"%inp;sys.exit()

def bounds(rr,npts,diff,lb,ub):
  # compute lower and ubber bounds if possible...
  if npts<oo.crit:lb=rr
  else:ub=rr  #lots of points and points are some distance from one another, back up a bit
  if not ub:rr=lb+oo.sinc
  elif not lb:rr=ub-oo.sinc
  else:rr=(ub+lb)/2.; diff=ub-lb
  print "lower bound = %s upper bound = %s"%(lb,ub)
  return (rr,npts,diff,lb,ub)

def makestruc(xx,yy,rot=''):
  poscar="graphene\n1\n%.8f %.8f 0\n0 %.8f 0\n0 0 27.94437798\n2\ndirect\n0 0 0 C\n0.33333333 0.66666667 0.0 C"
  d=1.41199; a=2*d*cos(30*pi/180)
  a1x,a1y,a2y=(sqrt(3)/2*a,-a/2,a)
  pos=structure(poscar%(a1x,a1y,a2y))
  if rot:pos.rot(rot)
  pos.axial_strain('%s %s 0'%(xx/100.,yy/100.))
  return pos

def pad(inp):
    if inp<10:jj='00'+str(inp)
    elif inp<100:jj='0'+str(inp)
    else:jj=str(inp)
    return jj  

if __name__ == "__main__":
  prr=20
  for theta in [0,10,20,30,40,50,55,60,65,70,80,90]:
      try:os.mkdir("theta_"+pad(theta))
      except: print "directory theta %s already exists... skipping"%theta;continue
      os.chdir("theta_"+pad(theta))
      rr=prr
      xx,yy=(rr*cos(theta*np.pi/180.),rr*sin(theta*np.pi/180.))
      lb=ub='';diff=2*oo.tol; nstep=0
      while diff>oo.tol and nstep<oo.maxiter:
        dn="rho_%s_%.2f"%(pad(nstep),rr)
        try:os.mkdir(dn)
        except: print "directory %s already exists... skipping"%pad(theta);break
        os.chdir(dn);print "RUNNING DIRECTORY %s"%dn
        pos=makestruc(xx,yy,oo.rot) # define the structure...
        npts=runvasp('relax') # relax first...
        print "Electron-Hole puncture from relax: %s"%npts
        if oo.mesh: npts+=runvasp('mesh'); print "Electron-Hole puncture from mesh: %s"%npts
        else:npts+=runvasp('static');print "Electron-Hole puncture from static: %s"%npts
        if dir(oo).count('pathb'):npts+=runvasp('pathb');print "Electron-Hole puncture from path_bands: %s"%npts
        rr,npts,diff,lb,ub=bounds(rr,npts,diff,lb,ub) # check bounds...
        xx=rr*cos(theta*pi/180.); yy=rr*sin(theta*pi/180.)
        nstep+=1; os.chdir("..")
      OUT=open("output",'w'); OUT.write("%s %s   %s %s   %s %s\n"%(theta,rr,rr*cos(theta*pi/180.),rr*sin(theta*pi/180.),lb,ub)); OUT.close()
      os.chdir("..");prr=rr
