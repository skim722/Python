#!/usr/bin/env python
import os,sys,re
from numpy import *
from numpy.linalg import *
from math import *

##########################################################################################
anumb={'V':23,'O':8,'Li':3,'Na':11,'Co':27,'La':57,'Ti':22,'Cd':48,'Au':79,'H':1,'He':2,'Be':4,'B':5,'C':6,'N':7,'F':9,'Ne':10,'Mg':12,'Al':13,'Si':14,'P':15,'S':16,'Cl':17,'Ar':18,'K':19,'Ca':20,'Sc':21,'Cr':24,'Mn':25,'Fe':26,'Co':27,'Ni':28,'Cu':29,'Zn':30,'Ga':31,'Ge':32,'As':33,'Se':34,'Br':35,'Kr':36,'In':49,'NULL':0}
class structure:
 """A general class to store, manipulate, and print a lattice."""
 def __init__(self,poscar): 
  """takes a poscar file as a list (ie. split by line), or as a string seperated by line breaks. 
  then the following info is stored:\
  \nvec=lattice vectors\natd=direct atomic coordinates\natc=cartesian atomic coordinates\
  \ntitle\nscale\nnatoms\nspecpoint=list with the # of atoms of each species\
  \natpoint=list with species number for each atom \nnames=name of each atom\nspecname=name of each species"""
  if type(poscar)==type(''):
    if os.path.isfile(poscar):poscar=open(poscar).readlines()
    else:poscar=poscar.split("\n")
  elif type(poscar)!=type([]): print "structure constructor takes string or list...";sys.exit()
  self.poscar=poscar
  self.title=poscar[0].strip()
  self.scale=float(poscar[1])
  # get lattice vectors
  self.vec=zeros((3,3),float)
  for i in range(3):self.vec[i]=array([float(s) for s in poscar[2+i].split()])
  # rescale the lattice vectors...
  self.vec=self.vec*self.scale
  # get the species-number pointer
  self.nspecies=len(poscar[5].split())
  self.specpoint=array([int(s) for s in poscar[5].split()])
  self.atpoint=[s+1 for s in range(self.nspecies) for j in range(self.specpoint[s])]
  # get total number of atoms
  self.natoms=sum(self.specpoint)
  # direct or cartesian?
  self.coord=poscar[6].strip().lower()[0]
  # check to make sure the poscar has proper number of lines...
  pl=[s for s in range(len(poscar)) if not poscar[s].split()]
  if not pl:pl=len(poscar)
  else: pl=pl[0]
  if pl!=self.natoms+7: print "POSCAR has wrong number of lines";sys.exit()
  # get basis atom coordinates....
  self.at=zeros((self.natoms,3),float)
  for i in range(self.natoms): self.at[i]=array([float(s) for s in poscar[i+7].split()[:3]])
  if self.coord=='c':self.atc=self.at*self.scale; self.atd=dot(self.atc,inv(self.vec)) ; self.at=self.atc
  elif self.coord=='d':self.atd=self.at; self.atc=dot(self.atd,self.vec)
  else: print "ERROR: coordinates must be direct or cartesian";sys.exit()
  self.scale=1.0 # the vec and atc were rescaled above... set scale to 1...
  self.cellinfo()
  # get the names of atoms if given... otherwise give defaults...
  self.names=[]; self.specname=[]
  for i in range(self.natoms): 
    if len(poscar[7+i].split())>3:self.names.append(poscar[7+i].split()[3]) 
    else: self.names.append("sp%s"%self.atpoint[i])
  for i in self.names:
   if not self.specname.count(i):self.specname.append(i)

 def copy(self): return structure(self.pposcar(True))
 
 def cellinfo(self):
  """Caclulates the cell volume and the vector lengths and anglges. """
  self.vol=dot(self.vec[2],cross(self.vec[0],self.vec[1]))
  # get lengths and angles of vectors...
  self.abc=[sqrt(dot(s,s)) for s in self.vec]; self.abg=[0]*3
  self.abg[0]=acos(dot(self.vec[1],self.vec[2])/(self.abc[1]*self.abc[2]))*180/pi
  self.abg[1]=acos(dot(self.vec[0],self.vec[2])/(self.abc[0]*self.abc[2]))*180/pi
  self.abg[2]=acos(dot(self.vec[0],self.vec[1])/(self.abc[0]*self.abc[1]))*180/pi
  # define reciprocal lattice vectors...
  self.rvec=zeros((3,3),float)
  self.rvec[0]=cross(self.vec[1],self.vec[2])/self.vol
  self.rvec[1]=cross(self.vec[2],self.vec[0])/self.vol
  self.rvec[2]=cross(self.vec[0],self.vec[1])/self.vol
  self.rvec=2*pi*self.rvec        # 2pi convention
  self.rvol=dot(self.rvec[2],cross(self.rvec[0],self.rvec[1]))
  # get lengths and angles of reciprocal vectors...
  self.rabc=[sqrt(dot(s,s)) for s in self.rvec]; self.rabg=[0]*3
  self.rabg[0]=acos(dot(self.rvec[1],self.rvec[2])/(self.rabc[1]*self.rabc[2]))*180/pi
  self.rabg[1]=acos(dot(self.rvec[0],self.rvec[2])/(self.rabc[0]*self.rabc[2]))*180/pi
  self.rabg[2]=acos(dot(self.rvec[0],self.rvec[1])/(self.rabc[0]*self.rabc[1]))*180/pi

  
 def __add__(self,other):
   if type(self)!=type(other): print "cannot add different data types"; sys.exit()
   if self.natoms!=other.natoms and self.nspecies != other.nspecies and self.coord!=other.coord: 
     print "different number of atoms or species or coordinate system";sys.exit()
   new=structure(self.poscar) 
   new.vec=self.vec+other.vec
   new.atc=self.atc+other.atc   
   new.atd=self.atd+other.atd   
   new.at=self.at+other.at   
   new.cellinfo()
   return new
    
 def __sub__(self,other):
   if type(self)!=type(other): print "cannot add different data types"; sys.exit()
   if self.natoms!=other.natoms and self.nspecies != other.nspecies and self.coord!=other.coord: 
     print "different number of atoms or species or coordinate system";sys.exit()
   new=structure(self.poscar) 
   new.vec=self.vec-other.vec
   new.atc=self.atc-other.atc   
   new.atd=self.atd-other.atd   
   new.at=self.at-other.at   
   #new.cellinfo()
   return new
    
 def __mul__(self,other):
   if type(other)!=type(1) and type(other)!=type(1.0): print "must give integer or float"; sys.exit()
   new=structure(self.poscar) 
   new.vec=self.vec*other
   new.atc=self.atc*other
   new.atd=self.atd*other
   new.at=self.at*other
   #new.cellinfo()
   return new
    
 def __div__(self,other):
   if type(other)!=type(1) and type(other)!=type(1.0): print "must give integer or float"; sys.exit()
   new=structure(self.poscar) 
   new.vec=self.vec/other
   new.atc=self.atc/other
   new.atd=self.atd/other
   new.at=self.at/other
   #new.cellinfo()
   return new

 def combine(self,other):
   """Combine two structures. """
   self.specpoint=array(self.specpoint.tolist()+other.specpoint.tolist())
   self.names=self.names+other.names
   self.nspecies+=len(other.specpoint)
   self.atpoint=[s+1 for s in range(self.nspecies) for j in range(self.specpoint[s])]
   self.natoms=sum(self.specpoint)
   self.specname=self.specname+other.specname
   #print self.specname
   #for i in self.names:
   #   if not self.specname.count(i):self.specname.append(i)
   self.atc=array(self.atc.tolist()+other.atc.tolist())
   self.atd=dot(self.atc,inv(self.vec))
   if self.coord=='d':self.at=self.atd.copy()
   elif self.coord=='c':self.at=self.atc.copy()

    
 def shift_atom(self,at,inp,coord='d'):
   """Shift an individual atom by a b c where a,b,c are given in either direct or cartesian coordinates. """
   shift=array([float(eval(s)) for s in inp.split() if s ])
   if len(shift)!=3 : print "origin shift must be three coordinates"; sys.exit()
   if coord=='d':
     self.atd[at]+=shift
     self.atc=dot(self.atd,self.vec)
   if coord=='c':
     self.atc[at]+=shift
     self.atd=dot(self.atc,inv(self.vec)) 
   if self.coord=='d':self.at=self.atd.copy()
   elif self.coord=='c':self.at=self.atc.copy()
    
 def kill_atom(self,at):
   """Remove atoms from cell. This is not properly generalized yet, but it will get the job done for rasmol/povray. """
   temp=[];name=[]
   for i in range(self.natoms):
     if not [int(s) for s in at.split()].count(i):temp.append(self.atd[i]);name.append(self.names[i])
   self.natoms=len(temp)
   self.atd=array(temp)
   self.names=[s for s in name];del name
   self.atc=dot(self.atd,self.vec)
   if self.coord=='d':self.at=self.atd.copy()
   elif self.coord=='c':self.at=self.atc.copy()

 def rescale_atoms(self,scale):
   """Rescales the atoms only. Can be useful in obscure situations."""
   self.atc=self.atc*scale; self.atd=dot(self.atc,inv(self.vec))
   if self.coord=='d':self.at=self.atd.copy()
   elif self.coord=='c':self.at=self.atc.copy()
 
 def rescale(self,scale):
   """Rescales the unit cell. """
   self.vec=self.vec*scale; self.atc=self.atc*scale; self.abc=[scale*s for s in self.abc]; self.cellinfo()
 
 def resetscale(self,scale):
   """Resets scale in structure constant (ie. structure NOT changed).  """
   self.scale=scale; self.vec=self.vec/scale; self.atc=self.atc/scale; self.abc=[s/scale for s in self.abc]; self.cellinfo()

 def shift_origin(self,inp):
   """Shifts the origin of the unit cell by a b c where a,b,c are given in direct coordinates. """
   shift=array([float(eval(s)) for s in inp.split() if s ])
   if len(shift)!=3 : print "origin shift must be three coordinates"; sys.exit()
   self.atd=array([s+shift for s in self.atd])
   self.atc=dot(self.atd,self.vec)
   if self.coord=='d':self.at=self.atd.copy()
   elif self.coord=='c':self.at=self.atc.copy()
 
 def pposcar(self,ret=''):
   """Prints the poscar."""
   poscar=[0]*(self.natoms+7)
   poscar[0]=self.title+'\n'
   poscar[1]="%s\n"%self.scale
   for i in range(3): poscar[2+i]=("%.8f   %.8f   %.8f\n"%tuple(self.vec[i])).replace(" -","-")
   poscar[5]=" ".join([str(s) for s in self.specpoint])+'\n'
   poscar[6]="%s\n"%self.coord
   for i in range(self.natoms): 
     poscar[7+i]=("%.8f   %.8f   %.8f "%tuple(self.at[i])).replace(" -","-")+" %s\n"%self.names[i]
   if ret==1 or ret==True or ret=='t': return "".join(poscar).rstrip()
   elif ret!='' and type(ret)==type(''):OUT=open(ret,'w');OUT.write("".join(poscar).rstrip());OUT.close()
   else: print "".join(poscar).rstrip()
 
 def dirtocar(self): 
   """Changes the print format from direct to cartesian or vice versa."""
   if self.coord=='d':self.at=self.atc; self.coord='c'
   elif self.coord=='c':self.at=self.atd; self.coord='d'

 def pinfo(self): 
   """Prints the volume and vector lengths and angles."""
   print "=========== Lattice Information ========================"
   print "vol=%.5f\na=%.5f b=%.5f c=%.5f\nalpha=%.5f beta=%.5f gamma=%.5f\n"%(tuple([self.vol])+tuple(self.abc+self.abg))
   print "========================================================"

 def prinfo(self): 
   """Prints the reciprocal lattice information."""
   print "=========== Reciprocal Lattice Information ============="
   print (("%.8f %.8f %.8f \n"*3)%tuple(list(self.rvec[0])+list(self.rvec[1])+list(self.rvec[2]))).rstrip()
   print "vol=%.5f\na=%.5f b=%.5f c=%.5f\nalpha=%.5f beta=%.5f gamma=%.5f"%(tuple([self.rvol])+tuple(self.rabc+self.rabg))
   print "========================================================"

 def dist(self,verb=0,cut=2.5,transx=2,transy=2,transz=2,output='print'):
   """Prints out the distances between atoms. There is four levels of verbosity (ie. verb=0-3).\
      \nCut specifies the maximum bond length displayed. The trans paramters should not need to be changed.\
      \nThese set how many translations are taken to look for neighboring atoms."""
   comp=[]; natoms=self.natoms; at=self.atc; v=self.vec; name=self.names; self.distance=[]
   # loop over atoms in the unit cell...
   for at1 in range(natoms):
     for tx in range(-transx,transx+1):
       for ty in range(-transy,transy+1):
         for tz in range(-transz,transz+1):
           for at2 in range(natoms):
             dv=at[at1]-(at[at2]+tx*v[0]+ty*v[1]+tz*v[2])
             dist=sqrt(dot(dv,dv)) ; dv=list(dv)
             # need to format dv for printing
             for i in range(3):
               if dv[i]<0:dv[i]="%.5f"%dv[i]
               else:dv[i]=" %.5f"%dv[i]
             if verb==0: comp=comp+[[dist,"%2s%s---%2s%s  | dist=%.5f"%(name[at1],at1+1,name[at2],at2+1,dist)]]
             elif verb==1: comp=comp+[[dist,"%2s%s---%2s%s  | t1=%2s t2=%2s t3=%2s |  dist=%.5f"%(name[at1],at1+1,name[at2],at2+1,tx,ty,tz,dist)]]
             elif verb==2: comp=comp+[[dist,"%2s%s---%2s%s  | d1=%s d2=%s d3=%s |  dist=%.5f"%(name[at1],at1+1,name[at2],at2+1,dv[0],dv[1],dv[2],dist)]]
             elif verb==3: comp=comp+[[dist,"%2s%s---%2s%s  | t1=%2s t2=%2s t3=%2s | d1=%2s d2=%2s d3=%2s |  dist=%.5f"\
                                %(name[at1],at1+1,name[at2],at2+1,tx,ty,tz,dv[0],dv[1],dv[2],dist)]]
     comp.sort()
     # store the distances for future use...
     self.distance+=comp
     # print if requested...
     if output=='print':
       #for i in comp:print i[1]
       for i in comp:
         if 0<i[0]<=cut:print i[1]
       print ""
     comp=[]
   # final storage of distance information...
   self.distance=[s for s in sorted(self.distance) if 0<s[0]<=cut]  
   
 
 def xyz(self,call=0,cols="blue green red yellow cyan orange purple",scr=0):
   """Print out the xyz file for viewing in rasmol."""
   if call=='rl':call='r';scr=1
   if call=='jl':call='j';scr=1
   home=os.getenv('HOME')
   cols=cols.split()
   comp='%s\n'%self.natoms
   comp+='%s\n'%self.title
   for i in range(self.natoms): comp+='%s %s\n'%(self.names[i],(" %.8f   %.8f   %.8f"%tuple(self.atc[i])).replace(" -","-"))
   OUT=open("%s/.ras.xyz"%home,"w") 
   print >>OUT, comp; OUT.close()
   # make custom scriptfile to set bond size/color and atom colors...
   OUT=open("%s/.rscript"%home,'w') 
   if scr==0:OUT.write("wireframe 30\nspacefill 110\ncolor bonds white\n#set axes true\n")
   elif scr==1:OUT.write("wireframe 2\nspacefill 10\ncolor bonds white\n")
   for i,ii in enumerate(self.specname):
     if not anumb.keys().count(ii):print "ooppss. please add %s to the anumd dictionary of periodic table."%ii;sys.exit()
     OUT.write("select elemno=%s\ncolor %s\n"%(anumb[ii],cols[i] if i<len(cols) else 'purple'))
   OUT.close()
   # prepare for rasmol...
   if call=='r': 
     if not os.path.exists(home+"/bin/rasmol"):print "Could not find rasmol in your bin...";sys.exit()
     print "rasmol -script %s/.rscript -xyz %s/.ras.xyz"%(home,home)
     os.system("rasmol -script %s/.rscript -xyz %s/.ras.xyz"%(home,home))
   # prepare for jmol...
   elif call=='j': 
     # option to make bonds and atoms really small...
     if not os.path.exists(home+"/bin/Jmol.jar"):print "Could not find Jmol.jar in your bin...";sys.exit()
     print "java -jar %s/bin/Jmol.jar %s/.ras.xyz -s %s/.rscript"%(home,home,home)
     os.system("java -jar %s/bin/Jmol.jar %s/.ras.xyz -s %s/.rscript"%(home,home,home))
   elif call=='0' or call==0:print comp.rstrip()
   else:print "call must be 'r', 'j', 'jl', or 0. stopping.";sys.exit()

 def allincell(self,tol=10):
   """Translate basis atoms such that they all exist within unit cell."""
   atnew=self.atd.copy()
   for i in range(len(atnew)):
     for j in range(len(atnew[0])):
       atnew[i,j]=sqrt((round(atnew[i,j],tol)%1)**2) # sq-sqrt is to get rid of -0... just annoyance.
       #if round(atnew[i,j],tol)<0:     atnew[i,j]+=1
       #elif round(atnew[i,j],tol)>=1:  atnew[i,j]-=1
   self.atd=atnew.copy()
   self.atc=dot(self.atd,self.vec)
   if self.coord=='d':self.at=self.atd.copy()
   else:self.at=self.atc.copy()

 def outofcell(self,tol=0.98):
   """Translate basis atoms with direct coordinates above tol outside the unit cell."""
   atnew=self.atd.copy();tol=float(tol)
   for i in range(len(atnew)):
     for j in range(len(atnew[0])):
       if atnew[i,j]>tol:     atnew[i,j]-=1
   self.atd=atnew.copy()
   self.atc=dot(self.atd,self.vec)
   if self.coord=='d':self.at=self.atd.copy()

 def completecellxyz(self,inp=''):
   """Complete the unit cell for viewing 
      \n(ie. add a translated atom for each component of any atom with a zero direct coordinate).
      \nThe input determines which directions one completes in. It should be understood that one is
      \nleft with a cluster and NOT a unit cell."""
   #self.allincell() # start with all atoms in unit cell
   if not inp:inp=[1,1,1]
   elif len(inp)==3 and [s for s in inp if s=='0' or s=='1']:inp=[int(s) for s in inp]
   else:print "improper input for completecellxyz... give 100,110,111 or permutations thereof.";sys.exit()
   tol=10**-6; out=[]; np=[]; name=self.names
   for i,at in enumerate(self.atd):
     #np.append(name[i]); out.append(at.tolist())
     for vv in [array(s) for s in [1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,1],[1,1,1]]: 
       for jj in [0,1]: # check for atoms near 0 and 1...
         if len([ss for s,ss in enumerate(at) if vv[s] and abs(ss-1*jj)<tol])==sum(vv):
           if not [s for s,ss in enumerate(vv) if ss>0 and inp[s]==0]:
             out.append((at+vv).tolist())
             np.append(name[i])

   np=[s for s in name+np]; out=[s for s in self.atd.tolist()+out]
   self.title+=' NOT A UNIT CELL: ONLY A CLUSTER'
   self.names=np
   self.atpoint=[self.specname.index(s)+1 for s in self.names]
   self.specpoint=[self.atpoint.count(s) for s in range(1,self.nspecies+1)]
   self.natoms=len(out)
   self.atd=array(out)
   self.atc=dot(array(out),self.vec)
   if self.coord=='c':self.at=self.atc
   elif self.coord=='d':self.at=self.atd


 def label(self,inp):
   "Labels the atoms in the unit cell. must give a name for each species as the string 'spec1 spec2...'."
   specname=inp.replace("_"," ").split()
   if len(specname)!=self.nspecies: print "improper number of names give"; sys.exit()
   self.names=[specname[s-1] for s in self.atpoint]; self.specname=[]
   self.specname=[s for s in specname]
   #for i in self.names:
   #  if not self.specname.count(i):self.specname.append(i)

 def rot(self,inp,trans=''):
   "Rotates the unit cell. user must supply 3x3 rotation matrix as the string 'r11 r12 r13 r21 r22 r23...'"
   inp=inp.replace("_"," ")
   if len(inp.split())!=9:print "must give a proper rotation matrix";sys.exit()
   rot=reshape(array([eval(s) for s in inp.split()]),(3,3))
   if trans:rot=transpose(rot)
   # make sure matrix has zero determinant...
   if abs(det(rot)-1)>=10**-6:
     sys.stderr.write("#"*20+"\nrotation matrix must have determinant of 1... normalizing rows.\n"+"#"*20+"\n")
     rot=array([s/sqrt(dot(s,s)) for s in rot])
   self.vec=transpose(dot(rot,transpose(self.vec)))
   self.atc=dot(self.atd,self.vec) # be sure to update the cartesian coordinates...
   if self.coord=='c':self.at=self.atc 
   self.cellinfo()

 def axial_strain(self,inp):
   "Strains the unit cell. user must supply x,y,z strains as the string 's1 s2 s3'"
   inp=inp.replace("_"," ").strip().split()
   if len(inp)!=3:print "must give a proper axial strain";sys.exit()
   strain=zeros((3,3),float)
   for i in range(3):strain[i,i]=float(inp[i])+1
   self.vec=dot(self.vec,strain)
   self.atc=dot(self.atd,self.vec) # be sure to update the cartesian coordinates...
   if self.coord=='c':self.at=self.atc 
   self.cellinfo()

 def old_supa(self,inp):
   "Creates a supercell. user must supply 3x3 integer supercell matrix as the string 'r11 r12 r13 r21 r22 r23...'"
   # this was written years ago and is very slow... replaced below...
   if len(inp.split())!=9:print "must give a proper supercell matrix";sys.exit()
   supa=reshape(array([float(s) for s in inp.split()]),(3,3))
   v=self.vec; natoms=self.natoms; nspecies=self.nspecies; at=self.atc
   # construct superlattice
   vsup=dot(supa,v)

   # find volume of superlattice 
   volsup=dot(vsup[2],cross(vsup[0],vsup[1]))
   # this is a total numpy/python bug! 
   # without the round you can get wrong answer...
   mult=abs(int(round(volsup/self.vol)))

   # find maximum supa entry
   max=0
   for i in range(3):
     for j in range(3):
       if supa[i,j] > max:
          max=supa[i,j]
   # i was missing atoms with the above algorithm... so mult by 2!
   max=int(max*2)

   # translate all basis atoms by multiples of the translations less than cuttoff.
   count=zeros(natoms)
   newat=[]
   for i in range(natoms):
       for k in range(-max,max+1):
         for l in range(-max,max+1):
           for m in range(-max,max+1):
             tr_at_car=at[i,:]+v[0,:]*k+v[1,:]*l+v[2,:]*m
             tr_at_dir=dot(tr_at_car,inv(vsup))
             # this next line is critical for avoiding errors...
             for pp in range(3): tr_at_dir[pp]=float("%12.10f"%tr_at_dir[pp])
             if 0<=tr_at_dir[0]<1 and 0<=tr_at_dir[1]<1 and 0<=tr_at_dir[2]<1:
                count[i]+=1
                newat.append(tr_at_dir)
   for i in range(natoms):
     if count[i] != mult:
        print "WARNING: looks like algorithm did not find proper number of atom %s"%i
   # update the attributes...
   self.atd=array(newat); self.vec=vsup; self.specpoint*=mult
   self.atc=dot(self.atd,self.vec); self.natoms*=mult; self.vol*=mult
   if self.coord=='c':self.at=self.atc
   elif self.coord=='d':self.at=self.atd
   self.atpoint=[s+1 for s in range(self.nspecies) for j in range(self.specpoint[s])]
   self.specname=[]
   for i in self.names:
     if not self.specname.count(i):self.specname.append(i)
   self.label(" ".join(self.specname))
   self.cellinfo()

 def supa(self,inp,allincell=False):
   "Creates a supercell. user must supply 3x3 integer supercell matrix as the string 'r11 r12 r13 r21 r22 r23...'"
   if len(inp.split())!=9:print "must give a proper supercell matrix";sys.exit()
   supa=reshape(array([int(s) for s in inp.split()]),(3,3))
   v=self.vec; natoms=self.natoms; nspecies=self.nspecies; at=self.atc
   # construct superlattice
   vsup=dot(supa,v)

   # find volume of superlattice 
   volsup=dot(vsup[2],cross(vsup[0],vsup[1]))
   # this is a total numpy/python bug! 
   # without the round you can get wrong answer...
   mult=abs(int(round(volsup/self.vol)))

   # find all prim translations within the supercell...
   trlist=[]
   for i in range(sort(supa[:,0])[0],sort(supa[:,0])[-1]+1):#+1 if sort(supa[:,0])[-1]<=0 else 0):
     for j in range(sort(supa[:,1])[0],sort(supa[:,1])[-1]+1):#+1 if sort(supa[:,0])[-1]<=0 else 0):
       for k in range(sort(supa[:,2])[0],sort(supa[:,2])[-1]+1):#+1 if sort(supa[:,0])[-1]<=0 else 0):
         temp=dot(dot((i,j,k),v),inv(vsup))
         #if len([ s for s in temp if 0.<=round(s,6)<1.])==3:trlist.append((sqrt(dot(dot(temp,vsup),dot(temp,vsup))),tuple(temp)))
         if len([ s for s in temp if 0.<=round(s,6)<1.])==3:trlist.append(tuple(temp))
   if len(trlist)!=mult:print "did not find all translations withing supercell... algorithm problem.";sys.exit()
   trlist.sort()
   # get all new atoms...
   newat=[]
   for a in at:
     for tt in trlist: newat.append(dot(a,inv(vsup))+tt) 
   # update the attributes...
   self.atd=array(newat); self.vec=vsup; self.specpoint=array([s*mult for s in self.specpoint])
   self.atc=dot(self.atd,self.vec); self.natoms*=mult; self.vol*=mult
   if self.coord=='c':self.at=self.atc
   elif self.coord=='d':self.at=self.atd
   self.atpoint=[s+1 for s in range(self.nspecies) for j in range(self.specpoint[s])]
   # it was completely unecessary to find specname again... cant change. 
   # leave commented for a while to be sure... delete eventually.
   #self.specname=[]
   #for i in self.names:
   #  if not self.specname.count(i):self.specname.append(i)
   self.label(" ".join(self.specname))
   self.cellinfo()
   if allincell==True:self.allincell()
   # for backwards compatibility with old slow supa routine... can delete eventually...
   for i in range(natoms):self.atd[mult*i:mult*i+mult,:]=array(sorted(self.atd[mult*i:mult*i+mult,:].tolist()))
   self.atc=dot(self.atd,vsup)
   if self.coord=='c':self.at=self.atc
   elif self.coord=='d':self.at=self.atd



 def ptart(self):
   "print out the structure in lmtart format"
   at=self.atd; v=self.vec; natoms=self.natoms

   # format the vectors and basis for output...
   vec=["%12.8f,%12.8f,%12.8f"%tuple(s) for s in v]
   vec=[s.lstrip() for s in vec]
   vec="\n".join(vec)
   bas=["%12.8f,%12.8f,%12.8f"%tuple(s) for s in at]
   bas=[s.lstrip() for s in bas]
   bas="\n".join(bas)

   print """<FILE=STRFILE,INPUT=MODERN>
***********************************************************
<SECTION=CTRS>          ! CONTROL STRUCTURE:
 nAtom =%s               ! # of atoms
 BtoA = 1.00000         ! b over a ratio
 CtoA = 1.00000         ! c over a ratio
 iCrd =2                ! controls coordinate system
<SECTION=TRAN>          ! PRIMITIVE TRANSLATIONS:
 %s
<SECTION=BASS>          ! BASIS ATOMS :
%s"""%(natoms,vec,bas)

 def pstut(self):
   """Print out the structure in Stuttgart format."""
   v=self.vec; natoms=self.natoms; nspecies=self.nspecies; at=self.atc
   ctrl="""HEADER    %s\
   \nVERS      LMASA-47\
   \nSTRUC     ALAT=1.89\
   \n          PLAT=%.8f %.8f %.8f\
   \n               %.8f %.8f %.8f\
   \n               %.8f %.8f %.8f"""%(tuple([self.title])+tuple(v[0])+tuple(v[1])+tuple(v[2]))

   ctrl+="\nCLASS     ATOM=%s   Z=%s   "%(self.specname[0],anumb[self.specname[0]])
   for i in range(1,nspecies):ctrl+="\n          ATOM=%s   Z=%s   "%(self.specname[i],anumb[self.specname[i]])
   ctrl+="\nSITE      ATOM=%s   POS= %.8f %.8f %.8f"%(tuple([self.names[0]])+tuple(at[0]))
   for i in range(1,natoms):ctrl+="\n          ATOM=%s   POS= %.8f %.8f %.8f"%(tuple([self.names[i]])+tuple(at[i]))

   print ctrl

 def pdmftlab(self):
   "print out the structure in DMFTlab format"
   out="""\
lattice_unit=[
(%.8f,%.8f,%.8f),   
(%.8f,%.8f,%.8f),   
(%.8f,%.8f,%.8f)]
recip_lat=[
(%.8f,%.8f,%.8f),   
(%.8f,%.8f,%.8f),   
(%.8f,%.8f,%.8f)]
orbital_postitions=[
"""%tuple(list(ravel(self.vec))+list(ravel(self.rvec)))
   for i in self.at:
     out+="(%.8f,%.8f,%.8f),\n"%tuple(i) 
   out=out[:-2]+"]\n"  
   print out



 def povray(self,inp='',call=''):
   "create a povray input file from the current structure. "
   from mysub import update_dict
   dict={"cp":'0,0,-50',"cps":"","cd":'0,0,0',"lp":'',"ss":0.3,"br":"0,0,0","amb":0.5,"ph":0.8,'cz':0.0,
         'pl':'','drawcell':'','csize':'0.06','ccol':'0.3,0.3,0.3','camb':0.5,'bsize':'0.06','bcol':'0.9,0.9,0.9','bamb':0.5,
         'centercell':'t',"sky":'0,0,1','lt':'','ct':'orthographic angle','skysph':'','bonds':'','maxb':'2.5',"centercluster":'',
         'w':800,'h':600,'right':'1,0,0','rightsc':'4/3','up':'0,0,1','cmt':'bozo','turb':0.0,'addn_trmit':0.7,'addn_color':'',
         'screenx':'','out':'temp','strans':'0.0','btrans':'0','ctrans':0,'shift':'0,0,0 ',
         # now for the atom colors...
         "Na":"1,0,0","Na1":"1,0,0","Na2":"1,1,0","Co":"0,0,1","O":"0,1,0","Co1":"0,0,1","C":"0,0,1","La":"1,0,0.2",
         "H":"1,1,1","Pu":"0,0,1","a":"1,0,0","b":'0,0,1','c':'glass','j':'Jade',"Ni":"1,1,0","Al":"0,1,1","NULL":"0,0,0",
         "NULL1":"0.5,0.5,0.5",
         # now for orbital colors...
         'a1g':"0,0.8,0",'egp1':'0.8,0,0.1','egp2':'0.2,1,1.0',
         "eg1":"0,1,1",'eg2':"0,0.8,0",'a1g_cubic':"0,0.8,0",'egp1_cubic':'0.8,0,0.1','egp2_cubic':'0.2,1,1.0',
         'xy':'1,0.1,0.1','xz':'0.1,1,0','yz':'0.2,1.0,1','x2y2':'0.4,1,0.6','z2':"0,0.8,0",'test':'0.1,0,1.0',
         'xyz':'1,0.1,0.1','in_xyz':'0,1.0,0.1',
         }
   if os.path.isfile(inp):inp=open(inp).read() # read in commands from inp if file...
   dict=update_dict(dict,[s for s in inp.split() if s.lstrip()[0]!='#']) # update the dictionary
   for i in dict.keys():  # execute the dictionary
     if type(dict[i])!=type(''): exec "%s=%s"%(i,dict[i])
     else: exec "%s='%s'"%(i,dict[i])
   # be sure element has a colormap...
   for i in self.specname:
     if not dict.keys().count(i):print "color not given for element %s"%i; sys.exit()
   # a list of orbitals...
   orbs={'a1g':0.006,'egp1':0.004,'egp2':0.004,'eg1':0.004,'eg2':0.004,
         'a1g_cubic':0.004,'egp1_cubic':0.004,'egp2_cubic':0.004,'test':0.004,
         'px':'','py':'','pz':'','xy':0.004,'xz':0.004,'yz':0.004,'x2y2':0.004,'z2':0.006,'xyz':0.0025}

   atc=self.atc
   vec=self.vec
   mc=array([0.0,0.0,0.0])
   mc=array([float(s) for s in dict['shift'].split(',')])
   # out the center of mass of unit cell to zero if desired...
   if centercluster:self.clustercenter() 
   elif centercell: mc+=self.cellcenter(); atc=[s-mc for s in self.atc] 
   
   preamble="""\
   #include "colors.inc"
   #include "shapes.inc"
   #include "textures.inc"

   #declare map1 =    color_map {
          [0.00 color rgb <1, 0.25, 0.25>]
          [0.33 color rgb <1, 1, 0.25>]
          [0.50 color rgb <0.25, 1, 0.25>]
          [0.66 color rgb <0.25, 0.25, 1>]
          [1.00 color rgb <1, 0.25, 0.25>]
       }
 
    #declare pig1 = pigment { agate    color_map {map1} turbulence %(turb)s }
    #declare pig2 = pigment { bozo     color_map {map1} turbulence %(turb)s }
    #declare pig3 = pigment { granite  color_map {map1} turbulence %(turb)s }
    #declare pig4 = pigment { leopard  color_map {map1} turbulence %(turb)s }
    #declare pig5 = pigment { marble   color_map {map1} turbulence %(turb)s }
    #declare pig6 = pigment { onion    color_map {map1} turbulence %(turb)s }
    #declare pig7 = pigment { wood     color_map {map1} turbulence %(turb)s }
    #declare glass = pigment { color rgb <1, 0.25, 0.25> transmit 1.0 }

   """%dict


   background="""
   background {
     color rgb <%s>
   }
   """%br
   if cps:
     tt=[float(s) for s in cps.split(",")];tt[1]*=pi/180.;tt[2]*=pi/180.
     dict['cp']=cp="%s,%s,%s"%(tt[0]*cos(tt[1])*sin(tt[2]),tt[0]*sin(tt[1])*sin(tt[2]),tt[0]*cos(tt[2]))
   else:
     x,y,z=[float(s) for s in cp.split(",")]
     r=sqrt(x**2+y**2+z**2)
     if x==0 and y>0:theta=90
     elif x==0 and y<0:theta=-90
     elif x==0 and y==0: theta=0
     else:theta=atan(y/x)*180/pi
     if r==0:phi=0
     else:phi=acos(z/r)*180/pi
     print "(r,theta,phi)  cps=%s,%s,%s "%(r,theta,phi)
   if cz:cp="%s,%s,%s"%(tuple([float(s)*cz for s in cp.split(',')]))
   if ct=='p':dict['ct']=ct='perspective'
   camera="""\
   camera {
     %s
     location <%s>
     sky <%s>
     look_at <%s>
     right <%s>*%s
     up <%s>
   }
   """%(ct,cp,sky,cd,right,rightsc,up)

   if not lp:lp=cp
   if ct.count('orthographic') and not lt:lt='shadowless'
   light="""\
   light_source {
     <%s>
     color White
     %s
   }
   """%(lp,lt)

   sphere="""\
   sphere {
     <%(sph_pos)s>, %(ss)s
     pigment { %(sph_col)s %(sph_trmit)s transmit %(strans)s}
     finish { phong %(ph)s ambient %(amb)s   }
     //finish { ambient 0.5 diffuse 0.2 phong 0.5 reflection  {0.15} metallic 0.5 }
   }
   """

   cylinder="""\
   cylinder {
        <%s,%s,%s>, <%s,%s,%s>, %s
        pigment { color rgb <%s> transmit %s }
        finish  { ambient %s }
   }
   """
 
   plane="""\
   plane { // the floor
     z, %s  // along the x-z plane (y is the normal vector)
     pigment { checker color rgb <0.2,0.2,0.2> color Black } // checkered pattern
     finish { diffuse 0.25  ambient 0.6 }
   }
   """

   # this is the old primitive version which had to reload every time orbital was plotted
   # now we use the better version below...
   orbital="""\
   #declare TRANSLATE = transform {translate <%s,%s,%s> scale %s} ;
   #declare FRONT_TEXTURE = texture {
    pigment { color rgb <%s> }
    //normal { bumps 2.4 scale 0.5 }
    finish { ambient 0.5 diffuse 0.4 phong 1 reflection  {0.1} metallic 1 }
   }
   #declare BACK_TEXTURE = texture {
    pigment { color rgb <%s> }
    finish { ambient 0.5 diffuse 0.4 phong 1 reflection  {0.1} metallic 1 }
   }
   #include "orb_fine_%s.inc"
   """

   orbitalnew="""
   object{ %s 
   translate <%s,%s,%s> scale %s
   texture { 
       pigment { color rgb <%s> }
       finish { ambient 0.5 diffuse 0.4 phong 1 reflection  {0.1} metallic 1 }
       }
   interior_texture { 
       pigment { color rgb <%s> }
       finish { ambient 0.5 diffuse 0.4 phong 1 reflection  {0.1} metallic 1 }
       }

   }
   """

   skysphere="""\
  sky_sphere {
    pigment {
      gradient y
      color_map {
        [0.000 0.002 color rgb <1.0, 0.2, 0.0>
                     color rgb <1.0, 0.2, 0.0>]
        [0.002 0.200 color rgb <0.8, 0.1, 0.0>
                     color rgb <0.2, 0.2, 0.3>]
      }
      scale 2
      translate -1
    }
    rotate -135*x
  }
   """

   # concatenate the various commands...
   comp=preamble+background+camera+light
   if skysph:comp+=skysphere
   if pl:comp+=plane%pl
   # if there are any orbitals present, put in there include statement
   for i in self.specname:
     if orbs.keys().count(i): comp+='\n   #include "orb_fine_%s.inc"'%i
   # add in the orbitals and spheres...
   for i,jj in enumerate(atc): 
     if orbs.keys().count(self.names[i]): # this is an orbital... 
       ii=jj/orbs[self.names[i]] # rescale the translations because we will rescale the orbitals...
       scale=orbs[self.names[i]]; incol=dict[self.names[i]] ; excol=dict[self.names[i]]; name=self.names[i]
       if dict.has_key("in_"+self.names[i]):incol=dict["in_"+self.names[i]] # allow for different interior color...
       comp+=orbitalnew%(name,ii[0],ii[1],ii[2],scale,excol,incol)
     else:
       dict['sph_trmit']=" "
       dict['sph_pos']="%s,%s,%s"%tuple(jj) 
       if dict[self.names[i]].count(",")==2 :dict['sph_col']="color rgb <%s>"%dict[self.names[i]]
       else:dict['sph_col']=dict[self.names[i]]
       if self.names[i]=="NULL":temp=dict['strans'];dict['strans']=1.0
       if self.names[i]=="NULL1":temp1=dict['ss'];dict['ss']=0.1
       # add transparency modifications for addneighbor atoms... 
       if dir(self).count('newatoms'):
         if [s for s in self.newatoms if dot(s-self.atd[i],s-self.atd[i])==0]:
           if dict['addn_trmit']!=0:dict['sph_trmit']=" transmit %s  "%dict['addn_trmit']
           if dict['addn_color']:
             if dict['addn_color'].count(",")==2 :dict['sph_col']="color rgb <%s>"%dict['addn_color']
             else:dict['sph_col']=dict['addn_color']
       comp+=sphere%dict # this is a sphere...
       if self.names[i]=="NULL":dict['strans']=temp
       if self.names[i]=="NULL1":dict['ss']=temp1

   # draw the bonds if required...
   if bonds:
     bcol=bcol.split("-")
     maxb=maxb.split(',')
     nbtype=len(bonds.split(','))
     if nbtype != len(bcol) and len(bcol)==1:bcol=bcol*nbtype
     if nbtype != len(maxb) and len(maxb)==1:maxb=maxb*nbtype
     if len(bcol)!=nbtype:print "wrong # of bond colors";sys.exit()
     if len(maxb)!=nbtype:print "wrong # of bond cutoffs";sys.exit()
     for jj,j in enumerate(bonds.split(',')):
       for ii,i in enumerate(self.bonds(float(maxb[jj]),j)):
         comp+=cylinder%(i[0]-mc[0],i[1]-mc[1],i[2]-mc[2],i[3]-mc[0],i[4]-mc[1],i[5]-mc[2],bsize,bcol[jj],btrans,bamb)

   # draw the periphery of the unit cell...
   if drawcell: 
     corner=vec[0]+vec[1]+vec[2]
     orig=array([0.0,0.0,0.0])
     for i in vec:
       negc=corner-i
       comp+=cylinder%(orig[0]-mc[0],orig[1]-mc[1],orig[2]-mc[2],i[0]-mc[0],i[1]-mc[1],i[2]-mc[2],csize,ccol,ctrans,camb)
       comp+=cylinder%(corner[0]-mc[0],corner[1]-mc[1],corner[2]-mc[2],negc[0]-mc[0],negc[1]-mc[1],negc[2]-mc[2],csize,ccol,ctrans,camb)
       for j in vec: 
         if (i-j).any(): comp+=cylinder%(i[0]-mc[0],i[1]-mc[1],i[2]-mc[2],i[0]+j[0]-mc[0],i[1]+j[1]-mc[1],i[2]+j[2]-mc[2],csize,ccol,ctrans,camb)

   # print or call command...
   if call: 
       OUT=open('temp.pov','w'); print >>OUT,comp; OUT.close()
       if not rightsc.count("."):rightsc+='.'
       os.system('povray +O%s -Itemp.pov +P +d -W%s -H%s 2>/dev/null'%(out,w,w/eval(rightsc)))
   else:print comp


 def bonds(self,max,bond='all'):
   "find all bonds which are less than length max. bond is of format sp1-sp2 (ie. Co-O) or all." 

   def checkdist(dist):
     if dist<0.01:  # check for unphysical bonds...
       print "ERROR"
       print "at%s - %s"%(i,ii)
       print "at%s - %s"%(j,jj)
       print "bond with length <0.01... unphysical."; sys.exit()
   
   l1=[] ; l2=[] ; bondl=[] 
   if bond=='all': l1=self.atc; l2=self.atc 
   elif not bond.count('-'): print 'improper format for bond variable'; sys.exit()
   else:
     s1=bond.split('-')[0]; s2=bond.split('-')[1]
     if not self.specname.count(s1) or not self.specname.count(s2): print "improper bond label";sys.exit()
     for i,ii in enumerate(self.names): # get the atoms of each species in seperate lists...
       if ii==s1:l1.append(self.atc[i])
       if ii==s2:l2.append(self.atc[i])
       
   # calculate all bonds within max cuttoff...
   for i,ii in enumerate(l1):
     for j,jj in enumerate(l2):
       dist=sqrt(dot(ii-jj,ii-jj))
       if dist < max: # check to see if bond length is below max... 
         if bond=='all' or s1==s2: # if we seek all bonds or same species we dont want to double count...
           if i>j: bondl.append(list(ii)+list(jj)); checkdist(dist)
         else: 
           bondl.append(list(ii)+list(jj)); checkdist(dist) # otherwise we just take the bond...

   if [s for s in bondl if sqrt((s[0]-s[3])**2+(s[1]-s[4])**2+(s[2]-s[5])**2)<0.01]:
       print "bond with length <0.01... unphysical."; sys.exit()
   return bondl

 def cellcenter(self):
   "return the center of mass of the unit cell"
   vec=self.vec
   #return (vec[0]+vec[1]+vec[2]+  vec[0]+vec[1]+  vec[2]+vec[0]+  vec[2]+vec[1]+ vec[0]+vec[1]+vec[2])/8
   return (vec[0]+vec[1]+vec[2])/2
   

 def clustercenter(self):
   "center of the cluster of the unit cell"
   temp=array([0.0,0.0,0.0])
   for i in self.atc:temp+=i
   temp=temp/self.natoms
   self.atc=self.atc-temp
   self.atd=dot(self.atc,inv(self.vec))
   if self.coord=='c':self.at=self.atc*1
   elif self.coord=='d':self.at=self.atd*1

 def closestto(self,inp=0):
   "shift all atoms by translation vectors within a shell such that they are closest to atom inp"
   temp=self.atd[int(inp)]
   for i,ii in enumerate(self.atd):
     list=[]
     for j in range(-1,2):
       for k in range(-1,2):
         for l in range(-1,2):
           dist=dot(ii+array([j,k,l])-temp,ii+array([j,k,l])-temp)
           list.append((dist,(j,k,l)))
     list.sort()
     self.atd[i]=ii+array(list[0][1])
   self.atc=dot(self.atd,self.vec)
   if self.coord=='c':self.at=self.atc*1
   elif self.coord=='d':self.at=self.atd*1

 def addneighbor(self,bond,cut,trans=1):
   "add neighbors within a specific range to a given atom"
   try:name1=bond.split('-')[0];name2=bond.split('-')[-1]
   except:print "string format must be like Co-O";return
   if not [s for s in self.names if name1.count(s)] or not [s for s in self.names if name2.count(s)]:
     print "one of bond atoms (%s or %s) not present in structure... stopping."%(name1,name2); sys.exit()
   self.dist(3,cut,trans,trans,trans,output='')
   # select the appropriate bonds from something like this...
   # [1.83411 , "O28---Co10  | t1= 0 t2=-1 t3= 0 | d1=-0.00000 d2= 1.60540 d3=-0.88694 |  dist=1.83411" ]
   newlist=[ [s[1].split("---")[1].split()[0], [int(k) for k in s[1].split("|")[1].replace("="," ").split() if k[0]!='t']  ] 
           for s in self.distance if s[1].split("---")[0].count(name1) and s[1].split("---")[-1].split()[0].count(name2)]
   # remove atoms are already in unit cell (000) and any duplicates... 
   newlist=[[k.split()[0],[int(m) for m in k.split()[1:] ]] 
             for k in set(["%s %s %s %s"%(s[0],s[1][0],s[1][1],s[1][2]) 
             for s in newlist if not s[1][0]==s[1][1]==s[1][2]==0 ]) ]
   self.newatoms=[]
   for ii in newlist:
     atnum=int("".join([s for s in ii[0] if s.isdigit()]))-1
     #atname="".join([s for s in ii[0] if not s.isdigit()])
     #print atnum, atname, self.names[atnum],self.atd[atnum],ii[1]
     self.newatoms+=[self.atd[atnum]+array(ii[1])]

   # add new atoms to the cell and update relevant variables...
   self.atd=array(list(self.atd[:sum(self.specpoint[:(self.specname.index(name2)+1)])])+self.newatoms+ \
                  list(self.atd[sum(self.specpoint[:(self.specname.index(name2)+1)]):]))
   self.atc=dot(self.atd,self.vec)
   if self.coord=='d':self.at=self.atd
   elif self.coord=='c':self.at=self.atc
   self.specpoint[self.specname.index(name2)]+=len(self.newatoms)
   self.atpoint=[s+1 for s in range(self.nspecies) for j in range(self.specpoint[s])]
   self.natoms=sum(self.specpoint)
   self.names=[self.specname[s-1] for s in self.atpoint]

########################################################################################################
# below are function to convert from other code formats to vasp's poscar format...
def pwscfin_to_poscar(file):
  if type(file)==type([]):print "please give filename or file as string";sys.exit()
  if os.path.isfile(file): file=open(file).read()
  nat=int(re.compile("nat *= *").split(file)[-1].replace(',',' ').split()[0])
  #ntyp=re.compile("ntyp *= *").split(file)[-1].replace(',',' ').split()[0]
  vec="\n".join(file.split("CELL_PARAMETERS")[-1].lstrip().split("\n")[0:3] )
  try:atd=re.compile("ATOMIC_POSITIONS *\{ *crystal *\}").split(file)[-1].lstrip().split("\n")[0:nat]
  except:print "problem with atomic coord. be sure to give direct coordinates in pwscf.";sys.exit()
  names=[s.split()[0] for s in atd]
  # get distinct atomic species in order...
  types=[ ss for s,ss in enumerate(names) if ss!=names[s-1]]
  # get number of each specie...
  numbset=" ".join([str(len([ k for k in names if k==s])) for s in types])
  atd="\n".join([" ".join(ss.split()[1:])+" "+names[s] for s,ss in enumerate(atd)])
  poscar=("title \n1.0\n"+vec+'\n'+numbset+'\n'+'direct\n'+atd).split('\n')
  return poscar

def pwscfout_to_poscar(file):
  if type(file)==type([]):print "please give filename or file as string";sys.exit()
  if os.path.isfile(file): file=open(file).read()
  natom=int(file.split("number of atoms/cell      =")[-1].split()[0])
  vec=file.split("a(1) = (")[-1].split("reciprocal axes: (cart")[0].\
      replace("a(2) = (","").replace("a(3) = (","").replace(")","").rstrip()
  atc=[s.split() for s in file.split("positions (a_0 units)")[-1].lstrip().split("\n")[:natom]]
  names=[s[1] for s in atc]; atc="\n".join([" ".join(ss[6:9])+" "+names[s] for s,ss in enumerate(atc)])
  # get distinct atomic species in order...
  types=[ ss for s,ss in enumerate(names) if ss!=names[s-1]]
  # get number of each specie...
  numbset=" ".join([str(len([ k for k in names if k==s])) for s in types])
  poscar=("title \n1.0\n"+vec+'\n'+numbset+'\n'+'cartesian\n'+atc).split('\n')
  return poscar

def abinit_to_poscar(file):
  INP=file.split('\n')
  keyw={"acell":"","rprim":"","natom":"","typat":"","xred":""}
  # get line numbers...
  for i,ii in enumerate(INP):
    for k in keyw.keys(): 
      if ii.split('#')[0].count(k):keyw[k]=i
  # collect info...
  acell=INP[keyw['acell']].split("acell")[-1].split("*")[-1].split()[0]
  natom=INP[keyw['natom']].split("natom")[-1].split()[0]
  typat=INP[keyw['typat']].split("typat")[-1].split("#")[0]
  rprim=xred=typptr=''
  # convert typat...
  for i in sorted(list(set(typat.split()))):typptr+=str(typat.count(i))+" "
  # get vectors...
  for i in range(int(keyw['rprim']),len(INP)):
    tt=INP[i].replace("rprim",'').split("#")[0].strip() 
    if tt:rprim=rprim+'\n '+tt
    if len(rprim.split())==9:break
  # get atoms...
  for i in range(int(keyw['xred']),len(INP)):
    tt=INP[i].replace("xred",'').split("#")[0].strip() 
    if tt:xred=xred+'\n '+tt
    if len(xred.split())==3*int(natom):break
  poscar=("title\n"+acell+rprim+'\n'+typptr+'\n'+'direct'+xred).split('\n')
  return poscar 

##########################   Main Code   ###############################################################

if __name__=="__main__":
  from mysub import *

  fileh=[s for s in sys.argv[1:] if os.path.isfile(s)]
  if fileh:
    file=open(fileh[0])
    if len(fileh)==2:file2=fileh[1]
    else: file2=''
  else: file=sys.stdin

  # check to see if we need any file conversions...
  if   sys.argv.count('-pwscfin' ): file=pwscfin_to_poscar(file.read());  sys.argv.remove('-pwscfin')
  elif sys.argv.count('-pwscfout'): file=pwscfout_to_poscar(file.read()); sys.argv.remove('-pwscfout')
  elif sys.argv.count('-abinit'):   file=abinit_to_poscar(file.read()); sys.argv.remove('-abinit')
  else:file=file.readlines()

#  # check if file has any commands...
#  for i in file:
#    if i.count('=') and not i.count('#'): sys.argv=i.split()+sys.argv 
#  # check for breaks in data to avoid irrelevant data...
#  try:br=[s for s,ss in enumerate(file) if not ss.split()][0]+1
#  except:br=len(file)
#  file=[s for s in file[:br] if not s.count('=') and s.split() and not s[0]=='#']
    

  dict={'rescale':0.0,'resetscale':0.0,'dirtocar':'','pinfo':'',
        'dist':'','verb':0,'cut':2.5,'xyz':'','allincell':'','tol':8,
        'compxyz':'','label':'','rot':'','supa':'','ptart':'','pstut':'','pdmftlab':'',
        'shift':'','trans':1,'pov':'','call':'r','closestto':'',
        'execute':'','save':'','outofcell':'','addneighbor':'','dagger':'','strain':'','center':''}
  dict=update_dict(dict,sys.argv) # update the dictionary
  for i in dict.keys():  # execute the dictionary
    if type(dict[i])!=type(''): exec "%s=%s"%(i,dict[i])
    else: exec "%s='%s'"%(i,dict[i])

  pos=structure(file) # instantiate the structure class...

  if file2:        pos.combine(structure(file2))
  if label:        pos.label(label)
  if rescale:      pos.rescale(rescale)
  if resetscale:   pos.resetscale(resetscale)
  if shift:        pos.shift_origin(shift) 
  if dirtocar:     pos.dirtocar()
  if rot:          pos.rot(rot,dagger)
  if strain:       pos.axial_strain(strain)
  if supa:         pos.supa(supa)
  if compxyz:      pos.completecellxyz(compxyz)
  if allincell:    pos.allincell(tol)
  if addneighbor:  pos.addneighbor(addneighbor,cut)
  if outofcell:    pos.outofcell(outofcell)
  if center:       pos.clustercenter()
  if closestto:    pos.closestto(closestto)
  if execute:      exec(execute)

  if ptart:        pos.ptart()
  elif pstut:      pos.pstut()
  elif pdmftlab:   pos.pdmftlab()
  elif pinfo:      pos.pinfo();pos.prinfo()
  elif dist:       pos.dist(verb,cut,trans,trans,trans)
  elif xyz:        pos.xyz(xyz)
  elif pov:        pos.povray(pov.replace(":","="),call)
  elif save and fileh: 
    os.system('cp %s /tmp/'%fileh[0])
    OUT=open(fileh[0],'w')
    print >>OUT,pos.pposcar(1); OUT.close()
  else:          pos.pposcar()
