#!/usr/bin/env python 
import os,sys
class plot:
  "this class is use to create grace plots in python"
  
  # this is the directory where all scratch files will be stored...
  ddir=os.getenv('HOME') # use home if defined...
  if not ddir:ddir='./'  # else use curent directory...

  def __init__(self,mode=''): 
    """if run in interactive mode, the constructor creates the named pipe used 
    to communicate with grace (ie. linux command mkfifo) and starts grace."""
    if mode=='interactive': 
      self.mode='interactive'
      try:
        if not os.path.exists("%s/.gracepipe"%self.ddir):os.system("mkfifo %s/.gracepipe"%self.ddir)
        os.system("xmgrace -npipe %s/.gracepipe &"%self.ddir)
      except: print "error in creating named pipe %s/.gracepipe"%self.ddir;return
    else: self.mode=''
    self.fileset=[]   # a list with the number of sets in each file
    self.filenames=[] # a list with the name of each data file
    self.file_com=[]  # a list with a possible comment for each set in each file
    self.batch=''     # a compilation of grace commands which are issued
    self.error=0      # a variable to track errors
    self.funcs=0      # a variable to track number of defined function files
    self.gracefuncs=0 # a variable to track number of defined grace function files
    self.gfsave=''    # a variable to track the x-mesh for the grace function files
    # a list of the possible output devices...
    self.devlist=['PostScript','EPS','JPEG','PNG','SVG','null']
    self.devdict={'PostScript':'ps','EPS':'eps','JPEG':'jpg','PNG':'png','SVG':'svg','null':''}
    self.picname=''

  def data(self,dlist,space_delimited='',pcol=[]):
    """plot a data file. dlist is either a filename as a string or a list of
    strings. all columns in the file will be treated as data. space_delimited
    determines wether or not a empty row within a file should be treated as
    demarcating a new set (default=False not yet totally robust). pcol allows one 
    to plot only certain columns of files having multiple column sets."""
    if not type(dlist)==type([]):dlist=[s for s in dlist.split() if s]
    if pcol and not type(pcol)==type([]):pcol=[pcol]
    if len(pcol)>1 and len(dlist)!=len(pcol):print "problem with pcol(-u) option";self.error+=1;return
    elif len(pcol)==1:pcol=pcol*len(dlist)
    for ii,i in enumerate(dlist):
      if type(i) != type(''): print "must give string for filename: %s not correct"%i
      # check to see if filename has a local overide for the pcol command (ie. file::1:2:3,4:5 )
      if i.count("::"):lpcol=i.split("::")[-1]; i=i.split("::")[0]
      else:lpcol=''
      if not os.path.isfile(i): print "file %s does not exist!"%i;self.error+=1;return
      self.filenames.append(i)
      if not pcol and not lpcol:
        self.put('read nxy "%s"'%i)
        if not space_delimited:
          for j in open(i): # get the number of sets in the file...
            if not j.count('#') and j.split():self.fileset.append(len(j.split())-1);break        
        # new sets defined by blank lines in addition to extra columns...
        else: self.fileset.append( sum( 
              [int(k.split()[-1]) for k in 
              " ".join([str(len(s.split())-1) for s in open(i) ]).split("-1") if k.split() ]))
      elif pcol or lpcol:
        self.put('read block "%s"'%i); temp=0 
        if lpcol: list=lpcol.split(",")
        else:     list=pcol[ii].split(",")
        for kk in list:
            if not kk.count(":"):print "problem with pcol entry %s"%kk;self.error+=1;return
            for jj in kk.split(":")[1:]:self.put('block xy "%s:%s"'%(kk.split(":")[0],jj));temp+=1
        self.fileset.append(temp)

  def function(self,inp1,inp2=[]):
    """plot a list/array from python. if only a single list is given, it should
    contain x,y1,y2,etc. if two lists are given, the first should be the x data
    and the second is y1,y2,y3,etc."""
    ext=str(self.funcs) #str(random.randint(1000,9999))
    OUT=open('%s/.ppfuncfile'%self.ddir+ext,'w')
    if not len(inp2):
      for i in inp1:
        for j in i:print >>OUT,j,
        print >>OUT,""
      OUT.close()
    else:
      if len(inp1)!=len(inp2):print "improper input arrays: lengths differ.";return
      for i,ii in enumerate(inp1):
        print >>OUT,ii,
        for j in str(inp2[i]).split():print >>OUT,j,
        print >>OUT,""
      OUT.close()
    self.data('%s/.ppfuncfile'%self.ddir+ext);self.funcs+=1


  def put(self,inp): 
    "send command to the grace pipe and/or append to batch variable"
    self.batch+="%s\n"%inp #.replace("'",'"')
    if self.mode=="interactive":
      for i in inp.split('\n'):os.system("echo '%s' > %s/.gracepipe"%(inp,self.ddir)) 
      os.system("echo 'redraw' > %s/.gracepipe"%self.ddir)

  def save(self,filename=''): 
    "save the entire session as a grace file"
    if not filename:filename='%s/.ppgraceinput'%self.ddir
    self.put('saveall "%s"'%filename)

  def newsubplot(self): 
    "create a new subplot"
    if not dir(self).count('subplot'):self.subplot=0
    self.subplot+=1
    self.put("g%s on"%self.subplot)
    self.put("focus g%s"%self.subplot)
    # reset the fileset list
    self.fileset=[]


  def show(self,device=''):
    """if interactive mode is not being used, this can be used to call grace
    AFTER all of the desired commands have been enetered. also, this method
    will create a hardcopy if PostScript, EPS, JPEG, or PNG is specified in the
    constructor, even in interactive mode."""
    OUT=open("%s/.ppbatchfile"%self.ddir,'w'); print >>OUT,self.batch;OUT.close()
    if type(device)==type([]) and device:device=device[0] # allow input to be a list... take 1st element...
    if not device and self.mode!="interactive": 
      os.system("xmgrace -noask -nosafe -batch %s/.ppbatchfile"%self.ddir)
    elif self.devlist.count(device):
      if not device=='null':
        os.system("xmgrace -nosafe -noask -hardcopy -hdevice %s -batch %s/.ppbatchfile"%(device,self.ddir))
        if self.picname:os.system("mv Untitled.%s %s.%s"%(self.devdict[device],self.picname,self.devdict[device]))
    else: print "must enter a proper device: \n%s"%(" ".join(self.devlist))


  def filestyle(self,inp=''):
    """define set attributes for all sets in a file. for example, when plotting
    two files with many sets, one want all sets from each file colored the
    same.  the method should be given an input string such as "lc=[1,2]" which
    means color all the lines in the first file color 1, etc (ie. for a plot
    made from 2 files). defaults are provided for lc,sc,sfc if nothing given.
    s=auto will give an automatic symbol for all sets in a given file. leg=auto
    will create a legend with the filename for each set within a file.  """
    # process the input string...
    if type(inp)==type(''):inp=inp.split()
    # provide defaults for line colors, symbol colors, and symbol fill colors if not given...
    # also provide filename legend and automatic symbols if leg=auto, s=auto...
    color='';legend='';symbol=''
    for i,ii in enumerate(self.fileset): 
      for j in range(ii):
        if i<=14:  color+="%s,"%(i+1)
        else :     color+="%s,"%(str(i)[-1])
        if i<10:   symbol+="%s,"%(i+1)
        else :     symbol+="%s,"%(str(i)[-1])
        legend+="%s,"%self.filenames[i]
    if inp.count("leg=auto"): 
      self.quickin("ll=0.8,1 lbf=0",use_default=True)
      for p,pp in enumerate(legend.split(",")):
        if pp!=legend.split(",")[p-1]: self.put('s%s legend "%s"'%(p,pp))
    if inp.count("s=auto"):   self.quickin("s=[%s]"%symbol[:-1],sum(self.fileset),use_default=True)
    inp=[s for s in inp if not s.count("=auto") ] # remove auto entry...
    if not [s for s in inp if s.count("lc=") or s.count("sc=") or s.count("sfc=")]:  
      self.quickin("lc=[%s] sc=[%s] sfc=[%s]"%tuple([color[:-1]]*3),sum(self.fileset),use_default=True)

    # process general input...
    nf=''
    for i in inp:
      key,val=i.split('=')[0],"=".join(i.split('=')[1:]); temp=''
      if val.count('['):val=val.replace("[","").replace("]","").split(',')
      else:val=[val]
      if len(val)==1:    # if only one option, store for later...
        #self.quickin("%s=%s"%(key,val[0]),use_default=True)
        nf+=" %s=%s "%(key,val[0])
      elif len(val)!=len(self.fileset): # check for proper number of items given...
        print self.fileset
        print "number of options for %s dont match number of files"%key;self.error+=1;return
      else:              # construct the appropriate input array...
        if key!='leg':
          for i,ii in enumerate(self.fileset): 
            for j in range(ii):
              temp+="%s,"%(val[i])
          self.quickin("%s=[%s]"%(key,temp[:-1]),use_default=True)
        else:
          for p,pp in enumerate(val):
            if pp!=val[p-1]: self.put('s%s legend "%s"'%(sum(self.fileset[:p]),pp.replace("__"," ")))
    self.quickin(nf) # now execute all single options together...

        
  def color_grade(self,filestyle=False):
    # make the line color grade from color to color...
    nset=sum(self.fileset)
    if filestyle:nset=len(self.fileset) # grade line color by number of files instead...
    # two different color schemes...
    if not (nset-3)%2:div1=div2=(nset-3)/2
    else:div1=(nset-3)/2; div2=div1+1
    col=[[0,255,0]]
    col+=[[0,255,255*s/(div1+1)] for s in range(1,div1+1)]
    col+=[[0,255,255]]
    col+=[[0,255-255*s/(div2+1),255] for s in range(1,div2+1)]
    col+=[[0,0,255]]
    
    col=[[255,0,0]]
    col+=[[255,255*s/(div1+1),0] for s in range(1,div1+1)]
    col+=[[255,255,0]]
    col+=[[255-255*s/(div2+1),255,0] for s in range(1,div2+1)]
    col+=[[0,255,0]]
    col.reverse()
    for i in range(nset): 
      self.put('map color %s to (%s, %s, %s), "c%s"'%(i+100,col[i][0],col[i][1],col[i][2],i))
    # apply the colors using quickin method...
    temp="lc=["
    if not filestyle:
      for i in range(nset): temp+="%s,"%(i+100)
    else:
      for i,ii in enumerate(self.fileset): 
        for j in range(ii): temp+="%s,"%(i+100)
    temp=temp[:-1]+"]"
    self.quickin(temp)

  def auto_symbol(self):
    """this method automatically gives a unique symbol to each set. however, there are only ten symbols
    defined in xmgrace so the symbols will repeat after ten sets."""
    self.quickin("s=[%s]"%",".join([str(int(str(s)[0])+1) for s in range(sum(self.fileset))]))



  def auto_leg(self):
    """this method creates a legend where each set is labeled by it's filename and set number.""" 
    legend=''
    for i,ii in enumerate(self.fileset): 
      for j in range(ii): legend+="%s__set%s,"%(self.filenames[i],j+1)
    self.quickin("ll=0.8,1 lbf=0 leg=[%s]"%legend,use_default=True)

  def gracefunc(self,inp):
    """this method plots a function using grace's built-in capabilities. all we need to do is provide the x-mesh.
    inp should be given as a string or list of xmin,xmax,xinc,func or simply "xmin,xmax,xinc func1 func2 etc" or
    even "xmin,xmax,xinc func1  xmin,xmax,xinc func2 func3 etc". convention should be clear. if xinc is chosen
    to be less than zero, it is assumed to be the number of points and not the increment. enter "help" to print
    selected grace functions.
    """
    if inp.count('help'): print \
    """selected grace functions. see manual for more details. \
    \nabs(x) \nacos(x) \nacosh(x) \nasin(x) \nasinh(x) \natan(x) \natanh(x) \
    \nceil(x) \ncos(x) \ncosh(x) \nexp(x) \nfac(n) \nfloor(x) \nirand(n) \nln(x) \
    \nlog10(x) \nlog2(x) \nmod(x,y) \npi \nrand \nrint(x) \nrsum(x) \nsin(x) \
    \nsinh(x) \nsqr(x) \nsqrt(x) \ntan(x) \ntanh(x) \nerf(x) \nerfc(x) \nfdtrc(x) \
    \nfdtri(x) \ngdtr(a,b,x) \ngdtrc(a,b,x) \nndtr(x) \nndtri(x) \nnorm(x) \
    \npdtr(k,m) \npdtrc(k,m) \npdtri(k,y) \nstdtr(k,t) \nstdtri(k,p) \nai(x) \
    \nbi(x) \nbeta(x) \nchi(x) \nci(x) \ndawsn(x) \nellie(phi, \nellik(phi, \
    \nellpe(m) \nellpk(m) \nexpn(n,x) \nfresnlc(x) \nfresnls(x) \ngamma(x) \
    \nhyp2f1(a,b,c,x) \nhyperg(a,b,x) \ni0e(x) \ni1e(x) \nigam(a,x) \nigamc(a,x) \
    \nigami(a,p) \nincbet(a,b,x) \nincbi(a,b,y) \niv(v,x) \njv(v,x) \nk0e(x) \
    \nk1e(x) \nkn(n,x) \nlbeta(x) \nlgamma(x) \npsi(x) \nrgamma(x) \nshi(x) \nsi(x) \
    \nspence(x) \nstruve(v,x) \nyv(v,x) \nzetac(x) \nselected grace functions. see manual for more details.  
    """; self.error+=1;return
    if type(inp)!=type([]):inp=inp.split()
    from math import pi
    for i in inp:
      i=i.split(",");
      if len(i)>=3 or self.gfsave=='':
        self.gfsave=self.gracefuncs # only print out a new xmesh file if required (gfsave tracks this).
        if len(i)>=3:xmin=float(eval(i[0]));xmax=float(eval(i[1]));xinc=float(eval(i[2]))
        else:        xmin=-10;xmax=10;xinc=0.01
        if xinc<0:xinc=(xmax-xmin)/xinc*-1
        if xinc>(xmax-xmin): print "problem with xinc for gracefunc";self.error+=1;return
        xval=[(str(xmin+s*xinc)+" ")*2 for s in range(int(round((xmax-xmin)/xinc))) if (xmin+s*xinc)<xmax ]+["%s %s \n"%(xmax,xmax)]
        OUT=open("%s/.ppgracefunc%s"%(self.ddir,self.gfsave),'w');OUT.write("\n".join(xval));del xval;OUT.close()
      if len(i)==1 or len(i)==4:
        func=i[-1].replace("[","(").replace("]",")").replace("exp","EXP") # x is replaced below... so take care of exp.
        self.data("%s/.ppgracefunc%s"%(self.ddir,self.gfsave))
        self.put("s%s.y=%s; autoscale"%(self.gracefuncs,func.replace("x","s%s.x"%self.gracefuncs)))
        self.gracefuncs+=1

  def operate_on_set(self,input):
    """operate on the datasets and modify the x and/or y values. x4..6:*2.5"""
    if not type(input)==type([]):input=input.split()
    for inp in input:
      set=inp[0]
      if not set=='x' and not set=='y':print "must give x or y as first character in string";self.error+=1;return
      if not inp.count(":"):print ": must delimit set and operation";self.error+=1;return
      setnum=[int(s) for s in inp.split(":")[0][1:].split("..") if s];setnum.sort()
      if len(setnum)==2:setnum=range(setnum[0],setnum[1]+1)
      elif len(setnum)==0:setnum=range(sum(self.fileset))
      elif len(setnum)>=3:print "must give a range of sets, 1 set, or no sets means all sets";self.error+=1;return
      oper=inp.split(":")[-1]
      for i in setnum:
        if i<0:i+=sum(self.fileset)
        self.put("S%s.%s=S%s.%s%s"%(str(i),set,str(i),set,oper))
      if not self.batch.count('xmax') and not self.batch.count('ymax'):self.put("autoscale")

  def contour(self,files):
    """this method allows for a primitive contour map. three columns must be specified for each set
       with the third column being the point color. only three columns per file. mulitple sets can be given 
       with delimiting sets by a (horizontal) space."""
    # create a color map...
    div1=div2=50
    col=[[255,0,0]]
    col+=[[255,255*s/(div1+1),0] for s in range(1,div1+1)]
    col+=[[255,255,0]]
    col+=[[255-255*s/(div2+1),255,0] for s in range(1,div2+1)]
    col+=[[0,255,0]]
    col.reverse()
    #col=[[int(255*x/100.),int(255*(1-x/100.)),0] for x in range(101)]
    for i in range(len(col)): 
      self.put('map color %s to (%s, %s, %s), "c%s"'%(i+100,col[i][0],col[i][1],col[i][2],i+100))
    if type(files)!=type([]):files=[files]
    for i in files:
      # find number of sets in file
      inp=[s for s,ss in enumerate(open(i).readlines()) if not ss.split() or ss.lstrip()[0]=='#']
      nsets=sum([1 for s,ss in enumerate(inp) if inp[s-1]!=inp[s]-1])+1
      self.put('read xycolor "%s"'%i)
      self.fileset.append(nsets)
      self.filenames.append(i) 
    #self.quickin("sfp=1 sp=0")
    self.quickin("sfp=1 sls=0 s=dia ")

  def postprocess(self):
    """this method does any postprocessing"""
    if not (self.batch.count("tick major") or self.batch.count("tick minor")):self.put("autoticks")

  def comment_to_legend(self):
    """takes comments from input file and makes these the legend label. comment must be
       of the form '# comment1 comment2 ...', and this must only be on the line before the set starts.
       this will likely be used in a space delimited file, so remeber to set space_delimited flag (ie. -s).  """
    for i,ii in enumerate(self.filenames):
      self.file_com.append(os.popen(r"grep \# %s"%ii).read().replace(r'#','').split())
      while len(self.file_com[-1])<self.fileset[i]:self.file_com[-1]+=['']
    for i,ii in enumerate(self.fileset):
      for j in range(ii): 
        comm=self.file_com[i][j].replace("__"," ").strip()
        if comm: self.put('s%s legend "%s"'%(sum(self.fileset[:i])+j,comm))


  def quickin(self,modparam,nset=0,use_default=True,debug=''):
    """this method allows most grace commands with a shortcut syntax. a
    dictionary is used to store the shortcut names as keys, and the values are
    a two item list with the first item being the grace command and the second
    item the deafault value (optiotnal) of the grace command. in the case of
    commands corresonding to data set quantities (ie. line color), the value of
    the grace command is given as a list. this method takes in a string such as
    "t=plot1 xl=eV" and will look through the dictionary to find the proper
    grace command. any text will automatically be quoted. set quantities are
    specified as "s=[1,2,3]", which correspond to a plot with three sets, and
    assigns symbol type one for set one, etc.  one may also specify plot
    defaults by manually modifying the dictionary below."""

    gdict=self.gracedict # define the dictionary with all the grace commands...

    if not nset:nset=sum(self.fileset)
    if type(modparam)!=type([]):
      temp=[]
      for i in modparam.split():
        if i.count("="):temp.append(i)
        else:temp[-1]+=" "+i
      modparam=[s for s in temp]
    if not [s for s in modparam if s.count("=")] and modparam:print "all quickin entries must have =";self.error+=1;return
    # convert to list of list by splitting on first =
    modparam=[[s.split("=")[0],"=".join(s.split("=")[1:])] for s in modparam]

    #################################
    # handle shortcuts...
    # this takes in a string as str=hi,0,0 instead of strzl=0,0 strz=hi
    tt=[s for s in modparam if s[0]=="str"]
    if tt:
      for i in tt:
        # changed line below so string could contain commas...
        modparam.insert(modparam.index(i)+1,["strz",i[-1].replace(","+",".join(i[-1].split(',')[-2:]),"")])
        modparam.insert(modparam.index(i)+2,["strzl",",".join(i[-1].split(',')[-2:])])
        modparam.pop(modparam.index(i))
    # customizations for linel and strzl...
    if [s for s in modparam if s[0]=="linel" or s[0]=="strzl"]:
      import re
      # this allows one to use xmin/etc in linel if defined...
      p=re.compile("(xmin|xmax|ymin|ymax)")
      tt=[s for s in modparam if s[0]=="linel" or s[0]=="strzl" and p.findall(s[-1])]
      if tt:
        for i in tt:
          mm=p.findall(i[-1])
          ind=modparam.index(i)
          for j in mm:
            vv=[s for s in modparam if s[0]==j]
            if not vv: print "cannot sub %s in linel because it is not defined"%j;sys.exit()
            modparam[ind]=[i[0],modparam[ind][-1].replace(j,vv[-1][-1])]

      # this allows one to use arithmetic in defining numbers...
      #p=re.compile("(\*|\-|\+|\/)")
      tt=[s for s in modparam if (s[0]=="linel" or s[0]=="strzl")] #and p.findall(s[-1])]
      if tt:
        for i in tt:
            modparam[modparam.index(i)]=[i[0],",".join([ str(eval(s)) for s in   i[1].split(",")])]
    #################################

    # check to be sure all items are in dictionary...
    bad=[s[0] for s in modparam if not gdict.keys().count(s[0])]
    if bad:
      for i in sorted(gdict.keys()):
        print "    %s%s<-->    %s%s  default=%s"%(i," "*(10-len(i)),gdict[i][0]," "*(30-len(gdict[i][0])),gdict[i][1])
      for i in bad: print "there is no quickinp keyword '%s'. see above printout for proper shortcut.\n"%i
      self.error+=1;return
    
    # define a list of the specified default values... unless modparam overrides...
    if use_default: default=[[s,gdict[s][1]] for s in gdict.keys() if gdict[s][1] and not [k[0] for k in modparam].count(s)]
    else:default=[]
    
    # construct the batchfile... be sure to see if input should be converted to an array
    header=""
    for i in default+modparam: 
      # if a scalar is given for an array variable, convert it
      if type(gdict[i[0]][1]) is type([]) and not i[1].count("["): i[1]=[i[1]]*nset
      # if an array is given, be sure it is correct size
      if type(gdict[i[0]][1]) is type([]) and i[1].count("["):
        i[1]=i[1].replace('[','').replace(']','').replace(',',' ').split()
        if len(i[1])!=nset: print "improper array size given for input variable %s"%i[0];self.error+=1;return 
      kw=i[0]; gstr=gdict[kw][0]; gval=i[1]; quo=gdict[kw][2]  # rename these for easier syntax below...
      # only print dict entries which have been given a value
      if type(gval) is not type([]): 
        # see if we need to replace any shortcuts for colors, symbols, etc...
        if self.styledict.keys().count(gval):gval=self.styledict[gval]
        # see if quotes are needed...
        if not quo: header=header+"%s %s \n"%(gdict[kw][0],gval) 
        else: header=header+'%s "%s" \n'%(gstr,gval) 
      else:
        for j in range(len(gval)): 
          # see if we need to replace any shortcuts for colors, symbols, etc...
          if self.styledict.keys().count(gval[j]):gval[j]=self.styledict[gval[j]]
          # see if quotes are needed...
          if not quo:
              header=header+"s%s %s %s \n"%(j,gstr,gval[j])
          else:
              header=header+'s%s %s "%s" \n'%(j,gstr,gval[j])
    # replace .sp. or __ by a space to avoid quoting at command line...
    header=header.replace("__"," ").replace(".sp."," ")
    # the method is done... return the batch or send it to grace...
    if debug:return header
    else: self.put(header)

  gracedict={
# shortcut            grace command                    value  quote?  
  "ps"           :[  "page size  "                   ,   ""   ,0],  # laptop "800,450"  full page 1000,600
  "pr"           :[  "page resize  "                 ,   ""   ,0],   
  "ft"           :[  "frame type "                   ,   ""   ,0],  # 0 1
  "fl"           :[  "frame linestyle "              ,   ""   ,0],  
  "flw"          :[  "frame linewidth "              ,   ""   ,0],  
  "fc"           :[  "frame color "                  ,   ""   ,0],  
  "fbc"          :[  "frame background color "       ,   ""   ,0],  
  "fp"           :[  "frame pattern "                ,   ""   ,0],  
  "fbp"          :[  "frame background pattern "     ,   ""   ,0],  
  "xmin"         :[  "world xmin "                   ,   ""   ,0],
  "xmax"         :[  "world xmax "                   ,   ""   ,0],
  "ymin"         :[  "world ymin "                   ,   ""   ,0],
  "ymax"         :[  "world ymax "                   ,   ""   ,0],
  "v"            :[  "view "                         ,   ""   ,0],  # 0.15,0.15,1.15,0.85
  "t"            :[  "title "                        ,   ""   ,1],
  "tf"           :[  "title font "                   ,   ""   ,0],
  "ts"           :[  "title size "                   ,   "2.5",0],
  "tc"           :[  "title color "                  ,   ""   ,0],
  "st"           :[  "subtitle    "                  ,   ""   ,1],
  "sts"          :[  "subtitle size "                ,   ""   ,0],
  "xa"           :[  "xaxis  "                       ,   ""   ,0], # on off
  "ax"           :[  "autoscale xaxes  "             ,   ""   ,0], 
  "xl"           :[  "xaxis  label "                 ,   ""   ,1],
  "xls"          :[  "xaxis  label char size "       ,   "2"  ,0],
  "xlf"          :[  "xaxis  label font "            ,   ""   ,0],
  "xt"           :[  "xaxis  tick major "            ,   ""   ,0],
  "xtm"          :[  "xaxis  tick minor ticks "      ,   ""   ,0],
  "xtl"          :[  "xaxis  ticklabel "             ,   ""   ,0], # on off
  "xtls"         :[  "xaxis  ticklabel char size "   ,   "1.5",0],
  "xtlf"         :[  "xaxis  ticklabel font "        ,   ""   ,0],
  "xtlc"         :[  "xaxis  ticklabel color "       ,   ""   ,0],
  "xtlo"         :[  "xaxis  ticklabel offset spec \
                      \n  xaxis ticklabel offset  "  ,   ""   ,0], # 0.0,0.0
  "xtp"          :[  "xaxis  tick place"             ,   ""   ,0], # both normal opposite
  "ya"           :[  "yaxis  "                       ,   ""   ,0], # on off
  "ay"           :[  "autoscale yaxes  "             ,   ""   ,0], 
  "yl"           :[  "yaxis  label "                 ,   ""   ,1],
  "yls"          :[  "yaxis  label char size "       ,   "2"  ,0],
  "ylf"          :[  "yaxis  label font "            ,   ""   ,0],
  "yt"           :[  "yaxis  tick major "            ,   ""   ,0],
  "ytm"          :[  "yaxis  tick minor ticks "      ,   ""   ,0],
  "ytl"          :[  "yaxis  ticklabel "             ,   ""   ,0], # on off
  "ytls"         :[  "yaxis  ticklabel char size "   ,   "1.5",0],
  "ytlf"         :[  "yaxis  ticklabel font "        ,   ""   ,0],
  "ytlc"         :[  "yaxis  ticklabel color "       ,   ""   ,0],
  "ytlo"         :[  "yaxis  ticklabel offset spec \
                      \n  yaxis ticklabel offset  "  ,   ""   ,0], # 0.0,0.0
  "ytp"          :[  "yaxis  tick place"             ,   ""   ,0], # both normal opposite
  "ll"           :[  "legend "                       ,   ""   ,0], #"0.85,0.8"
  "llt"          :[  "legend loctype "               ,   ""   ,0], # view or world
  "lcs"          :[  "legend char size "             ,   ""   ,0],
  "lbf"          :[  "legend box fill pattern "      ,   ""   ,0],
  "lbl"          :[  "legend box linestyle"          ,   "0"  ,0], # 1 0 
  "lvg"          :[  "legend vgap "                  ,   ""   ,0],
  "lhg"          :[  "legend hgap "                  ,   ""   ,0],
  "ys"           :[  "yaxes scale  "                 ,   ""   ,0], # Normal Logarithmic
  "xs"           :[  "xaxes scale  "                 ,   ""   ,0], # Normal Logarithmic
  "strz"         :[  "with string \
                      \n  string on \
                      \n  string def"                ,   ""   ,1], # string definition
  "strzl"        :[  "string "                       ,   ""   ,0], # 0.85,0.8 ! need to use this !
  "strc"         :[  "string color "                 ,   ""   ,0],  
  "strr"         :[  "string rot "                   ,   ""   ,0],  
  "strf"         :[  "string font "                  ,   ""   ,0],  
  "strj"         :[  "string just "                  ,   ""   ,0],  
  "strcs"        :[  "string char size "             ,   "2"   ,0],  
  "strlt"        :[  "string loctype "               ,   "world"   ,0], # view or world 
  "linel"        :[  "with line \
                      \n  line on \
                      \n  line def \
                      \n  line"                      ,   ""   ,0], # user drawn line (ie. NOT data)
  "linec"        :[  "line color        "            ,   ""   ,0],  
  "linelt"       :[  "line loctype      "            ,   "world"   ,0], # view or world 
  "linew"        :[  "line linewidth    "            ,   ""   ,0],  
  "lines"        :[  "line linestyle    "            ,   ""   ,0],  
  "linea"        :[  "line arrow        "            ,   ""   ,0], # 0 
  "lineat"       :[  "line arrow type   "            ,   ""   ,0], # 0 
  "linealen"     :[  "line arrow length "            ,   ""   ,0], # 1 
  "linealay"     :[  "line arrow layout "            ,   ""   ,0], # 1,1 
  "s"            :[  "symbol "                       ,   []   ,0],
  "ss"           :[  "symbol size "                  ,   []   ,0],
  "sc"           :[  "symbol color "                 ,   []   ,0],
  "sfc"          :[  "symbol fill color "            ,   []   ,0],
  "sfp"          :[  "symbol fill pattern "          ,   []   ,0],
  "sp"           :[  "symbol pattern "               ,   []   ,0],
  "sls"          :[  "symbol linestyle "             ,   []   ,0],
  "slw"          :[  "symbol linewidth "             ,   []   ,0],
  "sch"          :[  "symbol char "                  ,   []   ,0],
  "schf"         :[  "symbol char font  "            ,   []   ,0],
  "leg"          :[  "legend "                       ,   []   ,1],
  "l"            :[  "line type "                    ,   []   ,0],
  "ls"           :[  "line linestyle "               ,   []   ,0],
  "lw"           :[  "line linewidth "               ,   []   ,0],
  "lc"           :[  "line color "                   ,   []   ,0]}

  # this is a dictionary which contains abbreviated names for all the colors/symbols numbers...
  # be sure not to have any repeats!
  styledict={
  # abreviated short color names...
  "wh":0,"bl":1,"re":2,"gr":3,"bu":4,"ye":5,"br":6,"gy":7,"vi":8,"cy":9,"ma":10,"or":11,"in":12,"pu":13,"tq":14,"grd":15,
  # abreviated full color names...
  "white":0,"blue":1,"red":2,"green":3,"blue":4,"yellow":5,"brown":6,"grey":7,"violet":8,"cyan":9,"magenta":10,
  "orange":11,"indigo":12,"purple":13,"turquoise":14,"greendark":15,
  # abreviated symbol names...
  "o":1,"sq":2,"dia":3,"tu":4,"tl":5,"td":6,"tr":7,"+":8,"x":9,"*":10,
  # these shortcuts allow one to avoid quoting an entry... and = because this defines shortcut command above...
  "__":" ",".sp.":" ",".eq.":"=",
  }
  

#########################################################################
# the above tools may directly be used in shell mode.
# simply give a data file on the command line or pipe one. 
# additionally, one may both pipe AND give a data file.
# by default, the pipe will come last in the list of files. this may be 
# changed by instead giving the word "pipe" on the command line.
# any command line item with "=" will be given to the quickin method.
# i only support one plot per page in shell mode because i find
# that doing multiple plots per page on the command line is cluttered.
#########################################################################
if __name__=="__main__":
  mp=plot()       # instantiate the plot...
  # check for file with input commands...
  ifile=[sys.argv[s+1] for s,ss in enumerate(sys.argv) if ss.count("-i") and os.path.exists(sys.argv[s+1])]
  for i in ifile:
    ii=sys.argv.index(i)-1
    sys.argv=sys.argv[0:1]+sys.argv[1:ii]+" ".join([ m for m in [k.strip() for k in open(i) if k.strip()] if m[0]!="#"]).split()+sys.argv[ii:]  
  # get all the quickin commands...
  commands=[s for s in sys.argv[1:] if s.count('=') and s[0]!='+']
  # get the input data files and allow manual listing of pipe
  file=[s for s in sys.argv[1:] if (os.path.isfile(s.split("::")[0]) or s=='pipe') and not ifile.count(s) ]
  # automatically detect any pipes...
  if not sys.stdin.isatty() and not file.count('pipe'):file+=['pipe'];apipe=True
  # get any hardcopy options...
  hardcopy=[s for s in sys.argv if mp.devlist.count(s)]+['']
  # check for other options
  options=[s for s in sys.argv if ['-fs','-cg','-cgfs','-s','-u','-op','-save','-con','-com','-i'].count(s)]
  # check for plotting specific columns...
  pcol=[sys.argv[s+1] for s,ss in enumerate(sys.argv) if ss.count("-u") and sys.argv[s+1].count(":")]
  if len(pcol)==1:pcol=pcol*len(file)
  else:
    temp=['']*len(file)
    for i,ii in enumerate(file):
      for j in pcol:
        if sys.argv.index(ii)>sys.argv.index(j):temp[i]=j
    pcol=[s for s in temp if s]
  # check for operations on data sets...
  oper=[sys.argv[s+1] for s,ss in enumerate(sys.argv) if ss.count("-op") and sys.argv[s+1].count(":")]
  # check for contour files...
  if sys.argv.count("-con"):
    try: confile=[sys.argv[s+1] for s,ss in enumerate(sys.argv) if ss.count("-con") and  
             (os.path.isfile(sys.argv[s+1]) or sys.argv[s+1]=='pipe')]
    except:confile=[]
    if dir().count("apipe") and sys.argv.count("-con") and not confile:confile=["pipe"]
    file=[s for s in file if not confile.count(s)] # remove contour files from file list...
  else: confile=[]
  # check for raw xmgrace commands...
  put=[s for s in sys.argv if s[0:2]=='++' ]
  # check for raw python commands...
  pyput=[s for s in sys.argv if s[0]=='+' and not put.count(s)]
  for i in pyput:exec(i[1:])
  # check for user provided grace functions to plot...
  gracefunc=[s for s in sys.argv if s[0]=="@"]
  # check to be sure no extraneous info was given... 
  bad=set(sys.argv[1:])-set(commands+file+hardcopy+options+put+pcol+gracefunc+oper+confile+ifile+pyput)
  if bad:print "the following options are not recognized:\n"," ".join(list(bad));sys.exit()
  # deal with piped data...
  if file.count('pipe') or confile.count('pipe'):
    OUT=open('%s/.ppdatafile'%mp.ddir,'w')
    for i in sys.stdin: print >>OUT,i.strip()
    OUT.close() 
    if file.count('pipe'):      file.append('%s/.ppdatafile'%mp.ddir) 
    elif confile.count('pipe'): confile.append('%s/.ppdatafile'%mp.ddir)
    file=[s for s in file if s!='pipe']
    confile=[s for s in confile if s!='pipe']

  # add any user input grace functions...
  if gracefunc:mp.gracefunc([s[1:] for s in gracefunc])
  
  # now we can add the data to the plot...  
  if options.count("-s"):      mp.data(file,space_delimited=True) 
  elif options.count("-u"):    mp.data(file,pcol=pcol) 
  else:                        mp.data(file)
  if options.count("-con"):    mp.contour(confile) 

  # check for automatic symbol and legend option... filestyle does seperately...
  if commands.count("s=auto") and not options.count('-fs'): 
    mp.auto_symbol(); del(commands[commands.index("s=auto")]) 
  if commands.count("leg=auto") and not options.count('-fs'): 
    mp.auto_leg(); del(commands[commands.index("leg=auto")]) 

  # check for filestyle option or the regular case...
  if options.count('-fs'):mp.filestyle(commands)   # this applies quickin commands per file instead of per set...
  else: mp.quickin(commands,use_default=True)     # execute the shortcut commands... 

  # the color grade and filestyle colorgrade options...
  if options.count("-cg"):mp.color_grade()
  if options.count("-cgfs"):mp.color_grade(True)

  # convert comments to legend if requested...
  if options.count("-com"):mp.comment_to_legend()

  # execute any user defined operations...
  if oper:mp.operate_on_set(oper)

  # execute any raw xmgrace commands if given and then postprocess...
  for i in put:mp.put(i[2:])                
  mp.postprocess()

  # save gracefile if requested...
  if options.count("-save"):mp.save()

  # print out the resutling batch file that xmgrace will load and show plot...
  if not mp.error: print mp.batch; mp.show(hardcopy[0])     
  else: print "ERROR: problem with plot... see above."
