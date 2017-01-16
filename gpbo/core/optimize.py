# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.
from __future__ import print_function
xrange=range
import pickle
import scipy as sp
import os
import time
import logging
import copy
import pandas as pd
import gpbo

logger = logging.getLogger(__name__)

class optstate:
    def __init__(self):
        self.x = []
        self.ev = []
        self.y = []
        self.c = []
        self.C = 0
        self.n = 0
        self.Cfull=0.
        self.aqtime=[]
        self.aux=None
        self.localdone=False
        self.startlocal= None
        return
    
    def update(self,x,ev,y,c,taq):
        self.x.append(x)
        self.ev.append(copy.copy(ev))
        self.y.append(y)
        self.c.append(c)
        self.C +=c
        self.Cfull+=c+taq
        self.n+=1
        self.aqtime.append(taq)
        return 


class optimizer:
    def __init__(self,dirpath,name,aqpara,aqfn,stoppara,stopfn,reccpara,reccfn,ojf,ojfchar,choosefn,choosepara):
        self.dirpath = dirpath
        self.name = name
        self.setaq(aqpara,aqfn)
        self.setstopcon(stoppara,stopfn)
        self.setojf(ojf)
        self.setrecc(reccpara,reccfn)
        self.setchoose(choosepara,choosefn)
        self.ojfchar = ojfchar
        self.dx = ojfchar['dx']
        self.dev = ojfchar['dev']
        return
    
    def setaq(self,aqpara,aqfn):
        self.aqfn = aqfn
        self.aqpara = aqpara
        self.aqpersist = [None]*len(aqfn)
        return
    
    def setrecc(self,reccpara,reccfn):
        self.reccpara = reccpara
        self.reccfn = reccfn
        self.reccpersist = [None]*len(reccfn)
        return

    def setchoose(self,choosepara,choosefn):
        self.choosepara = choosepara
        self.choosefn =choosefn
        self.choosepersist = None
        return
    def setstopcon(self,stoppara,stopfn):
        self.stoppara = stoppara
        self.stopfn=stopfn
        return
    
    def setojf(self,ojf):
        self.ojf = ojf
        return
    
    def run(self):
        logger.info('startopt:')
        print( self.aqpara)
        self.stoppara['t0']=time.clock()
        lf = open(os.path.join(self.dirpath,self.name),'wb',0)
        lf.write(''.join(['n, ']+['x'+str(i)+', ' for i in xrange(self.dx)]+[i+', ' for i in self.aqpara[0]['ev'].keys()]+['y, c, ']+['rx'+str(i)+', ' for i in xrange(self.dx)]+['truey at xrecc, taq, tev, trc, realtime, aqauxdata'])+'\n')
        self.state = optstate()
        stepn=0
        while not self.stopfn(self.state,**self.stoppara):
            stepn+=1
            #print self.choosepara
            #print self.choosefn
            mode,self.choosepersist,chooseaux = self.choosefn(self.state,self.choosepersist,**self.choosepara)

            logger.info("---------------------\nstep {}:".format(stepn))
            
            t0 = time.clock()
            x,ev,self.aqpersist[mode],aqaux = self.aqfn[mode](self.state,self.aqpersist[mode],**self.aqpara[mode])
            t1 = time.clock()
            self.state.aux = aqaux
            logger.info("AQ returned {} : {}    aqtime: {}\nevaluate:".format(x,ev,t1-t0))
            
            y,c,ojaux  = self.ojf(x,**ev)
            t2 = time.clock()
            self.state.update(x,ev,y,c,t1-t0)
            logger.info("EV returned {} : {}     evaltime: {}".format(y,c,t2-t1))
            rx,self.reccpersist[mode],reaux = self.reccfn[mode](self.state,self.reccpersist[mode],**self.reccpara[mode])
            t3 = time.clock()
            
            if self.reccpara[mode]['check']:
                #logger.info("checkin {} : {}".format(rx,self.aqpara['ev']))
                checkpara=copy.copy(self.aqpara[mode]['ev'])
                checkpara['s']=1e-99
                checkpara['cheattrue']=True
                checky,checkc,checkojaux  = self.ojf(rx,**checkpara)
                #logger.info("checkout {} : {} : {}".format(checky,checkc,checkojaux))
            else:
                checky=sp.NaN

            
            logger.info("RC returned {}     recctime: {}\n".format(rx,t3-t2))
            logstr = ''.join([str(stepn)+', ']+[str(xi)+', ' for xi in x]+[str(evi[1])+', ' for evi in ev.items()]+[str(y)+', ']+[str(c)+', ']+[str(ri)+', ' for ri in rx]+[str(checky)+',']+[str(i)+', ' for i in [t1-t0,t2-t1,t3-t2]]+[time.strftime('%H:%M:%S  %d-%m-%y')])+','+''.join([str(k)+' '+str(aqaux[k]).replace(',',' ').replace('\n',';').replace('\r',';')+' ,' for k in aqaux.keys()])[:-1]+'\n'
            lf.write(logstr)

        import pickle
        obj = [self.reccpersist, self.aqpersist]
        pickle.dump(obj, open('dbout/persists.p', 'wb'))
        logger.info('endopt')

        return rx,reaux
    
def norlocalstopfn(optstate,**para):
    return nstopfn(optstate,**para) or localstopfn(optstate,**para)

def nstopfn(optstate,**para):
    return optstate.n >= para['nmax']

def localstopfn(optstate,**para):
    return optstate.localdone

def cstopfn(optstate,cmax = 1,includeaq=False):
    if not includeaq:
        logger.info('Used {} of {} evaluation budget.'.format(optstate.C,cmax))
        return optstate.C >= cmax
    else:
        logger.info('Used {} of {} evaluation budget.'.format(optstate.Cfull, cmax))
        return optstate.Cfull >= cmax

def totaltstopfn(optstate,**para):
    tused = sum(optstate.aqtime)+optstate.C

    if tused>=para['tmax']:
        logger.info('Time limit reached')
        return True
    else:
        hu=int(tused)/3600
        mu=(int(tused)%3600)/60
        su=int(tused)%60
        ht=int(para['tmax'])/3600
        mt=(int(para['tmax'])%3600)/60
        st=int(para['tmax'])%60
        logger.info('Used {}h {}m {}s of {}h {}m {}s budget \n of which {} acquisition {} evaluation'.format(hu,mu,su,ht,mt,st,(tused-optstate.C)/(1e-9+tused),optstate.C/(tused+1e-9)))
        return False


def search(optconfig):
    if not hasattr(optconfig,'fname'):
        optconfig.fname='traces.csv'
    multi=False
    if hasattr(optconfig,'multimode'):
        if optconfig.multimode:
            multi=True
    if not multi:
        O = optimizer(optconfig.path, optconfig.fname, [optconfig.aqpara], [optconfig.aqfn], optconfig.stoppara,
                                     optconfig.stopfn, [optconfig.reccpara], [optconfig.reccfn], optconfig.ojf,
                                     optconfig.ojfchar,gpbo.core.choosers.always0,dict())
    else:
        O = optimizer(optconfig.path, optconfig.fname, optconfig.aqpara, optconfig.aqfn, optconfig.stoppara,
                                     optconfig.stopfn, optconfig.reccpara, optconfig.reccfn, optconfig.ojf,
                                     optconfig.ojfchar,optconfig.chooser,optconfig.choosepara)

    return O.run()



def readoptdata(fname,includetaq=False):
    df = pd.DataFrame()

    with open(fname, 'r') as f:
        for line in f:
            df = pd.concat([df, pd.DataFrame([tuple(line.strip().split(','))])], ignore_index=True)

    j = 0
    for i in xrange(df.shape[1]):
        if not isinstance(df[i][0], str):
            df[i][0] = 'augdata{}'.format(j)
            j += 1
        else:
           df[i][0] = df[i][0].replace(' ', '')
    df.columns = df.iloc[0]
    df.drop(df.index[[0]], inplace=True)
    df.reset_index(inplace=True)
    l = len(df['c'])
    df['cacc'] = pd.Series(sp.empty(l), index=df.index)
    df['accE'] = pd.Series(sp.empty(l), index=df.index)
    df['accEA'] = pd.Series(sp.empty(l), index=df.index)
    for c in df.columns:
        try:
            df[c] = df[c].astype(float)  #
        except ValueError:
            pass


    #df['accEA'][0] = df.loc[0, ('c')]+df.loc[0, ('taq')]
    df.loc[0,('accEA')] = df.loc[0, ('c')]+df.loc[0, ('taq')]
    for i in xrange(1, l):
        df.loc[i,('accEA')] = df.loc[i - 1, 'accEA'] + df.loc[i, 'c']+df.loc[i, ('taq')]

    df.loc[0,('accE')] = df.loc[0, ('c')]
    for i in xrange(1, l):
        df.loc[i,('accE')] = df.loc[i - 1, 'accE'] + df.loc[i, 'c']


    if includetaq:
        df['cacc']=df['accEA']
    else:
        df['cacc']=df['accE']
    return df