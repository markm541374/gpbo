# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

import pickle
import scipy as sp
import os
import time
import logging
import copy
import pandas as pd

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
    def __init__(self,dirpath,name,aqpara,aqfn,stoppara,stopfn,reccpara,reccfn,ojf,ojfchar):
        print aqpara

        self.dirpath = dirpath
        self.name = name
        self.setaq(aqpara,aqfn)
        self.setstopcon(stoppara,stopfn)
        self.setojf(ojf)
        self.setrecc(reccpara,reccfn)
        
        self.ojfchar = ojfchar
        self.dx = ojfchar['dx']
        self.dev = ojfchar['dev']
        return
    
    def setaq(self,aqpara,aqfn):
        self.aqfn = aqfn
        self.aqpara = aqpara
        self.aqpersist = None
        return
    
    def setrecc(self,reccpara,reccfn):
        self.reccpara = reccpara
        self.reccfn = reccfn
        self.reccpersist = None
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
        
        lf = open(os.path.join(self.dirpath,self.name),'wb',0)
        lf.write(''.join(['n, ']+['x'+str(i)+', ' for i in xrange(self.dx)]+[i+', ' for i in self.aqpara['ev'].keys()]+['y, c, ']+['rx'+str(i)+', ' for i in xrange(self.dx)]+['truey at xrecc, taq, tev, trc, realtime, aqauxdata'])+'\n')
        self.state = optstate()
        stepn=0
        while not self.stopfn(self.state,**self.stoppara):
            stepn+=1
            logger.info("---------------------\nstep {}\naquisition:".format(stepn))
            
            t0 = time.time()
            x,ev,self.aqpersist,aqaux = self.aqfn(self.state,self.aqpersist,**self.aqpara)
            t1 = time.time()
            self.state.aux = aqaux
            logger.info("{} : {}    aqtime: {}\nevaluate:".format(x,ev,t1-t0))
            
            y,c,ojaux  = self.ojf(x,**ev)
            t2 = time.time()
            self.state.update(x,ev,y,c,t1-t0)
            logger.info("{} : {}     evaltime: {}\nreccomend:".format(y,c,t2-t1))
            rx,self.reccpersist,reaux = self.reccfn(self.state,self.reccpersist,**self.reccpara)
            t3 = time.time()
            
            if self.reccpara['check']:
                #logger.info("checkin {} : {}".format(rx,self.aqpara['ev']))
                checkpara=copy.copy(self.aqpara['ev'])
                checkpara['s']=1e-99
                checkpara['cheattrue']=True
                checky,checkc,checkojaux  = self.ojf(rx,**checkpara)
                #logger.info("checkout {} : {} : {}".format(checky,checkc,checkojaux))
            else:
                checky=sp.NaN

            
            logger.info("{}     recctime: {}\n".format(rx,t3-t2))
            logstr = ''.join([str(stepn)+', ']+[str(xi)+', ' for xi in x]+[str(evi[1])+', ' for evi in ev.items()]+[str(y)+', ']+[str(c)+', ']+[str(ri)+', ' for ri in rx]+[str(checky)+',']+[str(i)+', ' for i in [t1-t0,t2-t1,t3-t2]]+[time.strftime('%H:%M:%S  %d-%m-%y')])+','+''.join([str(k)+' '+str(aqaux[k]).replace(',',' ').replace('\n',';').replace('\r',';')+' ,' for k in aqaux.keys()])[:-1]+'\n'
            lf.write(logstr)

        import pickle
        obj = [self.reccpersist, self.aqpersist]
        pickle.dump(obj, open('dbout/persists.p', 'wb'))
        logger.info('endopt')

        return rx,reaux
    

def nstopfn(optstate,nmax = 1):
    return optstate.n >= nmax

def cstopfn(optstate,cmax = 1,includeaq=False):
    if not includeaq:
        logger.info('Used {} of {} evaluation budget.'.format(optstate.C,cmax))
        return optstate.C >= cmax
    else:
        logger.info('Used {} of {} evaluation budget.'.format(optstate.Cfull, cmax))
        return optstate.Cfull >= cmax



def search(optconfig):
    if not hasattr(optconfig,'fname'):
        optconfig.fname='traces.csv'
    O = optimizer(optconfig.path, optconfig.fname, optconfig.aqpara, optconfig.aqfn, optconfig.stoppara,
                                     optconfig.stopfn, optconfig.reccpara, optconfig.reccfn, optconfig.ojf,
                                     optconfig.ojfchar)

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

    for c in df.columns:
        try:
            df[c] = df[c].astype(float)  #
        except ValueError:
            pass
    if includetaq:
        df['cacc'][0] = df.loc[0, ('c')]+df.loc[0, ('taq')]
        for i in xrange(1, l):
            df['cacc'][i] = df.loc[i - 1, 'cacc'] + df.loc[i, 'c']+df.loc[i, ('taq')]
    else:
        df['cacc'][0] = df.loc[0, ('c')]
        for i in xrange(1, l):
            df['cacc'][i] = df.loc[i - 1, 'cacc'] + df.loc[i, 'c']
    return df