# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.
from __future__ import print_function
xrange=range
import scipy as sp
from scipy.optimize import minimize
import os
import sys
import time
from gpbo.core.optutils import silentdirect as direct
from gpbo.core.optutils import geteffectiveoverhead
import logging
import copy
from gpbo.core import GPdc
from gpbo.core import PES
import ESutils
#start with random
import objectives

import costs
logger = logging.getLogger(__name__)

def randomaq(optstate,persist,**para):
    logger.info('randomaq')
    q = sp.random.uniform(size=len(para['lb']))
    return [l+x*(u-l) for l,u,x in zip(para['lb'],para['ub'],q)],para['ev'],persist,dict()


# and grid

def bruteaq(optstate,persist,**para):
    para = copy.deepcopy(para)
    if persist==None:
        persist = {'pwr':0,'idx':0,'d':len(para['ub'])}

    
    pwr = persist['pwr']
    idx = persist['idx']
    d = persist['d']
    k=2**pwr
    q=[0]*d
    logger.info('bruteaq griddiv={}'.format(k))
    for j in xrange(d):
        
        a,b = divmod(idx,k**(d-j-1))
        idx=b
        q[j]=(2*a+1)/float(2*k)
    
    
    if persist['idx']+1>= k**d:
        persist['pwr']+=1
        persist['idx']=0
    else:
        persist['idx']+=1
    return [l+x*(u-l) for l,u,x in zip(para['lb'],para['ub'],q)],para['ev'],persist,dict()


#EIMAP
def EIMAPaq(optstate,persist,**para):
    ev=para['ev']
    ub = para['ub']
    lb = para['lb']
    nrandinit = para['nrandinit']
    mprior = para['mprior']
    sprior = para['sprior']
    kindex = para['kindex']
    maxf = para['maxf']

    #para = copy.deepcopy(para)
    if persist==None:
        persist = {'n':0,'d':len(ub)}
    n = persist['n']
    d = persist['d']
    if n<nrandinit:
        persist['n']+=1
        return randomaq(optstate,persist,ev=ev,lb=lb,ub=ub)
    logger.info('EIMAPaq')

    if para['overhead']=='predict':
       overhead = geteffectiveoverhead(optstate,nrandinit)
    #logger.debug(sp.vstack([e[0] for e in optstate.ev]))
    #raise
    x=sp.vstack(optstate.x)
    y=sp.vstack(optstate.y)
    s= sp.vstack([e['s'] for e in optstate.ev])
    dx=[e['d'] for e in optstate.ev]
    MAP = GPdc.searchMAPhyp(x, y, s, dx, mprior, sprior, kindex)
    logger.info('found MAPHYP {}'.format(MAP))

    G = GPdc.GPcore(x, y, s, dx, GPdc.kernel(kindex, d, MAP))
    def directwrap(xq,y):
        xq.resize([1,d])
        a = G.infer_lEI(xq,[ev['d']])
        return (-a[0,0],0)



    [xmin,ymin,ierror] = direct(directwrap,lb,ub,user_data=[], algmethod=1, maxf=maxf, logfilename='/dev/null')

    logger.info('DIRECT found max EI at {} {}'.format(xmin,ierror))
    #logger.debug([xmin,ymin,ierror])
    persist['n']+=1
    return [i for i in xmin],ev,persist,{'MAPHYP':MAP,'logEImin':ymin,'DIRECTmessage':ierror}


#PES with fixed s ev
def PESfsaq(optstate,persist,**para):
    t0=time.clock()
    para = copy.deepcopy(para)
    if persist==None:
        persist = {'n':0,'d':len(para['ub']),'overhead':0.}
    n = persist['n']
    d = persist['d']
    if n<para['nrandinit']:
        persist['n']+=1
        
        return randomaq(optstate,persist,**para)
    logger.info('PESfsaq')
    #logger.debug(sp.vstack([e[0] for e in optstate.ev]))
    #raise
    x=sp.vstack(optstate.x)
    y=sp.vstack(optstate.y)
    s= sp.vstack([e['s'] for e in optstate.ev])
    dx=[e['d'] for e in optstate.ev]
    
    pesobj = PES.PES(x,y,s,dx,para['lb'],para['ub'],para['kindex'],para['mprior'],para['sprior'],DH_SAMPLES=para['DH_SAMPLES'],DM_SAMPLES=para['DM_SAMPLES'], DM_SUPPORT=para['DM_SUPPORT'],DM_SLICELCBPARA=para['DM_SLICELCBPARA'],mode=para['SUPPORT_MODE'],noS=para['noS'],DM_DROP=para['drop'])


    [xmin,ymin,ierror] = pesobj.search_pes(para['ev']['s'],maxf=para['maxf'])


    logger.info('DIRECT found max PES at {} {}'.format(xmin,ierror))
    lhyp = sp.log10([k.hyp for k in pesobj.G.kf])
    lhmean = sp.mean(lhyp, axis=0)
    lhstd = sp.sqrt(sp.var(lhyp, axis=0))
    lhmin = lhyp.min(axis=0)
    lhmax = lhyp.max(axis=0)
    logger.debug('loghyperparameters:\nmean {}\nstd {}\nmin {}\nmax {}'.format(lhmean,lhstd,lhmin,lhmax))

    persist['overhead']=time.clock()-t0
    return [i for i in xmin],para['ev'],persist,{'logHYPstats':{'mean':lhmean,'std':lhstd,'min':lhmin,'max':lhmax},'HYPdraws':[k.hyp for k in pesobj.G.kf],'mindraws':pesobj.Z,'DIRECTmessage':ierror,'PESmin':ymin,'kindex':para['kindex'],}



#PES with variable s ev give costfunction
def PESvsaq(optstate,persist,**para):
    t0=time.clock()
    para = copy.deepcopy(para)
    if persist==None:
        persist = {'n':0,'d':len(para['ub']),'overhead':0.}
    n = persist['n']
    d = persist['d']
    if n<para['nrandinit']:
        persist['n']+=1
        para2=copy.deepcopy(para)
        para2['ev']['s']=10**(para['logsu'])
        return randomaq(optstate,persist,**para2)
    logger.info('PESvsaq')
    #logger.debug(sp.vstack([e[0] for e in optstate.ev]))
    #raise
    x=sp.vstack(optstate.x)
    y=sp.vstack(optstate.y)
    s= sp.vstack([e['s'] for e in optstate.ev])
    dx=[e['d'] for e in optstate.ev]
    
    pesobj = PES.PES(x,y,s,dx,para['lb'],para['ub'],para['kindex'],para['mprior'],para['sprior'],DH_SAMPLES=para['DH_SAMPLES'],DM_SAMPLES=para['DM_SAMPLES'], DM_SUPPORT=para['DM_SUPPORT'],DM_SLICELCBPARA=para['DM_SLICELCBPARA'],mode=para['SUPPORT_MODE'],noS=para['noS'])

    if para['overhead']=='last':
        over=persist['overhead']
    else:
        over=0.
    [xmin,ymin,ierror] = pesobj.search_acq(para['cfn'],para['logsl'],para['logsu'],maxf=para['maxf'],over=over)
    
    logger.debug([xmin,ymin,ierror])
    para['ev']['s']=10**xmin[-1]
    xout = [i for i in xmin[:-1]]


    lhyp = sp.log10([k.hyp for k in pesobj.G.kf])
    lhmean = sp.mean(lhyp, axis=0)
    lhstd = sp.sqrt(sp.var(lhyp, axis=0))
    lhmin = lhyp.min(axis=0)
    lhmax = lhyp.max(axis=0)
    logger.debug('loghyperparameters:\nmean {}\nstd {}\nmin {}\nmax {}'.format(lhmean, lhstd, lhmin, lhmax))

    persist['overhead']=time.clock()-t0
    return xout,para['ev'],persist,{'logHYPstats':{'mean':lhmean,'std':lhstd,'min':lhmin,'max':lhmax},'HYPdraws':[k.hyp for k in pesobj.G.kf],'kindex':para['kindex'],'mindraws':pesobj.Z,'DIRECTmessage':ierror,'PESmin':ymin}



def PESbsaq(optstate,persist,**para):
    t0=time.clock()
    para = copy.deepcopy(para)
    if persist==None:
        persist = {'n':0,'d':len(para['ub']),'overhead':0.}
    n = persist['n']
    d = persist['d']
    if n<para['nrandinit'] and para['startmode']=='full':
        persist['n']+=1
        para['ev']['xa'] = sp.random.uniform(para['xal'],para['xau'])
        return randomaq(optstate,persist,**para)
    elif n<para['nrandinit'] and para['startmode']=='inline':

        r=persist['n']%len(para['initpoints'])
        if r==0:
            _x,_par,_per,_d=randomaq(optstate, persist, **para)
            persist['_x']=_x
            persist['_par'] = _par
            persist['_per'] = _per
            persist['_d'] = _d
        else:
            _x = persist['_x']
            _par = persist['_par']
            _per = persist['_per']
            _d = persist['_d']

        persist['n'] += 1


        _par['xa'] = para['initpoints'][r]
        return _x,_par,_per,_d
    elif n < para['nrandinit']:
        raise
    else:
        pass
    logger.info('PESbsaq')
    
    x=sp.hstack([sp.vstack([e['xa'] for e in optstate.ev]),sp.vstack(optstate.x)])
    
    y=sp.vstack(optstate.y)
    s= sp.vstack([e['s'] for e in optstate.ev])
    dx=[e['d'] for e in optstate.ev]
    #print "\n at pesinplane x {} axis 0".format(x)
    pesobj = PES.PES_inplane(x,y,s,dx,[para['xal']]+para['lb'],[para['xau']]+para['ub'],para['kindex'],para['mprior'],para['sprior'],0,0,DH_SAMPLES=para['DH_SAMPLES'], DM_SAMPLES=para['DM_SAMPLES'], DM_SUPPORT=para['DM_SUPPORT'],DM_SLICELCBPARA=para['DM_SLICELCBPARA'],mode=para['SUPPORT_MODE'],DH_CHAINS=para['hyp_chains'])
    if para['traincfn']:#
        #print "XXXXXXXXXXXXXXx"
        cx=sp.vstack([e['xa'] for e in optstate.ev])
        cc=sp.vstack([e for e in optstate.c])
        if para['traincfn']=='llog1d':
            cfn = costs.traincfn1dll(cx,cc)
        elif para['traincfn']=='llogfull':
            cfn = costs.traincfnfull(x,cc)
        elif para['traincfn']=='predictive1d':
            cfn = costs.predictive1d(cx,cc,sp.array(optstate.aqtime),para['nrandinit'],para['cmax']-optstate.C)
        else:
            #default is 1d nieve gp
            cfn = costs.traincfn1d(cx,cc)
    else:
        cfn = para['cfn']
    if para['overhead']=='last':
        over=persist['overhead']
    elif para['overhead']=='predict':
        over=geteffectiveoverhead(optstate,para['nrandinit'])
    else:
        over=0.
    [xmin,ymin,ierror] = pesobj.search_acq(cfn,lambda s:para['ev']['s'],maxf=para['maxf'],over=over)
    logger.debug([xmin,ymin,ierror])
    para['ev']['xa']=xmin[0]
    xout = [i for i in xmin[1:]]
    try:
        logger.debug('Predicted overhead {}'.format(cfn(xout,**{'xa':xmin[0]})))
    except e as e:
        print(e)
    lhyp = sp.log10([k.hyp for k in pesobj.G.kf])
    lhmean = sp.mean(lhyp, axis=0)
    lhstd = sp.sqrt(sp.var(lhyp, axis=0))
    lhmin = lhyp.min(axis=0)
    lhmax = lhyp.max(axis=0)
    logger.debug('loghyperparameters:\nmean {}\nstd {}\nmin {}\nmax {}'.format(lhmean, lhstd, lhmin, lhmax))
    persist['overhead']=time.clock()-t0
    return xout,para['ev'],persist,{'logHYPstats':{'mean':lhmean,'std':lhstd,'min':lhmin,'max':lhmax},'HYPdraws':[k.hyp for k in pesobj.G.kf],'kindex':para['kindex'],'mindraws':pesobj.Z,'DIRECTmessage':ierror,'PESmin':ymin}

def choiceaq(optstate,persist,**para):
    para = copy.deepcopy(para)
    if persist==None:
        persist = [None,[None]*len(para['aqoptions'])]
    aqn,choosepersist,transfer = para['chooser'](optstate,persist[0],**para['choosepara'])
    persist[0]=choosepersist
    para['aqoptions'][aqn][1]['transfer']=transfer
    logger.debug('choose to use aquisition {}'.format(para['aqoptions'][aqn][0].__name__))
    x,ev,pers,aux = para['aqoptions'][aqn][0](optstate,persist[1][aqn],**para['aqoptions'][aqn][1])
    persist[1][aqn]=pers
    return x,ev,persist,aux

def splocalaq(optstate,persist,**para):

    if persist==None:
        persist={'n':0,'y':[],'x':[],'done':False}
    else:
        persist['y'].append(optstate.y[-1])

    global count
    count=0

    if not 'start' in persist.keys():
        try:
            persist['start']=optstate.startlocal
        except:
            persist['start']=para['start']
    logger.info('splocalaq from {} step {}'.format(persist['start'],persist['n']))
    def fwrap(x):
        global count

        if count>=persist['n']:
            raise KeyError([i for i in x])
        else:
            assert sp.all(x==persist['x'][count])
            #print 'fwrap {} count {} y {}'.format(x,count,persist['y'][count])
            count+=1
            return persist['y'][count-1]
    try:
        minimize(fwrap,persist['start'],method='l-bfgs-b',bounds=[(para['lb'][i],para['ub'][i])for i in range(len(para['ub']))])
        persist['done']=True
        optstate.localdone=True
        return persist['x'][-1],para['ev'],persist,{'msg':'localopt is complete'}
    except KeyError as k:
        x=k.args[0]
    persist['x'].append(x)
    persist['n']+=1
    #print 'xtoev {}'.format(x)
    return x,para['ev'],persist,{'msg':'localopt' }