# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.
import scipy as sp
import os
import time
import DIRECT
import logging
import copy
import GPdc
import PES
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
def EIMAPaq(optstate,persist,ev=None, ub = None, lb=None, nrandinit=None, mprior=None,sprior=None,kindex = None,volper=None):
    #para = copy.deepcopy(para)
    if persist==None:
        persist = {'n':0,'d':len(ub)}
    n = persist['n']
    d = persist['d']
    if n<nrandinit:
        persist['n']+=1
        return randomaq(optstate,persist,ev=ev,lb=lb,ub=ub)
    logger.info('EIMAPaq')
    #logger.debug(sp.vstack([e[0] for e in optstate.ev]))
    #raise
    x=sp.vstack(optstate.x)
    y=sp.vstack(optstate.y)
    s= sp.vstack([e['s'] for e in optstate.ev])
    dx=[e['d'] for e in optstate.ev]
    MAP = GPdc.searchMAPhyp(x, y, s, dx, mprior, sprior, kindex)
    logger.info('MAPHYP {}'.format(MAP))

    G = GPdc.GPcore(x, y, s, dx, GPdc.kernel(kindex, d, MAP))
    def directwrap(xq,y):
        xq.resize([1,d])
        a = G.infer_lEI(xq,[ev['d']])
        return (-a[0,0],0)
    
    [xmin,ymin,ierror] = DIRECT.solve(directwrap,lb,ub,user_data=[], algmethod=1, volper = volper, logfilename='/dev/null')
    #logger.debug([xmin,ymin,ierror])
    persist['n']+=1
    return [i for i in xmin],ev,persist,{'MAPHYP':MAP,'logEImin':ymin,'DIRECTmessage':ierror}


#PES with fixed s ev
def PESfsaq(optstate,persist,**para):
    para = copy.deepcopy(para)
    if persist==None:
        persist = {'n':0,'d':len(para['ub'])}
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
    
    pesobj = PES.PES(x,y,s,dx,para['lb'],para['ub'],para['kindex'],para['mprior'],para['sprior'],DH_SAMPLES=para['DH_SAMPLES'],DM_SAMPLES=para['DM_SAMPLES'], DM_SUPPORT=para['DM_SUPPORT'],DM_SLICELCBPARA=para['DM_SLICELCBPARA'],mode=para['SUPPORT_MODE'],noS=para['noS'])


    [xmin,ymin,ierror] = pesobj.search_pes(para['ev']['s'],volper=para['volper'])

    lhyp = sp.log10([k.hyp for k in pesobj.G.kf])
    lhmean = sp.mean(lhyp, axis=0)
    lhstd = sp.sqrt(sp.var(lhyp, axis=0))
    lhmin = lhyp.min(axis=0)
    lhmax = lhyp.max(axis=0)
    logger.debug('loghyperparameters:\nmean {}\nstd {}\nmin {}\nmax {}'.format(lhmean,lhstd,lhmin,lhmax))

    return [i for i in xmin],para['ev'],persist,{'logHYPstats':{'mean':lhmean,'std':lhstd,'min':lhmin,'max':lhmax},'HYPdraws':[k.hyp for k in pesobj.G.kf],'mindraws':pesobj.Z,'DIRECTmessage':ierror,'PESmin':ymin,'kindex':para['kindex'],}



#PES with variable s ev give costfunction
def PESvsaq(optstate,persist,**para):
    para = copy.deepcopy(para)
    if persist==None:
        persist = {'n':0,'d':len(para['ub'])}
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
    [xmin,ymin,ierror] = pesobj.search_acq(para['cfn'],para['logsl'],para['logsu'],volper=para['volper'])
    
    logger.debug([xmin,ymin,ierror])
    para['ev']['s']=10**xmin[-1]
    xout = [i for i in xmin[:-1]]

    lhyp = sp.log10([k.hyp for k in pesobj.G.kf])
    lhmean = sp.mean(lhyp, axis=0)
    lhstd = sp.sqrt(sp.var(lhyp, axis=0))
    lhmin = lhyp.min(axis=0)
    lhmax = lhyp.max(axis=0)
    logger.debug('loghyperparameters:\nmean {}\nstd {}\nmin {}\nmax {}'.format(lhmean, lhstd, lhmin, lhmax))

    return xout,para['ev'],persist,{'logHYPstats':{'mean':lhmean,'std':lhstd,'min':lhmin,'max':lhmax},'HYPdraws':[k.hyp for k in pesobj.G.kf],'kindex':para['kindex'],'mindraws':pesobj.Z,'DIRECTmessage':ierror,'PESmin':ymin}



def PESbsaq(optstate,persist,**para):
    para = copy.deepcopy(para)
    if persist==None:
        persist = {'n':0,'d':len(para['ub'])}
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
    pesobj = PES.PES_inplane(x,y,s,dx,[para['xal']]+para['lb'],[para['xau']]+para['ub'],para['kindex'],para['mprior'],para['sprior'],0,0,DH_SAMPLES=para['DH_SAMPLES'], DM_SAMPLES=para['DM_SAMPLES'], DM_SUPPORT=para['DM_SUPPORT'],DM_SLICELCBPARA=para['DM_SLICELCBPARA'],mode=para['SUPPORT_MODE'])
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
        
    [xmin,ymin,ierror] = pesobj.search_acq(cfn,lambda s:para['ev']['s'],volper=para['volper'])
    logger.debug([xmin,ymin,ierror])
    para['ev']['xa']=xmin[0]
    xout = [i for i in xmin[1:]]

    lhyp = sp.log10([k.hyp for k in pesobj.G.kf])
    lhmean = sp.mean(lhyp, axis=0)
    lhstd = sp.sqrt(sp.var(lhyp, axis=0))
    lhmin = lhyp.min(axis=0)
    lhmax = lhyp.max(axis=0)
    logger.debug('loghyperparameters:\nmean {}\nstd {}\nmin {}\nmax {}'.format(lhmean, lhstd, lhmin, lhmax))

    return xout,para['ev'],persist,{'logHYPstats':{'mean':lhmean,'std':lhstd,'min':lhmin,'max':lhmax},'HYPdraws':[k.hyp for k in pesobj.G.kf],'kindex':para['kindex'],'mindraws':pesobj.Z,'DIRECTmessage':ierror,'PESmin':ymin}
