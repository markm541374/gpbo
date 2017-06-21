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
#from gpbo.core.optutils import silentdirect as direct
from gpbo.core.optutils import geteffectiveoverhead
import logging
import copy
import gpbo
from gpbo.core import GPdc
from gpbo.core import PES
from gpbo.core.optutils import multilocal

try:
    from matplotlib import pyplot as plt
    from matplotlib import patches
    plots=True
    plt.style.use('seaborn-paper')
except ImportError:
    plots=False
    plt=None
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

    if persist==None:
        persist = {'n':0,'d':len(ub)}
    n = optstate.n
    d = persist['d']
    if n<nrandinit:
        persist['n']+=1
        return randomaq(optstate,persist,ev=ev,lb=lb,ub=ub)
    logger.info('EIMAPaq')

    if para['overhead']=='predict':
       overhead = geteffectiveoverhead(optstate,nrandinit)
    x=sp.vstack(optstate.x)
    y=sp.vstack(optstate.y)
    s= sp.vstack([e['s']+10**optstate.condition for e in optstate.ev])
    dx=[e['d'] for e in optstate.ev]
    MAP = GPdc.searchMAPhyp(x, y, s, dx, mprior, sprior, kindex)
    logger.info('found MAPHYP {}'.format(MAP))

    G = GPdc.GPcore(x, y, s, dx, GPdc.kernel(kindex, d, mprior))

    G.m[0].set_parameter_dict(MAP)
    def wrap(x):
        xq = copy.copy(x)
        xq.resize([1,d])
        a = G.infer_lEI(xq,[ev['d']])
        return -a[0,0]
    print(wrap([0.,0.]))
    xmin,ymin,ierror = gpbo.core.optutils.twopartopt(wrap,para['lb'],para['ub'],para['dpara'],para['lpara'])
    #logger.debug([xmin,ymin,ierror])
    logger.info('localrefine found max EI at {} {} {}'.format(xmin,sp.exp(ymin),ierror))
    m,v = G.infer_diag_post(xmin,[[sp.NaN]])
    PIatX = sp.stats.norm.cdf(min(y),loc=m[0,0],scale=sp.sqrt(v[0,0]))
    persist['n']+=1
    return [i for i in xmin],ev,persist,{'MAPHYP':MAP,'logEImin':ymin,'DIRECTmessage':ierror,'EImax':sp.exp(-ymin),'PIatX':PIatX}


def eihypaq(optstate,persist,**para):
    t0=time.clock()
    para = copy.deepcopy(para)
    if persist==None:
        persist = {'n':0,'d':len(para['ub']),'overhead':0.,'raiseS':False}
    n = optstate.n
    d = persist['d']
    if n<para['nrandinit']:
        persist['n']+=1

        return randomaq(optstate,persist,**para)
    logger.info('EIHYPaq')
    x=sp.vstack(optstate.x)
    y=sp.vstack(optstate.y)
    s= sp.vstack([e['s']+10**optstate.condition for e in optstate.ev])
    dx=[e['d'] for e in optstate.ev]
    presetH=False
    if 'choosereturn' in para.keys():
        if 'reuseH' in para['choosereturn'].keys():
            presetH = para['choosereturn']['reuseH']
    if not presetH:
        G = PES.makeG(x,y,s,dx,para['kindex'],para['mprior'],para['sprior'],para['DH_SAMPLES'])
    else:
        logger.info('reusing preselected hyperparameters')
        G =  GPdc.GPcore(x,y,s,dx, [GPdc.kernel(para['kindex'], x.shape[1], h) for h in presetH])
    fixEI=False
    fixVal=0.
    if 'choosereturn' in para.keys():
        if 'offsetEI' in para['choosereturn'].keys():
            fixEI=True
            fixVal = para['choosereturn']['offsetEI']
            logger.info('EIoffset by {}'.format(fixVal))
    def wrap(Q):
        x = sp.array([Q])
        v = G.infer_lEI_post(x,[[sp.NaN]],fixI=fixEI,I=fixVal)[0,0]
        return -v

    xmin,ymin,ierror = gpbo.core.optutils.twopartopt(wrap,para['lb'],para['ub'],para['dpara'],para['lpara'])


    logger.info('DIRECT found max EI at {} {}'.format(xmin,ierror))
    lhyp = sp.log10([k.hyp for k in G.kf])
    lhmean = sp.mean(lhyp, axis=0)
    lhstd = sp.sqrt(sp.var(lhyp, axis=0))
    lhmin = lhyp.min(axis=0)
    lhmax = lhyp.max(axis=0)
    logger.debug('loghyperparameters:\nmean {}\nstd {}\nmin {}\nmax {}'.format(lhmean,lhstd,lhmin,lhmax))
    m,v = G.infer_diag_post(xmin,[[sp.NaN]])
    PIatX = sp.stats.norm.cdf(min(y),loc=m[0,0],scale=sp.sqrt(v[0,0]))
    persist['overhead']=time.clock()-t0
    return [i for i in xmin],para['ev'],persist,{'logHYPstats':{'mean':lhmean,'std':lhstd,'min':lhmin,'max':lhmax},'HYPdraws':[k.hyp for k in G.kf],'DIRECTmessage':ierror,'EImax':sp.exp(-ymin),'kindex':para['kindex'],'PIatX':PIatX}

#PES with fixed s ev
def PESfsaq(optstate,persist,**para):
    t0=time.clock()
    para = copy.deepcopy(para)
    if persist==None:
        persist = {'n':0,'d':len(para['ub']),'overhead':0.,'raiseS':False}
    n = optstate.n
    d = persist['d']
    if n<para['nrandinit']:
        persist['n']+=1
        
        return randomaq(optstate,persist,**para)
    logger.info('PESfsaq')
    #logger.debug(sp.vstack([e[0] for e in optstate.ev]))
    #raise
    x=sp.vstack(optstate.x)
    y=sp.vstack(optstate.y)
    s= sp.vstack([e['s']+10**optstate.condition for e in optstate.ev])
    dx=[e['d'] for e in optstate.ev]
    presetH=False
    if 'choosereturn' in para.keys():
        if 'reuseH' in para['choosereturn'].keys():
            presetH = para['choosereturn']['reuseH']
    pesobj = PES.PES(x,y,s,dx,para['lb'],para['ub'],para['kindex'],para['mprior'],para['sprior'],DH_SAMPLES=para['DH_SAMPLES'],DM_SAMPLES=para['DM_SAMPLES'], DM_SUPPORT=para['DM_SUPPORT'],DM_SLICELCBPARA=para['DM_SLICELCBPARA'],mode=para['SUPPORT_MODE'],noS=para['noS'],DM_DROP=para['drop'],preselectH=presetH)
    [xmin,ymin,ierror] = pesobj.search_pes(para['ev']['s'],para)

    logger.info('DIRECT found max PES at {} {}'.format(xmin,ierror))
    lhyp = sp.log10([k.hyp for k in pesobj.G.kf])
    lhmean = sp.mean(lhyp, axis=0)
    lhstd = sp.sqrt(sp.var(lhyp, axis=0))
    lhmin = lhyp.min(axis=0)
    lhmax = lhyp.max(axis=0)
    logger.debug('loghyperparameters:\nmean {}\nstd {}\nmin {}\nmax {}'.format(lhmean,lhstd,lhmin,lhmax))

    m,v = pesobj.G.infer_diag_post(xmin,[[sp.NaN]])
    PIatX = sp.stats.norm.cdf(min(y),loc=m[0,0],scale=sp.sqrt(v[0,0]))
    persist['overhead']=time.clock()-t0
    return [i for i in xmin],para['ev'],persist,{'logHYPstats':{'mean':lhmean,'std':lhstd,'min':lhmin,'max':lhmax},'HYPdraws':[k.hyp for k in pesobj.G.kf],'mindraws':pesobj.Z,'DIRECTmessage':ierror,'PESmin':ymin,'kindex':para['kindex'],'PIatX':PIatX}



def vmaxaq(optstate,persist,**para):
    t0=time.clock()
    para = copy.deepcopy(para)
    if persist==None:
        persist = {'n':0,'d':len(para['ub']),'overhead':0.,'raiseS':False}
    n = optstate.n
    d = persist['d']
    if n<para['nrandinit']:
        persist['n']+=1

        return randomaq(optstate,persist,**para)
    logger.info('vmaxaq')
    #logger.debug(sp.vstack([e[0] for e in optstate.ev]))
    #raise
    x=sp.vstack(optstate.x)
    y=sp.vstack(optstate.y)
    s= sp.vstack([e['s']+10**optstate.condition for e in optstate.ev])
    dx=[e['d'] for e in optstate.ev]
    presetH=False
    if 'choosereturn' in para.keys():
        if 'reuseH' in para['choosereturn'].keys():
            presetH = para['choosereturn']['reuseH']
    if not presetH:
        G = PES.makeG(x,y,s,dx,para['kindex'],para['mprior'],para['sprior'],para['DH_SAMPLES'])
    else:
        logger.info('reusing preselected hyperparameters')
        G =  GPdc.GPcore(x,y,s,dx, [GPdc.kernel(para['kindex'], x.shape[1], h) for h in presetH])
    def wrap(Q):
        x = sp.array([Q])
        v = G.infer_diag_post(x,[[sp.NaN]])[1][0,0]
        return -v

    xmin,ymin,ierror = gpbo.core.optutils.twopartopt(wrap,para['lb'],para['ub'],para['dpara'],para['lpara'])
    vmax = -ymin


    logger.info('DIRECT found max PES at {} {}'.format(xmin,ierror))
    lhyp = sp.log10([k.hyp for k in G.kf])
    lhmean = sp.mean(lhyp, axis=0)
    lhstd = sp.sqrt(sp.var(lhyp, axis=0))
    lhmin = lhyp.min(axis=0)
    lhmax = lhyp.max(axis=0)
    logger.debug('loghyperparameters:\nmean {}\nstd {}\nmin {}\nmax {}'.format(lhmean,lhstd,lhmin,lhmax))

    persist['overhead']=time.clock()-t0
    if gpbo.core.debugoutput['datavis']:
        fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(45, 45))
        # plot the current GP
        if d==2:
            gpbo.core.optutils.gpplot(ax[0,0],ax[0,1],G,para['lb'],para['ub'],ns=60)
            ax[0,0].set_title('GP_post_mean')
            ax[0,1].set_title('GP_post_var')
            ax[1, 1].plot(x[:,0], x[:,1], 'ro')

        try:
            fname = 'vmaxaqplot' + time.strftime('%d_%m_%y_%H:%M:%S') + '.png'
            print('saving as {}'.format(fname))
            fig.savefig(os.path.join(gpbo.core.debugoutput['path'], fname))
        except BaseException as e:
            logger.error(str(e))
        fig.clf()
        plt.close(fig)
    return [i for i in xmin],para['ev'],persist,{'logHYPstats':{'mean':lhmean,'std':lhstd,'min':lhmin,'max':lhmax},'HYPdraws':[k.hyp for k in G.kf],'DIRECTmessage':ierror,'PESmin':ymin,'kindex':para['kindex'],}

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
        if para['sinitrand']:
            srel = sp.random.uniform()
            switch = sp.random.uniform()
            if srel>switch:
                srel=1-srel
            para2['ev']['s']=10**(para['logsu']+srel*(para['logsl']-para['logsu']))

        else:
            para2['ev']['s']=10**(para['logsu'])
        return randomaq(optstate,persist,**para2)
    logger.info('PESvsaq')
    #logger.debug(sp.vstack([e[0] for e in optstate.ev]))
    #raise
    x=sp.vstack(optstate.x)
    y=sp.vstack(optstate.y)
    s= sp.vstack([e['s']+10**optstate.condition for e in optstate.ev])
    dx=[e['d'] for e in optstate.ev]
    
    pesobj = PES.PES(x,y,s,dx,para['lb'],para['ub'],para['kindex'],para['mprior'],para['sprior'],DH_SAMPLES=para['DH_SAMPLES'],DM_SAMPLES=para['DM_SAMPLES'], DM_SUPPORT=para['DM_SUPPORT'],DM_SLICELCBPARA=para['DM_SLICELCBPARA'],mode=para['SUPPORT_MODE'],noS=para['noS'])
    if para['traincfn']:#
        #print "XXXXXXXXXXXXXXx"
        cx=sp.vstack([sp.log10(e['s']) for e in optstate.ev])
        cc=sp.vstack([e for e in optstate.c])
        if para['traincfn']=='llog1d':
            cfnlog = costs.traincfn1dll(cx,cc)
        elif para['traincfn']=='llogfull':
            cfnlog = costs.traincfnfull(x,cc)
        elif para['traincfn']=='predictive1d':
            cfnlog = costs.predictive1d(cx,cc,sp.array(optstate.aqtime),para['nrandinit'],para['cmax']-optstate.C)
        else:
            #default is 1d nieve gp
            cfnlog = costs.traincfn1d(cx,cc)
        def cfn(x,**ev):
            #print(x)
            #print(ev)
            return cfnlog(x,**{'xa':sp.log10(ev['s'])})
    else:
        cfn = para['cfn']

    if para['overhead']=='last':
        over=persist['overhead']
    elif para['overhead']=='predict':
        over=geteffectiveoverhead(optstate,para['nrandinit'])
    else:
        over=0.
    [xmin,ymin,ierror] = pesobj.search_acq(cfn,para['logsl'],para['logsu'],para,over=over)
    
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
    return xout,para['ev'],persist,{'overheadprediction':over,'logHYPstats':{'mean':lhmean,'std':lhstd,'min':lhmin,'max':lhmax},'HYPdraws':[k.hyp for k in pesobj.G.kf],'kindex':para['kindex'],'mindraws':pesobj.Z,'DIRECTmessage':ierror,'PESmin':ymin}



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
    s= sp.vstack([e['s']+10**optstate.condition for e in optstate.ev])
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
    [xmin,ymin,ierror] = pesobj.search_acq(cfn,lambda s:para['ev']['s'],para,over=over)
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
    #logger.error( str(persist))
    if persist==None:
        persist={'n':0,'y':[],'z':[],'done':False}
        for k in para['choosereturn'].keys():
            persist[k]=para['choosereturn'][k]

        if 'H' in persist.keys():
            R = sp.linalg.cholesky(persist['H']).T
            persist['R']=R
        else:
            logger.debug('no precondition provided')
            persist['R']=sp.eye(len(persist['start']))
    else:
        persist['y'].append(optstate.y[-1])

    global count
    count=0
    logger.info('splocalaq from {} ({}) step {}'.format(persist['start'],persist['R'].dot(persist['start']),persist['n']))
    def fwrap(z):
        global count

        if count>=persist['n']:
            raise KeyError([i for i in z])
        else:
            assert sp.all(z==persist['z'][count])
            #print 'fwrap {} count {} y {}'.format(x,count,persist['y'][count])
            count+=1
            return persist['y'][count-1]
    try:
        R=minimize(fwrap,persist['R'].dot(persist['start']),method='bfgs',options={'gtol':0.000001})
        persist['done']=True
        optstate.localdone=True
        logger.info('localopt finished with z: {} (x: {}) y: {} {}'.format(R.x,sp.linalg.solve(persist['R'],persist['z'][-1]),R.fun,R.message))
        return list(sp.linalg.solve(persist['R'],persist['z'][-1])),para['ev'],persist,{'msg':'localopt is complete {}'.format(str(R))}
    except KeyError as k:
        z=k.args[0]
    persist['z'].append(z)
    persist['n']+=1
    #print 'xtoev {}'.format(x)
    x = sp.linalg.solve(persist['R'],z)
    return list(x),para['ev'],persist,{'msg':'localopt' }