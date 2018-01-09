import numpy as np
import scipy as sp
import gpbo
import pandas as pd
from matplotlib import pyplot as plt
import pickle
import os
if False:
    #regenerate the data (implicitly saved as results/bounddata.p)
    gpbo.core.debugoutput['adaptive']=True
    #df = pd.read_csv('/home/mark/gpbo/exps/stopping/colville/results/s_3.csv')
    df = pd.read_csv('results/switchingfigdata.csv')
    dim=2
    X = np.vstack([df[' x{}'.format(i)].values for i in range(dim)]).T
    Y = np.array([df[' y'].values]).T
    print(Y.min())
    S = np.array([[0.]*Y.size]).T
    D=[[np.NaN]]*Y.size
    O = gpbo.core.optimize.optstate()
    n=70
    for i in range(n):
        O.update(X[i,:],{'d':[np.NaN],'s':0.},Y[i,0],1.,1.)
    #C = gpbo.core.config.eihypdefault(lambda x:0,4,10,0.,'','')
    C=gpbo.core.config.switchdefault(lambda x:0,dim,10,250,0.,'','')
    C.choosepara['rotate']=False
    persist = {'n':n,'d':dim,'flip':False,'raiseS':False,'R':np.eye(2)}

    for k in ['localsampleregret','Rexists','sampleregret','expectedregret','localrsam','localrest','expectedRbinorm']:
        persist[k]=[]
    #C.aqpara['choosereturn']={'offsetEI':0}
    #x,par,per,aux = gpbo.core.acquisitions.eihypaq(O,persist,**C.aqpara)
    #x1,par,per,aux = gpbo.core.acquisitions.eihypaq(O,persist,**C.aqpara)
    #persist['R']=R
    #xr,par,per,aux = gpbo.core.acquisitions.eihypaq(O,persist,**C.aqpara)
    x = gpbo.core.choosers.globallocalregret(O,persist,**C.choosepara)

def pltcdf(Y,C,ax,col):
    return ax.plot(sp.hstack([[i,i] for i in Y])[1:-1],sp.hstack([[i-C[0],i] for i in C])[1:-1],color=col,label='Sampled CDF')
def plotfig():
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    sup,mu,vvmax,mvmax,ymin,Yout,Cout = pickle.load(open('results/bounddata.p','r'))
    f2,a2 = plt.subplots(1)
    pltcdf(Yout,Cout,a2,colors[1])
    a2.set_yscale('logit')

    a2.plot(sup,sp.stats.norm.cdf(sup,loc=mu,scale=sp.sqrt(vvmax)),color=colors[0],linestyle='--', label='Approx Tail Upper Bound')
    a2.plot(sup,sp.stats.norm.cdf(sup,loc=mvmax,scale=sp.sqrt(vvmax)),color=colors[0],linestyle='-.',label='Lower Bound')
    a2.axvline(ymin,label='Posterior Mean Minimum',color='k',linestyle=':')
    a2.set_ylabel('CDF')
    a2.set_xlabel('y')
    from matplotlib.ticker import NullFormatter
    a2.yaxis.set_minor_formatter(NullFormatter())
    a2.spines['left']._adjust_location()

    a2.legend()
    plt.tight_layout()
    f2.savefig('figs/ends.pdf')
plotfig()