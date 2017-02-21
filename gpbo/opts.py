from __future__ import print_function
xrange=range

import gpbo
import copy
import os
from matplotlib import pyplot as plt
plt.rc('font',serif='Times')
import matplotlib.patches as mpatches
import scipy as sp
try:
    from gpbo.exps.thirdwrap.mtbowrap import optmtbo
    from gpbo.exps.thirdwrap.fabwrap import optfabolas
    from gpbo.exps.thirdwrap.fabwrap import optfabolas_mod
except:
    print('\n\ndidnt import robo!!!!!!\n\n')
def runexp(f,lb,ub,path,nreps,confs,indexoffset=0):
    for i_ in range(nreps):
        ii=i_+indexoffset
        for C in confs:
            if C[0][:5]=='eimle':
                C[1].path=path
                C[1].fname='{}_{}.csv'.format(C[0],ii)
                C[1].aqpara['lb']=lb
                C[1].aqpara['ub']=ub
                C[1].reccpara['lb']=lb
                C[1].reccpara['ub']=ub
                def wrap(x,**ev):
                    tmp = copy.deepcopy(ev)
                    tmp['xa']=0
                    return f(x,**tmp)
                C[1].ojf = wrap
                try:
                    out = gpbo.search(C[1])
                except:
                    pass
            elif C[0][:5]=='pesfs':
                C[1].path = path
                C[1].fname = '{}_{}.csv'.format(C[0],ii)
                C[1].aqpara['lb'] = lb
                C[1].aqpara['ub'] = ub
                C[1].reccpara['lb'] = lb
                C[1].reccpara['ub'] = ub
                def wrap(x, **ev):
                    tmp = copy.deepcopy(ev)
                    tmp['xa'] = 0
                    return f(x, **tmp)
                C[1].ojf = wrap
                #try:
                out = gpbo.search(C[1])
                #except:
                #    pass

            elif C[0][:5]=='pesbs':
                C[1].path = path
                C[1].fname = '{}_{}.csv'.format(C[0],ii)
                C[1].aqpara['lb'] = [i for i in lb]
                C[1].aqpara['ub'] = [i for i in ub]
                C[1].reccpara['lb'] = [i for i in lb]
                C[1].reccpara['ub'] = [i for i in ub]
                C[1].ojf=f

                #try:
                out = gpbo.search(C[1])
                #except:
                #    pass
            elif C[0][:4]=='mtbo':
                try:
                    optmtbo(f, lb, ub, 1.-(1./C[1]['lowtask']), C[1]['nsteps'], ninit=C[1]['ninit'],fpath=path,fname='{}_{}.csv'.format(C[0],ii),mod=C[1]['switchestimator'])
                except:
                    pass
            elif C[0][:7]=='fabolas':
                try:
                    optfabolas(f,lb,ub,C[1]['nsteps'],C[1]['ninit'],fname='{}_{}.csv'.format(C[0],ii), fpath=path)
                except:
                    pass
            elif C[0][:6]=='fabmod':
                try:
                    optfabolas_mod(f,lb,ub,C[1]['nsteps'],C[1]['ninit'],fname='{}_{}.csv'.format(C[0],ii), fpath=path,switchestimator=C[1]['switchestimator'],switchkernel=C[1]['switchkernel'])
                except:
                    pass
            else:
                print( "not an optimization method")
def plotquarts(a,data1,data2,col,lab):
    n=len(data1)
    mx=-sp.Inf
    mn=sp.Inf
    for i in xrange(n):
        mn = min(mn,min(data1[i]))
        mx = max(mx,max(data1[i]))
    xaxis = sp.linspace(mn,mx,200)
##    print( data1)
#    print( data2)
#    print( xaxis)
    low0, med0, upp0 = gpbo.core.ESutils.quartsirregular(data1,data2,xaxis)

    #        a.fill_between(xaxis, low0, upp0, facecolor='lightblue', edgecolor='lightblue', alpha=0.5)
    a.plot(xaxis, med0, color=col, label=lab)
    a.fill_between(xaxis,upp0,low0,edgecolor=None,facecolor=col,alpha=0.1)
    return



def plotall(confs,nreps,path,trueopt=False,logx=False,labelfn = lambda x:x):
    f=[]
    a=[]
    pmax=20
    for i in range(pmax):
        f_,a_ = plt.subplots(1)
        f.append(f_)
        a.append(a_)
    colorlist = ['b','r','g','purple','k','grey','orange','c','lightgreen','lightblue','pink']
    ci=-1
    for C in confs:
        ci+=1
        col = colorlist[ci]
        #collect the data
        data=[]
        for ii in range(nreps):
            data.append(gpbo.optimize.readoptdata(os.path.join(path,'{}_{}.csv'.format(C[0],ii))))

        if True:
            #first plot is all the opts per step

            for ii in range(nreps):
                a[0].plot(data[ii]['index'],data[ii]['trueyatxrecc'],color=col,label=labelfn(C[0]))
            #and averaged
            plotquarts(a[4],[data[k]['index'] for k in range(nreps)],[data[k]['trueyatxrecc'] for k in range(nreps)],col,labelfn(C[0]))

            #second is all the opts per evaluation cost
            for ii in range(nreps):
                a[1].plot(data[ii]['accE'],data[ii]['trueyatxrecc'],color=col,label=labelfn(C[0]))
            #and averaged
            plotquarts(a[5],[data[k]['accE'] for k in range(nreps)],[data[k]['trueyatxrecc'] for k in range(nreps)],col,labelfn(C[0]))

            #third is all the opts per evaluation + acquisition cost
            for ii in range(nreps):
                a[2].plot(data[ii]['accEA'],data[ii]['trueyatxrecc'],color=col,label=labelfn(C[0]))
            #and averaged
            plotquarts(a[6],[data[k]['accEA'] for k in range(nreps)],[data[k]['trueyatxrecc'] for k in range(nreps)],col,labelfn(C[0]))

            #fourth is evcost per step
            for ii in range(nreps):
                a[3].plot(data[ii]['index'],data[ii]['c'],color=col,label=labelfn(C[0]))
            #and averaged
            plotquarts(a[7],[data[k]['index'] for k in range(nreps)],[data[k]['c'] for k in range(nreps)],col,labelfn(C[0]))

             #fiifth is overhead clock time
            for ii in range(nreps):
                a[14].plot(data[ii]['index'],data[ii]['taq'],color=col,label=labelfn(C[0]))
            #and averaged
            plotquarts(a[15],[data[k]['index'] for k in range(nreps)],[data[k]['taq'] for k in range(nreps)],col,labelfn(C[0]))

            #fiifth is overhead clock time
            try:
                for ii in range(nreps):
                    a[16].plot(data[ii]['index'],data[ii]['xa'],color=col,label=labelfn(C[0]))
                #and averaged
                plotquarts(a[17],[data[k]['index'] for k in range(nreps)],[data[k]['xa'] for k in range(nreps)],col,labelfn(C[0]))
            except:
                pass
        if trueopt:
            #first plot is all the opts per step
            for ii in range(nreps):
                a[8].plot(data[ii]['index'],data[ii]['trueyatxrecc']-trueopt,color=col,label=labelfn(C[0]))
            #and averaged
            plotquarts(a[11],[data[k]['index'] for k in range(nreps)],[data[k]['trueyatxrecc']-trueopt for k in range(nreps)],col,labelfn(C[0]))

            #second is all the opts per evaluation cost
            for ii in range(nreps):
                a[9].plot(data[ii]['accE'],data[ii]['trueyatxrecc']-trueopt,color=col,label=labelfn(C[0]))
            #and averaged
            plotquarts(a[12],[data[k]['accE'] for k in range(nreps)],[data[k]['trueyatxrecc']-trueopt for k in range(nreps)],col,labelfn(C[0]))

            #third is all the opts per evaluation + acquisition cost
            for ii in range(nreps):
                a[10].plot(data[ii]['accEA'],data[ii]['trueyatxrecc']-trueopt,color=col,label=labelfn(C[0]))
            #and averaged
            plotquarts(a[13],[data[k]['accEA'] for k in range(nreps)],[data[k]['trueyatxrecc']-trueopt for k in range(nreps)],col,labelfn(C[0]))


    a[0].legend()
    a[0].set_xlabel('steps')
    a[0].set_ylabel('result')
    f[0].savefig(os.path.join(path,'out0.png'),bbox_inches='tight', pad_inches=0.1)

    a[1].legend()
    a[1].set_xlabel('Evaluation Cost')
    a[1].set_ylabel('result')
    if logx:
        a[1].set_xscale('log')
    a[1].axis('tight')
    f[1].savefig(os.path.join(path,'out1.png'),bbox_inches='tight', pad_inches=0.1)

    a[2].legend()
    a[2].set_xlabel('Evaluation+Acquisition Cost')
    a[2].set_ylabel('result')
    if logx:
        a[2].set_xscale('log')
    a[2].axis('tight')
    f[2].savefig(os.path.join(path,'out2.png'),bbox_inches='tight', pad_inches=0.1)

    a[3].legend()
    a[3].set_xlabel('steps')
    a[3].set_ylabel('EVcost')
    f[3].savefig(os.path.join(path,'out3.png'),bbox_inches='tight', pad_inches=0.1)

    a[4].legend()
    a[4].set_xlabel('steps')
    a[4].set_ylabel('result')
    f[4].savefig(os.path.join(path,'out4.png'),bbox_inches='tight', pad_inches=0.1)

    a[5].legend()
    a[5].set_xlabel('Evaluation Cost')
    a[5].set_ylabel('result')
    if logx:
        a[5].set_xscale('log')
    a[5].axis('tight')
    f[5].savefig(os.path.join(path,'out5.png'),bbox_inches='tight', pad_inches=0.1)

    a[6].legend()
    a[6].set_xlabel('Evaluation+Acquisition Cost')
    a[6].set_ylabel('result')
    if logx:
        a[6].set_xscale('log')
    a[6].axis('tight')
    f[6].savefig(os.path.join(path,'out6.png'),bbox_inches='tight', pad_inches=0.1)

    a[7].legend()
    a[7].set_xlabel('steps')
    a[7].set_ylabel('EVcost')
    f[7].savefig(os.path.join(path,'out7.png'),bbox_inches='tight', pad_inches=0.1)

    a[14].legend()
    a[14].set_xlabel('steps')
    a[14].set_ylabel('overhead clocktime')
    f[14].savefig(os.path.join(path,'out14.png'),bbox_inches='tight', pad_inches=0.1)

    a[15].legend()
    a[15].set_xlabel('steps')
    a[15].set_ylabel('overhead clocktime')
    f[15].savefig(os.path.join(path,'out15.png'),bbox_inches='tight', pad_inches=0.1)

    a[16].legend()
    a[16].set_xlabel('steps')
    a[16].set_ylabel('env Var')
    f[16].savefig(os.path.join(path,'out16.png'),bbox_inches='tight', pad_inches=0.1)

    a[17].legend()
    a[17].set_xlabel('steps')
    a[17].set_ylabel('env Var')
    f[17].savefig(os.path.join(path,'out17.png'),bbox_inches='tight', pad_inches=0.1)

    if trueopt:
        a[8].set_xlabel('steps')
        a[8].set_ylabel('regret')
        a[8].set_yscale('log')
        f[8].savefig(os.path.join(path,'out8.png'),bbox_inches='tight', pad_inches=0.1)

        a[9].set_xlabel('Evaluation Cost')
        a[9].set_ylabel('regret')
        a[9].set_yscale('log')
        f[9].savefig(os.path.join(path,'out9.png'),bbox_inches='tight', pad_inches=0.1)

        a[10].set_xlabel('Evaluation+Acquisition Cost')
        a[10].set_ylabel('regret')
        a[10].set_yscale('log')
        if logx:
            a[10].set_xscale('log')
        a[10].axis('tight')
        f[10].savefig(os.path.join(path,'out10.png'),bbox_inches='tight', pad_inches=0.1)

        a[11].legend()
        a[11].set_xlabel('steps')
        a[11].set_ylabel('regret')
        a[11].set_yscale('log')
        if logx:
            a[11].set_xscale('log')
        a[11].axis('tight')
        f[11].savefig(os.path.join(path,'out11.png'),bbox_inches='tight', pad_inches=0.1)

        a[12].legend()
        a[12].set_xlabel('Evaluation Cost')
        a[12].set_ylabel('regret')
        a[12].set_yscale('log')
        if logx:
            a[12].set_xscale('log')
        a[12].axis('tight')
        f[12].savefig(os.path.join(path,'out12.png'),bbox_inches='tight', pad_inches=0.1)

        a[13].legend()
        a[13].set_xlabel('Total Clock Time')
        a[13].set_ylabel('Median IR')
        a[13].set_yscale('log')
        if logx:
            a[13].set_xscale('log')
        a[13].axis('tight')
        f[13].savefig(os.path.join(path,'out13.png'),bbox_inches='tight', pad_inches=0.1)

