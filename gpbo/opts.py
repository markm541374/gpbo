from __future__ import print_function
xrange=range

import gpbo
import copy
import os
import matplotlib
from matplotlib import pyplot as plt
#plt.style.use('seaborn-paper')
#plt.rc('font',serif='Times')




import matplotlib.patches as mpatches
import scipy as sp
print('removed robo import in opts.py l17 due to theano errors')
#try:
#    from gpbo.exps.thirdwrap.mtbowrap import optmtbo
from gpbo.exps.thirdwrap.fabwrap import optfabolas
from gpbo.exps.thirdwrap.fabwrap import optfabolas_mod
#except:
#    print('\n\ndidnt import robo!!!!!!\n\n')
def runexp(f,lb,ub,path,nreps,confs,indexoffset=0):
    for i_ in range(nreps):
        ii=i_+indexoffset
        for C in confs:
            if C[0][:2]=='ei':
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
                try:
                    out = gpbo.search(C[1])
                except:
                    pass

            elif C[0][:5]=='pesbs':
                C[1].path = path
                C[1].fname = '{}_{}.csv'.format(C[0],ii)
                C[1].aqpara['lb'] = [i for i in lb]
                C[1].aqpara['ub'] = [i for i in ub]
                C[1].reccpara['lb'] = [i for i in lb]
                C[1].reccpara['ub'] = [i for i in ub]
                C[1].ojf=f

                try:
                    out = gpbo.search(C[1])
                except:
                    pass
            elif C[0][:9]=='switching':
                C[1].path = path
                C[1].fname = '{}_{}.csv'.format(C[0],ii)
                C[1].ojf=f
                try:
                    out = gpbo.search(C[1])
                except:
                    pass
            elif C[0][:5]=='pesvs':
                C[1].path = path
                C[1].fname = '{}_{}.csv'.format(C[0],ii)
                C[1].aqpara['lb'] = [i for i in lb]
                C[1].aqpara['ub'] = [i for i in ub]
                C[1].reccpara['lb'] = [i for i in lb]
                C[1].reccpara['ub'] = [i for i in ub]
                C[1].ojf=f
                try:
                    out = gpbo.search(C[1])
                except:
                    pass
            elif C[0][:4]=='mtbo':
                try:
                    optmtbo(f, lb, ub, 1.-(1./C[1]['lowtask']), C[1]['nsteps'], ninit=C[1]['ninit'],fpath=path,fname='{}_{}.csv'.format(C[0],ii),mod=C[1]['switchestimator'])
                except:
                    pass
            elif C[0][:7]=='fabolas':
                try:
                    optfabolas(f,lb,ub,C[1]['nsteps'],C[1]['ninit'],fname='{}_{}.csv'.format(C[0],ii), fpath=path)
                except:
                    raise
            elif C[0][:6]=='fabmod':
                try:
                    optfabolas_mod(f,lb,ub,C[1]['nsteps'],C[1]['ninit'],fname='{}_{}.csv'.format(C[0],ii), fpath=path,switchestimator=C[1]['switchestimator'],switchkernel=C[1]['switchkernel'],timelimit=C[1]['timelimit'])
                except:
                    raise
            else:
                print( "not an optimization method")
def plotquarts(a,data1,data2,col,line,lab,log=False):
    n=len(data1)
    mx=-sp.Inf
    mn=sp.Inf
    for i in xrange(n):
        mn = min(mn,min(data1[i]))
        mx = max(mx,max(data1[i]))
    if not log:
        xaxis = sp.linspace(mn,mx,200)
    else:
        xaxis = sp.logspace(sp.log10(mn),sp.log10(mx),200)+1e-9
##    print( data1)
#    print( data2)
#    print( xaxis)
    low0, med0, upp0 = gpbo.core.ESutils.quartsirregular(data1,data2,xaxis)

    #        a.fill_between(xaxis, low0, upp0, facecolor='lightblue', edgecolor='lightblue', alpha=0.5)
    a.plot(xaxis, med0, color=col, linestyle=line, label=lab)
    a.fill_between(xaxis,upp0,low0,edgecolor=col, linestyle=line,facecolor=col,lw=0.0,alpha=0.1)
    return


def getquarts(data1,data2,log=False):
    n=len(data1)
    mx=-sp.Inf
    mn=sp.Inf
    for i in xrange(n):
        mn = min(mn,min(data1[i]))
        mx = max(mx,max(data1[i]))
    if not log:
        xaxis = sp.linspace(mn,mx,200)
    else:
        xaxis = sp.logspace(sp.log10(mn),sp.log10(mx),200)

    low0, med0, upp0 = gpbo.core.ESutils.quartsirregular(data1,data2,xaxis)
    return xaxis,low0,med0,upp0

def getmvint(data1,data2,nstd=2,logx=False,logy=False):
    n=len(data1)
    mx=-sp.Inf
    mn=sp.Inf
    for i in xrange(n):
        mn = min(mn,min(data1[i]))
        mx = max(mx,max(data1[i]))
    if not logx:
        xaxis = sp.linspace(mn,mx,200)
    else:
        xaxis = sp.logspace(sp.log10(mn),sp.log10(mx),200)
    if logy:
        M,V = gpbo.core.ESutils.mvirregular(data1,[d.apply(sp.log) for d in data2],xaxis)
        S = nstd*sp.sqrt(V)
        return xaxis,sp.exp(M-S),sp.exp(M),sp.exp(M+S)
    else:
        M,V = gpbo.core.ESutils.mvirregular(data1,data2,xaxis)
        S = nstd*sp.sqrt(V)
        return xaxis,M-S,M,M+S

def plotquartsends(a,xdata_, ydata_,col,line,lab,log=False,mean=False):
    xdata = [sp.array(i) for i in xdata_]
    ydata = [sp.array(i) for i in ydata_]
    n = len(xdata)
    ints = []
    starts = sp.empty(n)
    ends = sp.empty(n)
    yends = sp.empty(n)
    for i in xrange(n):
        starts[i] = xdata[i][0]
        ends[i] = xdata[i][-1]
        yends[i] = ydata[i][-1]
    yendorder = sp.argsort(yends)

    Ydata = [sp.hstack([y[0], y, y[-1]]) for y in ydata]
    #the pad values are slightly outside the true range to so that exp(log(value)) stays in the interpolation range
    Xdata = [sp.hstack([0.999*min(starts), x, max(ends)*1.001]) for x in xdata]
    #print(min(starts),max(ends))
    for i in xrange(n):
        #print(Xdata[i][0],Xdata[i][-1])
        ints.append(sp.interpolate.interp1d(Xdata[i], Ydata[i]))
       # a.plot(Xdata[i], Ydata[i], 'lightblue')

    if log:
        x = sp.logspace(sp.log10(min(starts)), sp.log10(max(ends)), 200)
    else:
        x = sp.linspace(min(starts), max(ends), 200)

    #print(x)
    if mean:
        a.plot(x, map(lambda x: sp.mean([i(x) for i in ints]), x), color=col,label=lab)
    else:
        a.plot(x, map(lambda x: sp.percentile([i(x) for i in ints], 50), x), color=col,label=lab)
        #m = map(lambda x: sp.mean([i(x) for i in ints]), x)
        #v = map(lambda x: sp.mean([i(x) for i in ints]), x)
        #a.plot(x, map(lambda x: sp.mean([i(x) for i in ints]), x), color=col,label=lab)

        y25 = map(lambda x: sp.percentile([i(x) for i in ints], 25), x)
        y75 = map(lambda x: sp.percentile([i(x) for i in ints], 75), x)
        a.fill_between(x,y25,y75,edgecolor=col, facecolor=col,lw=0.0,alpha=0.1)
    #a.plot(ends[yendorder], yends[yendorder], '.',color=col ,linestyle=line)
    #print("endvalues: {}".format(yends))
    a2 = a.twinx()
    a2.grid(False)
    a2.plot(ends[sp.argsort(ends)],sp.linspace(1,0,n),color=col, linestyle='--',linewidth=0.4)
    a2.set_ylabel('fraction of optimizations still running')
    return

def plotall(confs,nreps,path,trueopt=False,logx=False,labelfn = lambda x:x,axisset=dict(),skipinit=False,sixylabel=False,thirteenylabel=False,allylabel=False,showends=False,needed=None,legend=True,forcefigsize=None,xmax=None,legendreplace=None):
    if showends:
        pq=plotquartsends
    else:
        pq=plotquarts
    if needed is None:
        needed=range(20)
    f=[]
    a=[]
    pmax=21
    class fakeax:
        def plot(self,*args,**kwargs):
            return
        def set_yscale(self,*args,**kwargs):
            return
        def set_xscale(self,*args,**kwargs):
            return
    for i in range(pmax):
        if i in needed:
            if forcefigsize is None:
                f_,a_ = plt.subplots(1)
            else:
                f_,a_ = plt.subplots(1,figsize=forcefigsize)

            f.append(f_)
            a.append(a_)
        else:
            f.append(None)
            a.append(fakeax())
    colorlist = plt.rcParams['axes.prop_cycle'].by_key()['color']
    #['b','r','g','purple','k','grey','orange','c','lightgreen','lightblue','pink','b','r','g','purple','k','grey','orange','c','lightgreen','lightblue','pink']
    ci=-1

    allmomin=sp.Inf
    for C in confs:
        if  C[0][:5]=='pesbs' or C[0][:3]=='fab':
            ninit=0
        else:
            ninit=0
        if not skipinit:
            ninit=0
        ci+=1
        col = colorlist[ci]
        line = 'solid'#lslist[ci]
        #collect the data
        data=[]
        for ii in range(nreps):
            thisdata = gpbo.optimize.readoptdata(os.path.join(path,'{}_{}.csv'.format(C[0],ii)))
            allmomin = min(allmomin,thisdata['trueyatxrecc'].values.min())
            data.append(thisdata)

        #if True:
            #first plot is all the opts per step
            #termregret = sp.mean([list(d['trueyatxrecc'])[-1] for d in data])
        if 0 in needed:
            for ii in range(nreps):
                a[0].plot(data[ii]['index'],data[ii]['trueyatxrecc'],color=col, linestyle=line,label=labelfn(C[0]))
            #and averaged
        if 4 in needed:
            pq(a[4],[data[k]['index'] for k in range(nreps)],[data[k]['trueyatxrecc'] for k in range(nreps)],col,line,labelfn(C[0]))

        if 1 in needed:
            #second is all the opts per evaluation cost
            for ii in range(nreps):
                a[1].plot(data[ii]['accE'],data[ii]['trueyatxrecc'],color=col, linestyle=line,label=labelfn(C[0]))
            #and averaged
        if 5 in needed:
            pq(a[5],[data[k]['accE'] for k in range(nreps)],[data[k]['trueyatxrecc'] for k in range(nreps)],col,line,labelfn(C[0]),log=logx)

        if 2 in needed:
            #third is all the opts per evaluation + acquisition cost
            for ii in range(nreps):
                a[2].plot(data[ii]['accEA'],data[ii]['trueyatxrecc'],color=col, linestyle=line,label=labelfn(C[0]))
        if 6 in needed:
            #and averaged
            pq(a[6],[data[k]['accEA'] for k in range(nreps)],[data[k]['trueyatxrecc'] for k in range(nreps)],col,line,labelfn(C[0]),log=logx)

        if 3 in needed:
            #fourth is evcost per step
            for ii in range(nreps):
                a[3].plot(data[ii]['index'],data[ii]['c'],color=col, linestyle=line,label=labelfn(C[0]))
        if 7 in needed:
            #and averaged
            pq(a[7],[data[k]['index'] for k in range(nreps)],[data[k]['c'] for k in range(nreps)],col,line,labelfn(C[0]))

        if 14 in needed:
             #fiifth is overhead clock time
            for ii in range(nreps):
                a[14].plot(data[ii]['index'],data[ii]['taq'],color=col, linestyle=line,label=labelfn(C[0]))
        if 15 in needed:
            #and averaged
            pq(a[15],[data[k]['index'] for k in range(nreps)],[data[k]['taq'] for k in range(nreps)],col,line,labelfn(C[0]))

        if (16 in needed) or (17 in needed):
            #fiifth is overhead clock time
            try:
                for ii in range(nreps):
                    a[16].plot(data[ii]['index'],data[ii]['xa'],color=col, linestyle=line,label=labelfn(C[0]))
                #and averaged
                pq(a[17],[data[k]['index'] for k in range(nreps)],[data[k]['xa'] for k in range(nreps)],col,line,labelfn(C[0]))
            except:
                for ii in range(nreps):
                    a[16].plot(data[ii]['index'],data[ii]['s'],color=col, linestyle=line,label=labelfn(C[0]))
                #and averaged
                pq(a[17],[data[k]['index'] for k in range(nreps)],[data[k]['s'] for k in range(nreps)],col,line,labelfn(C[0]))
                try:
                    a[16].set_yscale('log')
                    a[17].set_yscale('log')
                except:
                    pass
        if (18 in needed) or (19 in needed):
            #try:
            for ii in range(nreps):
                a[18].plot(data[ii]['index'],data[ii]['condition'],color=col, linestyle=line,label=labelfn(C[0]))
            #and averaged
            pq(a[19],[data[k]['index'] for k in range(nreps)],[data[k]['condition'] for k in range(nreps)],col,line,labelfn(C[0]))
            a[18].set_yscale('log')
            a[19].set_yscale('log')
            #except:
            #    pass
        if trueopt:
            if 8 in needed:
                #first plot is all the opts per step
                for ii in range(nreps):
                    a[8].plot(data[ii]['index'],data[ii]['trueyatxrecc']-trueopt,color=col, linestyle=line,label=labelfn(C[0]))
            if 11 in needed:
                #and averaged
                pq(a[11],[data[k]['n'] for k in range(nreps)],[data[k]['trueyatxrecc']-trueopt for k in range(nreps)],col,line,labelfn(C[0]))

            if 9 in needed:
                #second is all the opts per evaluation cost
                for ii in range(nreps):
                    a[9].plot(data[ii]['accE'],data[ii]['trueyatxrecc']-trueopt,color=col, linestyle=line,label=labelfn(C[0]))
            if 12 in needed:
                #and averaged
                pq(a[12],[data[k]['accE'] for k in range(nreps)],[data[k]['trueyatxrecc']-trueopt for k in range(nreps)],col,line,labelfn(C[0]),log=True)

            if 10 in needed:
                #third is all the opts per evaluation + acquisition cost
                for ii in range(nreps):
                    a[10].plot(data[ii]['accEA'],data[ii]['trueyatxrecc']-trueopt,color=col, linestyle=line,label=labelfn(C[0]))
            if 13 in needed:
                #and averaged
                pq(a[13],[data[k]['accEA'][ninit:] for k in range(nreps)],[data[k]['trueyatxrecc'][ninit:]-trueopt for k in range(nreps)],col,line,labelfn(C[0]),log=True)

            if 20 in needed:
                #and averaged
                pq(a[20],[data[k]['n'] for k in range(nreps)],[data[k]['trueyatxrecc']-trueopt for k in range(nreps)],col,line,labelfn(C[0]),mean=True)
    print('allmomin {}'.format(allmomin))
    if 0 in needed:
        a[0].legend()
        a[0].set_xlabel('Steps')
        a[0].set_ylabel('result')
        if 0 in axisset.keys():
            a[0].axis(axisset[0])
        f[0].savefig(os.path.join(path,'out0.pdf'),bbox_inches='tight', pad_inches=0.1)

    if 1 in needed:
        a[1].legend()
        a[1].set_xlabel('Evaluation Cost')
        a[1].set_ylabel('result')
        if logx:
            a[1].set_xscale('log')
        a[1].axis('tight')
        f[1].savefig(os.path.join(path,'out1.pdf'),bbox_inches='tight', pad_inches=0.1)

    if 2 in needed:
        #a[2].legend()
        a[2].set_xlabel('Evaluation+Acquisition Cost')
        a[2].set_ylabel('result')
        if logx:
            a[2].set_xscale('log')
        a[2].axis('tight')
        f[2].savefig(os.path.join(path,'out2.pdf'),bbox_inches='tight', pad_inches=0.1)

    if 3 in needed:
        a[3].legend()
        a[3].set_xlabel('Steps')
        a[3].set_ylabel('EVcost')
        f[3].savefig(os.path.join(path,'out3.pdf'),bbox_inches='tight', pad_inches=0.1)

    if 4 in needed:
        a[4].legend()
        a[4].set_xlabel('Steps')
        a[4].set_ylabel('result')
        f[4].savefig(os.path.join(path,'out4.pdf'),bbox_inches='tight', pad_inches=0.1)

    if 5 in needed:
        if legend:
            a[5].legend()
        a[5].set_xlabel('Evaluation Cost (s)')
        if sixylabel:
            a[5].set_ylabel(sixylabel)
        elif allylabel:
            a[5].set_ylabel(allylabel)
        else:
            a[5].set_ylabel('Result')
        if logx:
            a[5].set_xscale('log')
        a[5].axis('tight')
        if 5 in axisset.keys():
            a[5].axis(axisset[5])
        f[5].savefig(os.path.join(path,'out5.pdf'),bbox_inches='tight', pad_inches=0.1)

    if 6 in needed:
        if legend:
            a[6].legend()
        a[6].set_xlabel('Total Clock Time (s)')
        if sixylabel:
            a[6].set_ylabel(sixylabel)
        elif allylabel:
            a[6].set_ylabel(allylabel)
        else:
            a[6].set_ylabel('Result')
        if logx:
            a[6].set_xscale('log')
        a[6].axis('tight')
        if 6 in axisset.keys():
            a[6].axis(axisset[6])
        f[6].savefig(os.path.join(path,'out6.pdf'),bbox_inches='tight', pad_inches=0.1)

    if 7 in needed:
        a[7].legend()
        a[7].set_xlabel('Steps')
        a[7].set_ylabel('EVcost')
        f[7].savefig(os.path.join(path,'out7.pdf'),bbox_inches='tight', pad_inches=0.1)

    if 14 in needed:
        a[14].legend()
        a[14].set_xlabel('Steps')
        a[14].set_ylabel('Overhead Clocktime (s)')
        f[14].savefig(os.path.join(path,'out14.pdf'),bbox_inches='tight', pad_inches=0.1)

    if 15 in needed:
        #:a[15].legend()
        a[15].set_xlabel('Steps')
        a[15].set_ylabel('Overhead Clocktime (s)')
        if 15 in axisset.keys():
            a[15].axis(axisset[15],'tight')
        f[15].savefig(os.path.join(path,'out15.pdf'),bbox_inches='tight', pad_inches=0.1)

    if 16 in needed:
        a[16].legend()
        a[16].set_xlabel('Steps')
        a[16].set_ylabel('env Var')
        f[16].savefig(os.path.join(path,'out16.pdf'),bbox_inches='tight', pad_inches=0.1)

    if 17 in needed:
        a[17].legend()
        a[17].set_xlabel('Steps')
        a[17].set_ylabel('env Var')
        f[17].savefig(os.path.join(path,'out17.pdf'),bbox_inches='tight', pad_inches=0.1)

    if 18 in needed:
        a[18].legend()
        a[18].set_xlabel('Steps')
        a[18].set_ylabel('condition magnitude')
        f[18].savefig(os.path.join(path,'out18.pdf'),bbox_inches='tight', pad_inches=0.1)

    if 19 in needed:
        if legend:
            a[19].legend()
        a[19].set_xlabel('Steps')
        a[19].set_ylabel('Condition Number')
        f[19].savefig(os.path.join(path,'out19.pdf'),bbox_inches='tight', pad_inches=0.1)

    if trueopt:
        if 8 in needed:
            a[8].set_xlabel('Steps')
            a[8].set_ylabel('regret')
            a[8].set_yscale('log')
            f[8].savefig(os.path.join(path,'out8.pdf'),bbox_inches='tight', pad_inches=0.1)

        if 9 in needed:
            a[9].set_xlabel('Evaluation Cost (s)')
            a[9].set_ylabel('regret')
            a[9].set_yscale('log')
            f[9].savefig(os.path.join(path,'out9.pdf'),bbox_inches='tight', pad_inches=0.1)

        if 10 in needed:
            a[10].set_xlabel('Evaluation+Acquisition Cost')
            a[10].set_ylabel('Median Immediate Regret')
            a[10].set_yscale('log')
            if logx:
                a[10].set_xscale('log')
            a[10].axis('tight')
            f[10].savefig(os.path.join(path,'out10.pdf'),bbox_inches='tight', pad_inches=0.1)

        if 11 in needed:
            if legend:
                if legend=='outside':
                    a[11].legend(bbox_to_anchor=(1.14,1), loc="upper left")
                else:
                    a[11].legend()

            a[11].set_xlabel('Steps')
            a[11].set_ylabel('Median Immediate Regret')
            a[11].set_yscale('log')
            if logx:
                a[11].set_xscale('log')
            a[11].axis('tight')
            if 11 in axisset.keys():
                a[11].axis(axisset[11],'tight')
            if xmax:
                a[11].set_xlim(0,xmax)
            f[11].savefig(os.path.join(path,'out11.pdf'),bbox_inches='tight', pad_inches=0.1)

        if 12 in needed:
            if legend:
                a[12].legend()
            a[12].set_xlabel('Evaluation Cost')
            if allylabel:
                a[12].set_ylabel(allylabel)
            else:
                a[12].set_ylabel('Median Immediate Regret')
            a[12].set_yscale('log')
            if logx:
                a[12].set_xscale('log')
            a[12].axis('tight')
            if 12 in axisset.keys():
                a[12].axis(axisset[12])
            f[12].savefig(os.path.join(path,'out12.pdf'),bbox_inches='tight', pad_inches=0.1)

        if 13 in needed:
            if legend:
                a[13].legend()
            a[13].set_xlabel('Total Clock Time (s)')
            if thirteenylabel:
                a[13].set_ylabel(thirteenylabel)
            elif allylabel:
                a[13].set_ylabel(allylabel)
            else:
                a[13].set_ylabel('Median Immediate Regret')
            a[13].set_yscale('log')
            if logx:
                a[13].set_xscale('log')
            a[13].axis('tight')
            if 13 in axisset.keys():
                a[13].axis(axisset[13],'tight')
            f[13].savefig(os.path.join(path,'out13.pdf'),bbox_inches='tight', pad_inches=0.1)

        if 20 in needed:
            if legend=='override':
                a[20].legend(loc="lower left",handles=legendreplace)
            elif not legend:
                pass
            else:
                a[20].legend(loc="lower left")

            a[20].set_xlabel('Steps')
            a[20].set_ylabel('Mean Immediate Regret')
            a[20].set_yscale('log')
            if logx:
                a[20].set_xscale('log')
            a[20].axis('tight')
            if 20 in axisset.keys():
                a[20].axis(axisset[20],'tight')
            a[20].set_xlim(0,140)
            f[20].savefig(os.path.join(path,'out20.pdf'),bbox_inches='tight', pad_inches=0.1)

def plotprofile(confs,nreps,path,tol=0.9,target=1e-6):
    f=[]
    a=[]
    pmax=1
    for i in range(pmax):
        f_,a_ = plt.subplots(1)
        for item in ([a_.title, a_.xaxis.label, a_.yaxis.label] + a_.get_xticklabels() + a_.get_yticklabels()):
            item.set_fontsize(10)
        f.append(f_)
        a.append(a_)
    colorlist = ['b','r','g','purple','k','grey','orange','c','lightgreen','lightblue','pink','b','r','g','purple','k','grey','orange','c','lightgreen','lightblue','pink']
    lslist = ['solid' , 'dashed', 'dashdot', 'dotted','solid' , 'dashed', 'dashdot', 'dotted','solid' , 'dashed', 'dashdot', 'dotted','solid' , 'dashed', 'dashdot', 'dotted','solid' , 'dashed', 'dashdot', 'dotted']
    ci=-1

    for C in confs:
        print('plotting {}...'.format(C[0]))
        ci+=1
        col = colorlist[ci]
        line = lslist[ci]
        #collect the data
        data=[]
        xoverheads = sp.empty(nreps)
        xntotarget = sp.zeros(nreps)
        overheads = sp.empty(nreps)
        noverheads = [sp.empty(nreps) for i in range(len(C[1]['N']))]
        ninrun = sp.zeros(nreps)
        support = sp.logspace(-2,5,200)
        success = sp.zeros(nreps)
        for ii in range(nreps):
            D = gpbo.optimize.readoptdata(os.path.join(path,'{}_{}.csv'.format(C[0],ii)))
            A = sp.array(D['trueyatxrecc'].values)
            if A.min()>=target:
                xoverheads[ii]=sum(D['taq'])
                overheads[ii]=sum(D['taq'])
                for k,n in enumerate(C[1]['N']):
                    noverheads[k][ii] = sum(D['taq'][:n])
            else:
                success[ii]=1
                i = sp.argmax(A<=target)#while A[i]>=target:
                xoverheads[ii] = sum(D['taq'][:i+1])
                overheads[ii] = sum(D['taq'])
                xntotarget[ii] = i
                ninrun[ii]=len(D['taq'])
                for k,n in enumerate(C[1]['N']):
                    noverheads[k][ii] = sum(D['taq'].values[:n])
        if sp.mean(success)>=tol:
            if C[1]['oracle']:
                a[0].plot(support,sp.mean(xoverheads)+sp.mean(xntotarget)*support,col,label=C[0]+'oracle',linestyle='dashdot')
            if C[1]['full']:
                a[0].plot(support,sp.mean(overheads)+sp.mean(ninrun)*support,col,label=C[0]+'_all',linestyle='dashed')
            for k,n in enumerate(C[1]['N']):
                if sp.percentile(xntotarget,int(tol*100))<n:
                    a[0].plot(support,sp.mean(noverheads[k])+n*support,col,label=C[0]+str(n),linestyle='solid')
        else:
            print('{} only acchived target on {}'.format(C[0],sp.mean(success)))
        a[0].set_xscale('log')
        a[0].set_yscale('log')
        a[0].legend()

    f[0].savefig(os.path.join(path,'profile_{}.png'.format(sp.log10(target))),bbox_inches='tight', pad_inches=0.1)

