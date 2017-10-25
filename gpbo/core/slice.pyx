# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.


#slice sampling
#TODO need to check initial steps
from __future__ import print_function
xrange=range
import numpy as np
import scipy as sp
import emcee
from tqdm import tqdm
def slice_sample(loglike, init, iters, sigma, step_out=True,burn=20,subsam=4):
    """
    from http://isaacslavitt.com/2013/12/30/metropolis-hastings-and-
    slice-sampling/
    with some changes
    based on http://homepages.inf.ed.ac.uk/imurray2/teaching/09mlss/
    """

    # dist = joint_dist()
    
    # set up empty sample holder
    cdef int i,mn,D
    cdef double llh0,last_llh,llh_l,llh_r
    D = len(init)
    samples = sp.empty([int(iters+burn),D])
    
    # initialize
    xx = init.copy()
    pt = 0
    print( 'slice sample:')

    for i in tqdm(xrange(int(iters*subsam+burn*subsam))):
        mn=i-burn+1
        #print '\r Drawn %d    ' % mn,
        #sys.stdout.flush()

        perm = list(range(D))
        sp.random.shuffle(perm)
        
        last_llh = loglike(xx)

        for d in perm:
            llh0 = last_llh + sp.log(sp.random.rand())
            rr = sp.random.rand(1)
            x_l = xx.copy()
            x_l[d] = x_l[d] - rr * sigma[d]
            x_r = xx.copy()
            x_r[d] = x_r[d] + (1 - rr) * sigma[d]
            #z=0
            if step_out:
                
                llh_l = loglike(x_l)
                while llh_l > llh0:
                    # print x_l
                    x_l[d] = x_l[d] - sigma[d]
                    llh_l = loglike(x_l)
                    #z+=1
                llh_r = loglike(x_r)
                while llh_r > llh0:
                    x_r[d] = x_r[d] + sigma[d]
                    llh_r = loglike(x_r)
                    #z+=1
            #print "XXXXXXXXXXXXXXstepout counter "+str(z)
            x_cur = xx.copy()
            while True:
                xd = sp.random.rand() * (x_r[d] - x_l[d]) + x_l[d]
                x_cur[d] = xd.copy()
                last_llh = loglike(x_cur)
                if last_llh >= llh0:
                    xx[d] = xd.copy()
                    break
                elif xd > xx[d]:
                    x_r[d] = xd
                elif xd < xx[d]:
                    x_l[d] = xd
                else:
                    print( [xd,x_r[d],x_l[d],xx[d],last_llh,llh0])
                    print( [x_cur, xx, loglike(xx),llh0,last_llh])
                    
                    raise RuntimeError('Slice sampler shrank too far.')
        if i%subsam==0:
            
            samples[pt, :] = xx.copy().ravel()
            pt+=1
   
    return samples[burn:,:]

def sample_emcee(llk,p0,n,burn=100,subsam=17):
    nwalkers=10
    if isinstance(p0,float):
        D=1
    else:
        D = p0.size
    sampler = emcee.EnsembleSampler(nwalkers, D, llk ,args=[])
    pos,prob,state = sampler.run_mcmc(np.vstack([np.array(p0).flatten()]*nwalkers)+np.random.normal(size=[nwalkers,1])*1e-6, burn)
    sampler.reset()
    sampler.run_mcmc(pos,max(1,subsam*n//nwalkers))
    return sampler.flatchain[np.arange(n)*subsam]
