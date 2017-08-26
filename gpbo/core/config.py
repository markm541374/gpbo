import gpbo
xrange=range
from gpbo.core import GPdc as GPdc
import scipy as sp


class eimledefault():
    """
    fixed s, space is [-1,1]^D

    """
    def __init__(self,f,D,n,s,path,fname):
        self.aqfn = gpbo.core.acquisitions.EIMAPaq
        self.aqpara= {
            'ev': {'s': s, 'd': [sp.NaN]},
            'lb': [-1.]*D,
            'ub': [1.]*D,
            'nrandinit': 10,
            'mprior': sp.array([1.]+[-0.3]*D),
            'sprior': sp.array([2.]*(D+1)),
            'kindex': GPdc.MAT52,
            'maxf':500+100*D,
            'overhead':None,
            'dpara': {'user_data': [],
                      'algmethod': 1,
                      'maxf': 500+100*D,
                      'logfilename': '/dev/null'},
            'lpara': {'gtol': 0.00001,
                      'maxfun': 200}
        }

        self.stoppara = {'nmax': n}
        self.stopfn = gpbo.core.optimize.nstopfn

        self.reccfn = gpbo.core.reccomenders.gpmaprecc
        self.reccpara = {
                'ev':self.aqpara['ev'],
                'lb':self.aqpara['lb'],
                'ub':self.aqpara['ub'],
                'mprior':self.aqpara['mprior'],
                'sprior':self.aqpara['sprior'],
                'kindex':self.aqpara['kindex'],
                'maxf':500+100*D,
                'onlyafter':self.aqpara['nrandinit'],
                'check':True,
                'smode':'direct',
                'dpara':self.aqpara['dpara'],
                'lpara':self.aqpara['lpara'],
                'everyn':1
                }
        self.ojfchar = {'dx': len(self.aqpara['lb']), 'dev': len(self.aqpara['ev'])}
        self.ojf=f

        self.path = path
        self.fname = fname
        return

class eimlelearns():
    """
    fixed s, space is [-1,1]^D

    """
    def __init__(self,f,D,n,s,path,fname):
        self.aqfn = gpbo.core.acquisitions.EIMAPaq
        self.aqpara= {
            'ev': {'s': s, 'd': [sp.NaN]},
            'lb': [-1.]*D,
            'ub': [1.]*D,
            'nrandinit': 10,
            'mprior': sp.array([1.]+[0.]*D+[-2]),
            'sprior': sp.array([1.]*(D+1)+[5]),
            'kindex': GPdc.SQUEXPCS,
            'maxf':500+100*D,
            'dpara': {'user_data': [],
                      'algmethod': 1,
                      'maxf': 500+100*D,
                      'logfilename': '/dev/null'},
            'lpara': {'gtol': 0.00001,
                      'maxfun': 200}
        }

        self.stoppara = {'nmax': n}
        self.stopfn = gpbo.core.optimize.nstopfn

        self.reccfn = gpbo.core.reccomenders.gpmaprecc
        self.reccpara = {
                'ev':self.aqpara['ev'],
                'lb':self.aqpara['lb'],
                'ub':self.aqpara['ub'],
                'mprior':self.aqpara['mprior'],
                'sprior':self.aqpara['sprior'],
                'kindex':self.aqpara['kindex'],
                'maxf':500+100*D,
                'onlyafter':self.aqpara['nrandinit'],
                'check':True,
                'dpara':self.aqpara['dpara'],
                'lpara':self.aqpara['lpara'],
                'everyn':1
                }
        self.ojfchar = {'dx': len(self.aqpara['lb']), 'dev': len(self.aqpara['ev'])}
        self.ojf=f

        self.path = path
        self.fname = fname
        return

class eihypdefault():
    def __init__(self,f,D,n,s,path,fname):
        self.aqfn = gpbo.core.acquisitions.eihypaq
        self.aqpara= {
            'ev': {'s': s, 'd': [sp.NaN]},
            'lb': [-1.]*D,
            'ub': [1.]*D,
            'nrandinit': 10,
            #'maxf':500+100*D,
            'mprior': sp.array([1.]+[0.]*D),
            'sprior': sp.array([1.]*(D+1)),
            'kindex': GPdc.MAT52,
            'DH_SAMPLES': 16+6*D,
            'drop':True,
            'noS': False,
            'dpara': {'user_data': [],
                      'algmethod': 1,
                      'maxf': 500+100*D,
                      'logfilename': '/dev/null'},
            'lpara': {'gtol': 0.00001,
                      'maxfun': 200}
        }

        self.stoppara = {'nmax': n}
        self.stopfn = gpbo.core.optimize.nstopfn

        self.reccfn = gpbo.core.reccomenders.gphinrecc
        self.reccpara = {
            'ev':self.aqpara['ev'],
            'lb':self.aqpara['lb'],
            'ub':self.aqpara['ub'],
            'mprior':self.aqpara['mprior'],
            'sprior':self.aqpara['sprior'],
            'kindex':self.aqpara['kindex'],
            'maxf':500+100*D,
            'onlyafter':self.aqpara['nrandinit'],
            'check':True,
            'dpara':self.aqpara['dpara'],
            'lpara':self.aqpara['lpara'],
            'everyn':1
        }
        self.ojfchar = {'dx': len(self.aqpara['lb']), 'dev': len(self.aqpara['ev'])}
        self.ojf=f

        self.path = path
        self.fname = fname
        return
class pesfsdefault():
    def __init__(self,f,D,n,s,path,fname,ninit=10):
        self.aqfn = gpbo.core.acquisitions.PESfsaq
        self.aqpara= {
            'ev': {'s': s, 'd': [sp.NaN]},
            'lb': [-1.]*D,
            'ub': [1.]*D,
            'nrandinit': ninit,
            #'maxf':500+100*D,
            'mprior': sp.array([1.]+[0.]*D),
            'sprior': sp.array([1.]*(D+1)),
            'kindex': GPdc.MAT52,
            'DH_SAMPLES': 16+6*D,
            'DM_SAMPLES': 20+8*D,
            'DM_SUPPORT': 750+250*D,
            'SUPPORT_MODE': [gpbo.core.ESutils.SUPPORT_LAPAPROT],
            'DM_SLICELCBPARA': 12+4.*D,
            'drop':True,
            'overhead':'none',
            'noS': False,
            'dpara': {'user_data': [],
                      'algmethod': 1,
                      'maxf': 500+100*D,
                      'logfilename': '/dev/null'},
            'lpara': {'gtol': 0.00001,
                      'maxfun': 200}
        }

        self.stoppara = {'nmax': n}
        self.stopfn = gpbo.core.optimize.nstopfn

        self.reccfn = gpbo.core.reccomenders.gphinrecc
        self.reccpara = {
                'ev':self.aqpara['ev'],
                'lb':self.aqpara['lb'],
                'ub':self.aqpara['ub'],
                'mprior':self.aqpara['mprior'],
                'sprior':self.aqpara['sprior'],
                'kindex':self.aqpara['kindex'],
                'maxf':500+100*D,
                'onlyafter':self.aqpara['nrandinit'],
                'check':True,
                'dpara':self.aqpara['dpara'],
                'lpara':self.aqpara['lpara'],
                'everyn':1
                }
        self.ojfchar = {'dx': len(self.aqpara['lb']), 'dev': len(self.aqpara['ev'])}
        self.ojf=f

        self.path = path
        self.fname = fname
        return

class pesvsdefault():
    def __init__(self,f,cfn,D,n,lsl,lsu,path,fname):
        self.aqfn = gpbo.core.acquisitions.PESvsaq
        self.aqpara= {
            'ev': {'s': lsu, 'd': [sp.NaN]},
            'lb': [-1.]*D,
            'ub': [1.]*D,
            'nrandinit': 10,
            'maxf': 500+100*D,
            'mprior': sp.array([1.]+[0.]*D),
            'sprior': sp.array([1.]*(D+1)),
            'kindex': GPdc.MAT52,
            'DH_SAMPLES': 16+6*D,
            'DM_SAMPLES': 20+8*D,
            'DM_SUPPORT': 750+250*D,
            'SUPPORT_MODE': [gpbo.core.ESutils.SUPPORT_LAPAPROT],
            'DM_SLICELCBPARA': 12+4*D,
            'noS': False,
            'logsu': lsu,
            'logsl': lsl,
            'sinitrand':True,
            'overhead':'None',
            'cfn': cfn,
            'traincfn':'llog1d',
            'dpara': {'user_data': [],
                      'algmethod': 1,
                      'maxf': 500+100*D,
                      'logfilename': '/dev/null'},
            'lpara': {'gtol': 0.00001,
                      'maxfun': 200}
        }

        self.stoppara = {'nmax': n}
        self.stopfn = gpbo.core.optimize.nstopfn

        self.reccfn = gpbo.core.reccomenders.gphinrecc
        self.reccpara = {
                'ev':self.aqpara['ev'],
                'lb':self.aqpara['lb'],
                'ub':self.aqpara['ub'],
                'mprior':self.aqpara['mprior'],
                'sprior':self.aqpara['sprior'],
                'kindex':self.aqpara['kindex'],
                'maxf':500+100*D,
                'onlyafter':self.aqpara['nrandinit'],
                'check':True,
                'dpara':self.aqpara['dpara'],
                'lpara':self.aqpara['lpara'],
                'everyn':1
                }
        self.ojfchar = {'dx': len(self.aqpara['lb']), 'dev': len(self.aqpara['ev'])}
        self.ojf=f

        self.path = path
        self.fname = fname
        return

class pesbsdefault():
    def __init__(self,f,D,n,s,path,fname):
        self.aqfn = gpbo.core.acquisitions.PESbsaq
        self.aqpara = {
            'ev':{'s':s,'d':[sp.NaN],'xa':0.},
            'lb':[-1.] * D,
            'ub':[ 1.] * D,
            'maxf':500+250*(D+1),
            'mprior': sp.array([1.]+[0.]*(D+1)),
            'sprior': sp.array([1.]*(D+2)),
            'kindex':GPdc.MAT52,
            'DH_SAMPLES':16+6*D,
            'DM_SAMPLES':20+8*D,
            'hyp_chains':1,
            'DM_SUPPORT':750+250*D,
            'SUPPORT_MODE':[gpbo.core.ESutils.SUPPORT_LAPAPROT],
            'DM_SLICELCBPARA':12+4*D,
            'noS':False,
            'nrandinit':20,
            'cfn':lambda x,ev:42.,
            'traincfn':True,
            'xau':1.,
            'xal':0.,
            'startmode':'inline',
            'initpoints':[0.5,0.75,0.875],
            'dpara': {'user_data': [],
                      'algmethod': 1,
                      'maxf': 500+100*D,
                      'logfilename': '/dev/null'},
            'lpara': {'gtol': 0.00001,
                      'maxfun': 200}
            }

        self.stoppara = {'nmax': n}
        self.stopfn = gpbo.core.optimize.nstopfn
        self.reccfn = gpbo.core.reccomenders.gphinasrecc
        self.reccpara = {
            'ev': self.aqpara['ev'],
            'lb': self.aqpara['lb'],
            'ub': self.aqpara['ub'],
            'mprior': self.aqpara['mprior'],
            'sprior': self.aqpara['sprior'],
            'kindex': self.aqpara['kindex'],
            'maxf': 500+100*D, #10**(-min(12,max(6.,3*D))),
            'onlyafter': self.aqpara['nrandinit'],
            'check': True,
            'dpara':self.aqpara['dpara'],
            'lpara':self.aqpara['lpara'],
            'everyn': 1
        }
        self.ojfchar = {'dx': len(self.aqpara['lb']), 'dev': len(self.aqpara['ev'])}
        self.ojf=f
        self.path = path
        self.fname = fname
        return

class pesbslearns():
    def __init__(self,f,D,n,s,path,fname):
        self.aqfn = gpbo.core.acquisitions.PESbsaq
        self.aqpara = {
            'ev':{'s':s,'d':[sp.NaN],'xa':0.},
            'lb':[-1.] * D,
            'ub':[ 1.] * D,
            'maxf':500+100*D,
            'mprior': sp.array([1.]+[0.]*(D+1)+[-3]),
            'sprior': sp.array([1.]*(D+2)+[3]),
            'kindex':GPdc.MAT52CS,
            'DH_SAMPLES':16,
            'DM_SAMPLES':32,
            'DM_SUPPORT':1000,
            'SUPPORT_MODE':[gpbo.core.ESutils.SUPPORT_LAPAPROT],
            'DM_SLICELCBPARA':16.,
            'noS':False,
            'nrandinit':20,
            'cfn':lambda x,ev:42.,
            'traincfn':True,
            'xau':1.,
            'xal':0.,
            'startmode':'inline',
            'initpoints':[0.5,0.75,0.875],
            'dpara': {'user_data': [],
                      'algmethod': 1,
                      'maxf': 500+100*D,
                      'logfilename': '/dev/null'},
            'lpara': {'gtol': 0.00001,
                      'maxfun': 200}
            }

        self.stoppara = {'nmax': n}
        self.stopfn = gpbo.core.optimize.nstopfn
        self.reccfn = gpbo.core.reccomenders.gphinasrecc
        self.reccpara = {
            'ev': self.aqpara['ev'],
            'lb': self.aqpara['lb'],
            'ub': self.aqpara['ub'],
            'mprior': self.aqpara['mprior'],
            'sprior': self.aqpara['sprior'],
            'kindex': self.aqpara['kindex'],
            'maxf': 500+100*D,
            'onlyafter': self.aqpara['nrandinit'],
            'check': True,
            'dpara':self.aqpara['dpara'],
            'lpara':self.aqpara['lpara'],
            'everyn': 1
        }
        self.ojfchar = {'dx': len(self.aqpara['lb']), 'dev': len(self.aqpara['ev'])}
        self.ojf=f
        self.path = path
        self.fname = fname
        return

class switchdefault():
    """
        fixed s, space is [-1,1]^D

        """

    def __init__(self, f, D, ninit,nstop, s, path, fname):


        C = gpbo.core.config.pesfsdefault(f, D, 10, s, 'results', 'introspection.csv',ninit=ninit)
        aq0 = C.aqfn
        aq0para = C.aqpara

        aq1 = gpbo.core.acquisitions.splocalaq
        aq1para = {
            'ev': {'s': s, 'd': [sp.NaN]},
            'lb': [-1.] * D,
            'ub': [1.] * D,
            'start': [0.] * D
        }
        C2 = gpbo.core.config.eihypdefault(f, D, ninit, s, 'results', 'introspection.csv')
        aq2 = C2.aqfn
        aq2para = C2.aqpara

        self.chooser = gpbo.core.choosers.globallocalregret
        self.choosepara = {
            'ev': aq0para['ev'],
            'lb': aq0para['lb'],
            'ub': aq0para['ub'],
            'mprior': aq0para['mprior'],
            'sprior': aq0para['sprior'],
            'kindex': aq0para['kindex'],
            'nhyp' : aq0para['DH_SAMPLES'],
            'onlyafter': aq0para['nrandinit'],
            'check': True,
            'everyn': 1,
            'support': 2500,
            'draws': 25000,
            'regretswitch':1e-4,
            'dpara': {'user_data': [],
                      'algmethod': 1,
                      'maxf': 2000,
                      'logfilename': '/dev/null'},
            'lpara': {'gtol': 0.00001,
                      'maxfun': 400},
            'pvetol':1e-2,
            'lineSh':1e-6,
            'nlineS':150
        }

        self.aqfn = [aq0,aq1,aq2]
        self.aqpara = [aq0para,aq1para,aq2para]
        self.multimode = True

        self.stoppara = {'nmax': nstop}
        self.stopfn = gpbo.core.optimize.norlocalstopfn

        reccfn0 = C.reccfn
        reccpara0 = C.reccpara
        reccpara0['smode']='dthenl'
        reccfn1 = gpbo.core.reccomenders.argminrecc
        reccpara1 = {'check': True}

        self.reccfn = [reccfn0,reccfn1,reccfn0]
        self.reccpara = [reccpara0,reccpara1,reccpara0]

        self.ojfchar = {'dx': len(aq0para['lb']), 'dev': len(aq0para['ev'])}
        self.ojf = f

        self.path = path
        self.fname = fname
        return
class directdefault:
    def __init__(self, f, D, n, s, path, fname):
        self.aqfn = gpbo.core.acquisitions.directaq
        self.aqpara = {
            'ev': {'s': s, 'd': [sp.NaN]},
            'lb': [-1.] * D,
            'ub': [1.] * D
        }

        self.stoppara = {'nmax': n}
        self.stopfn = gpbo.core.optimize.norlocalstopfn

        self.reccfn = gpbo.core.reccomenders.argminrecc
        self.reccpara = {
            'ev': self.aqpara['ev'],
            'check': True
        }
        self.ojfchar = {'dx': len(self.aqpara['lb']), 'dev': len(self.aqpara['ev'])}
        self.ojf = f

        self.path = path
        self.fname = fname
        return

class cmaesdefault:
    def __init__(self, f, D, n, s, path, fname):
        self.aqfn = gpbo.core.acquisitions.cmaesaq
        self.aqpara = {
            'ev': {'s': s, 'd': [sp.NaN]},
            'lb': [-1.] * D,
            'ub': [1.] * D
        }

        self.stoppara = {'nmax': n}
        self.stopfn = gpbo.core.optimize.norlocalstopfn

        self.reccfn = gpbo.core.reccomenders.argminrecc
        self.reccpara = {
            'ev': self.aqpara['ev'],
            'check': True
        }
        self.ojfchar = {'dx': len(self.aqpara['lb']), 'dev': len(self.aqpara['ev'])}
        self.ojf = f

        self.path = path
        self.fname = fname
        return
