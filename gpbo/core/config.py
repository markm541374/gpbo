import gpbo
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
            'mprior': sp.array([1.]+[0.]*D),
            'sprior': sp.array([1.]*(D+1)),
            'kindex': GPdc.MAT52,
            'volper':1e-6
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
                'volper':1e-6,
                'onlyafter':self.aqpara['nrandinit'],
                'check':True,
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
            'volper':1e-6
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
                'volper':1e-6,
                'onlyafter':self.aqpara['nrandinit'],
                'check':True,
                'everyn':1
                }
        self.ojfchar = {'dx': len(self.aqpara['lb']), 'dev': len(self.aqpara['ev'])}
        self.ojf=f

        self.path = path
        self.fname = fname
        return

class pesfsdefault():
    def __init__(self,f,D,n,s,path,fname):
        self.aqfn = gpbo.core.acquisitions.PESfsaq
        self.aqpara= {
            'ev': {'s': s, 'd': [sp.NaN]},
            'lb': [-1.]*D,
            'ub': [1.]*D,
            'nrandinit': 10,
            'volper': 1e-6,
            'mprior': sp.array([1.]+[0.]*D),
            'sprior': sp.array([1.]*(D+1)),
            'kindex': GPdc.MAT52,
            'DH_SAMPLES': 16,
            'DM_SAMPLES': 64,
            'DM_SUPPORT': 800,
            'SUPPORT_MODE': [gpbo.core.ESutils.SUPPORT_LAPAPROT],
            'DM_SLICELCBPARA': 1.,
            'noS': False,
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
                'volper':1e-6,
                'onlyafter':self.aqpara['nrandinit'],
                'check':True,
                'everyn':1
                }
        self.ojfchar = {'dx': len(self.aqpara['lb']), 'dev': len(self.aqpara['ev'])}
        self.ojf=f

        self.path = path
        self.fname = fname
        return

class pesfslearns():
    def __init__(self,f,D,n,s,path,fname):
        self.aqfn = gpbo.core.acquisitions.PESfsaq
        self.aqpara= {
            'ev': {'s': s, 'd': [sp.NaN]},
            'lb': [-1.]*D,
            'ub': [1.]*D,
            'nrandinit': 10,
            'volper': 1e-6,
            'mprior': sp.array([1.]+[0.]*D+[-2]),
            'sprior': sp.array([1.]*(D+1)+[3]),
            'kindex': GPdc.MAT52CS,
            'DH_SAMPLES': 16,
            'DM_SAMPLES': 16,
            'DM_SUPPORT': 800,
            'SUPPORT_MODE': [gpbo.core.ESutils.SUPPORT_LAPAPROT],
            'DM_SLICELCBPARA': 1.,
            'noS': False,
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
                'volper':1e-6,
                'onlyafter':self.aqpara['nrandinit'],
                'check':True,
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
            'volper': 1e-6,
            'mprior': sp.array([1.]+[0.]*D),
            'sprior': sp.array([1.]*(D+1)),
            'kindex': GPdc.MAT52,
            'DH_SAMPLES': 16,
            'DM_SAMPLES': 16,
            'DM_SUPPORT': 800,
            'SUPPORT_MODE': [gpbo.core.ESutils.SUPPORT_LAPAPROT],
            'DM_SLICELCBPARA': 1.,
            'noS': False,
            'logsu': lsu,
            'logsl': lsl,
            'cfn': cfn
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
                'volper':1e-6,
                'onlyafter':self.aqpara['nrandinit'],
                'check':True,
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
            'volper':1e-6,
            'mprior': sp.array([1.]+[0.]*(D+1)),
            'sprior': sp.array([1.]*(D+2)),
            'kindex':GPdc.MAT52,
            'DH_SAMPLES':16,
            'DM_SAMPLES':16,
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
            'initpoints':[0.5,0.75,0.875]
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
            'volper': 1e-6,
            'onlyafter': self.aqpara['nrandinit'],
            'check': True,
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
            'volper':1e-6,
            'mprior': sp.array([1.]+[0.]*(D+1)+[-3]),
            'sprior': sp.array([1.]*(D+2)+[3]),
            'kindex':GPdc.MAT52CS,
            'DH_SAMPLES':16,
            'DM_SAMPLES':16,
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
            'initpoints':[0.5,0.75,0.875]
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
            'volper': 1e-6,
            'onlyafter': self.aqpara['nrandinit'],
            'check': True,
            'everyn': 1
        }
        self.ojfchar = {'dx': len(self.aqpara['lb']), 'dev': len(self.aqpara['ev'])}
        self.ojf=f
        self.path = path
        self.fname = fname
        return