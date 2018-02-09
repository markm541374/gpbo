
#import gpbo
#rpath='results'
#gpbo.figs.stoppingplots(rpath, ['switchingp_-2', 'switchingp_-4', 'pesfs_-2','pesfs_-4', 'eihyp_-10','eihyp_-12'], 16,fname='branin',title='Branin')
#gpbo.figs.overhead(rpath, 'pesfs_-4',['switchingp_-2','switchingp_-4'], 16,fname='branin',title='',legendnames=['PES',['switching_2']])
import gpbo
print('plotting gphyp stopping')
rpath='results'
fpath='figs'
gpbo.figs.stoppingplots(rpath, ['switchingp_-2', 'switchingp_-4', 'pesfs_-2','pesfs_-4', 'eihyp_-10','eihyp_-12'], 16,offset=10.0114506889,legendnames= ['BLOSSOM: $R_{global}=10^{-2}$', 'BLOSSOM: $R_{global}=10^{-4}$', 'PES: $AQ_{stop} = 10^{-2}$','PES: $AQ_{stop}=10^{-4}$', 'EI: $PI_{stop}=10^{-10}$','EI: $PI_{stop}=10^{-12}$'],fname='GPhyp',title='GP hyperparameters',fpath=fpath,logy=True,showlegend=True)

#gpbo.figs.overhead(rpath, 'pesfs_-4',['switchingp_-2','switchingp_-4'], 16,fname='branin',title='',legendnames=['PES',['Switching: $R_{global}=10^{-2}$']],fpath=fpath)
