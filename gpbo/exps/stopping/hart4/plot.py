
#import gpbo
#rpath='results'
#this has an offset!!!!
#gpbo.figs.stoppingplots(rpath, ['switchingp_-2', 'switchingp_-4', 'pesfs_-2','pesfs_-4', 'eihyp_-10','eihyp_-12'], 16,fname='hart4',title='Hartmann 4')
#gpbo.figs.overhead(rpath, 'pesfs_-24',['switchingp_-2','switchingp_-4'], 16,fname='hart4',title='',legendnames=['PES',['switching_2']])

import gpbo
print('plotting hart4 stopping')
rpath='results'
fpath='figs'
gpbo.figs.stoppingplots(rpath, ['switchingp_-2', 'switchingp_-4', 'pesfs_-2','pesfs_-4', 'eihyp_-10','eihyp_-12'], 16,legendnames= ['BLOSSOM: $R_{global}=10^{-2}$', 'BLOSSOM: $R_{global}=10^{-4}$', 'PES: $AQ_{stop} = 10^{-2}$','PES: $AQ_{stop}=10^{-4}$', 'EI: $PI_{stop}=10^{-10}$','EI: $PI_{stop}=10^{-12}$'],fname='hart4',title='Hartmann 4',fpath=fpath)

#gpbo.figs.overhead(rpath, 'pesfs_-24',['switchingp_-2','switchingp_-4'], 16,fname='hart4',title='',legendnames=['PES',['switching 2']],fpath=fpath)
