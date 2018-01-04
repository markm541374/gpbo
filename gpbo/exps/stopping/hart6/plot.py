#import gpbo
#rpath='results'
#this has an offset!!!!
#gpbo.figs.stoppingplots(rpath, ['switchingp_-2', 'switchingp_-4', 'pesfs_-2','pesfs_-4', 'eihyp_-10','eihyp_-12'], 10,fname='hart6',title='Hartmann 6')
#gpbo.figs.overhead(rpath, 'pesfs_-4',['switchingp_-2','switchingp_-4'], 10,fname='hart6',title='',legendnames=['PES',['switching_2']])

import gpbo
print('plotting hart6 stopping')
rpath='results'
fpath='figs'
gpbo.figs.stoppingplots(rpath, ['switchingp_-2', 'switchingp_-4', 'pesfs_-2','pesfs_-4', 'eihyp_-10','eihyp_-12'], 10,legendnames= ['Switching: $R_{global}=10^{-2}$', 'Switching: $R_{global}=10^{-4}$', 'PES: $AQ_{stop} = 10^{-2}$','PES: $AQ_{stop}=10^{-4}$', 'EI: $PI_{stop}=10^{-10}$','EI: $PI_{stop}=10^{-12}$'],fname='hart6',title='Hartmann 6',fpath=fpath,showlegend=True)

#gpbo.figs.overhead(rpath, 'pesfs_-4',['switchingp_-2','switchingp_-4'], 16,fname='hart6',title='',legendnames=['PES',['switching 2']],fpath=fpath)
