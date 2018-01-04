
#import gpbo
#rpath='results'
##gpbo.figs.stoppingplots(rpath, ['switchingp_-2', 'switchingp_-4', 'pesfs_-2','pesfs_-4', 'eihyp_-10','eihyp_-12'], 12,fname='camel6',title='Camel 6')
#gpbo.figs.overhead(rpath, 'pesfs_-24',['switchingp_-2','switchingp_-4'], 16,fname='camel6',title='Camel 6')
import gpbo

print('plotting camel6 stopping')
rpath='results'
fpath='figs'
gpbo.figs.stoppingplots(rpath, ['switchingp_-2', 'switchingp_-4', 'pesfs_-2','pesfs_-4', 'eihyp_-10','eihyp_-12'], 16,legendnames= ['Switching: $R_{global}=10^{-2}$', 'Switching: $R_{global}=10^{-4}$', 'PES: $AQ_{stop} = 10^{-2}$','PES: $AQ_{stop}=10^{-4}$', 'EI: $PI_{stop}=10^{-10}$','EI: $PI_{stop}=10^{-12}$'],fname='camel6',title='Camel 6',fpath=fpath)

#gpbo.figs.overhead(rpath, 'pesfs_-24',['switchingp_-2','switchingp_-4'], 16,fname='camel6',title='',legendnames=['PES',['switching 2']],fpath=fpath)
