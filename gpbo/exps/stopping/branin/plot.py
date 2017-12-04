
import gpbo
rpath='results'
gpbo.figs.stoppingplots(rpath, ['switchingp_-2', 'switchingp_-4', 'pesfs_-2','pesfs_-4', 'eihyp_-10','eihyp_-12'], 16,fname='branin',title='Branin')
gpbo.figs.overhead(rpath, 'pesfs_-4',['switchingp_-2','switchingp_-4'], 16,fname='branin',title='',legendnames=['PES',['switching_2']])
