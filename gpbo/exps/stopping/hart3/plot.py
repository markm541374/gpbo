
import gpbo
rpath='results'
#this has an offset!!!!
#gpbo.figs.stoppingplots(rpath, ['switchingp_-2', 'switchingp_-4', 'pesfs_-2','pesfs_-4', 'eihyp_-10','eihyp_-12'], 16,fname='hart3',title='Hartmann 3')
gpbo.figs.overhead(rpath, 'pesfs_-4',['switchingp_-2','switchingp_-4'], 16,fname='hart3',title='',legendnames=['PES',['switching_2']])
