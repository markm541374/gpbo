
import gpbo
rpath='results'
#this has an offset!!!!
#gpbo.figs.stoppingplots(rpath, ['switchingp_-2', 'switchingp_-4', 'pesfs_-2','pesfs_-4', 'eihyp_-10','eihyp_-12'], 16,fname='hart4',title='Hartmann 4')
gpbo.figs.overhead(rpath, 'pesfs_-24',['switchingp_-2','switchingp_-4'], 16,fname='hart4',title='',legendnames=['PES',['switching_2']])
