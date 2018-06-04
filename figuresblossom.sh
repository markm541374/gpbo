#!/usr/bin/env bash
export PYTHONPATH="/home/mark/Dropbox/workspace/gpbo:${PYTHONPATH}"

#stopping performance plots
if false; then
    cp matplotlibrcblossom gpbo/exps/stopping/hart3/matplotlibrc
    (cd gpbo/exps/stopping/hart3 ; python2 plot.py)
    cp gpbo/exps/stopping/hart3/figs/stopping_hart3.pdf blossomfigs/localstop/

    cp matplotlibrcblossom gpbo/exps/stopping/hart4/matplotlibrc
    (cd gpbo/exps/stopping/hart4 ; python2 plot.py)
    cp gpbo/exps/stopping/hart4/figs/stopping_hart4.pdf blossomfigs/localstop/

    cp matplotlibrcblossom gpbo/exps/stopping/camel3/matplotlibrc
    (cd gpbo/exps/stopping/camel3 ; python2 plot.py)
    cp gpbo/exps/stopping/camel3/figs/stopping_camel3.pdf blossomfigs/localstop/

    cp matplotlibrcblossom gpbo/exps/stopping/camel6/matplotlibrc
    (cd gpbo/exps/stopping/camel6 ; python2 plot.py)
    cp gpbo/exps/stopping/camel6/figs/stopping_camel6.pdf blossomfigs/localstop/
fi
if false; then
    cp matplotlibrcblossom gpbo/exps/stopping/hart6/matplotlibrc
    (cd gpbo/exps/stopping/hart6 ; python2 plot.py)
    cp gpbo/exps/stopping/hart6/figs/stopping_hart6.pdf blossomfigs/localstop/
fi
#branin separater from teh others because I aslo want the overhead plot
if false; then
    cp matplotlibrcblossom gpbo/exps/stopping/branin/matplotlibrc
    (cd gpbo/exps/stopping/branin ; python2 plot.py)
    cp gpbo/exps/stopping/branin/figs/stopping_branin.pdf blossomfigs/localstop/
    cp gpbo/exps/stopping/branin/figs/overhead_branin.pdf blossomfigs/localstop/

fi
if false; then
    cp matplotlibrcblossom gpbo/exps/stopping/branin/matplotlibrc
    (cd gpbo/exps/stopping/branin ; python2 regretbounds.py)
    cp gpbo/exps/stopping/branin/figs/ends.pdf blossomfigs/localstop/

fi
if false; then
    cp matplotlibrcblossom gpbo/exps/stopping/draws/matplotlibrc
    (cd gpbo/exps/stopping/draws ; python2 exp0.py)
    cp gpbo/exps/stopping/draws/results2d/out11.pdf blossomfigs/localstop/
    mv figures/localstop/out11.pdf blossomfigs/localstop/d2methods.pdf
fi
if $1; then
    cp matplotlibrcblossom gpbo/exps/stopping/gphyp/matplotlibrc
    (cd gpbo/exps/stopping/gphyp ; python2 plot.py)
    cp gpbo/exps/stopping/gphyp/figs/stopping_GPhyp.pdf blossomfigs/localstop/
fi
#theory chap demo figs for gp and bayesopt
if false; then
    cp matplotlibrcblossom gpbo/exps/demo/matplotlibrc
    (cd gpbo/exps/demo ; python2 gpfigs.py)
    cp gpbo/exps/demo/figs/gpdemo.pdf figures/theory/
fi
if false; then
    cp matplotlibrcblossom gpbo/exps/demo/matplotlibrc
    (cd gpbo/exps/demo ; python2 bayesoptfig.py)
    cp gpbo/exps/demo/figs/bayesoptdemo.pdf figures/theory/
fi
if false; then
    cp matplotlibrcblossom gpbo/exps/demo/matplotlibrc
    (cd gpbo/exps/demo ; python2 stoppingdemo.py)
    cp gpbo/exps/demo/figs/stoppingdemo.pdf blossomfigs/localstop/
fi
if false; then
    cp matplotlibrcblossom gpbo/exps/stopping/hart4/noisedemo/matplotlibrc
    (cd gpbo/exps/stopping/hart4/noisedemo ; python2 jitterplot.py)
    cp gpbo/exps/stopping/hart4/noisedemo/results/out11.pdf blossomfigs/localstop/
    cp gpbo/exps/stopping/hart4/noisedemo/results/out19.pdf blossomfigs/localstop/
    mv figures/localstop/out11.pdf blossomfigs/localstop/conditionregret.pdf
    mv figures/localstop/out19.pdf blossomfigs/localstop/conditionmagnitude.pdf
fi
if false; then
    cp matplotlibrcblossom gpbo/exps/stopping/hart4/noisedemo/matplotlibrc
    (cd gpbo/exps/stopping/hart4/noisedemo ; python2 jitterplotthin.py)
    cp gpbo/exps/stopping/hart4/noisedemo/results/out11.pdf blossomfigs/localstop/
    cp gpbo/exps/stopping/hart4/noisedemo/results/out19.pdf blossomfigs/localstop/
    mv figures/localstop/out11.pdf blossomfigs/localstop/conditionregretthin.pdf
    mv figures/localstop/out19.pdf blossomfigs/localstop/conditionmagnitudethin.pdf
fi
if false; then
    cp matplotlibrcblossom gpbo/exps/stopping/draws/matplotlibrc
    (cd gpbo/exps/stopping/draws ; python2 meanregret.py)
    cp gpbo/exps/stopping/draws/results2d/out20.pdf blossomfigs/localstop/
    mv figures/localstop/out20.pdf blossomfigs/localstop/draw2mean.pdf
fi
