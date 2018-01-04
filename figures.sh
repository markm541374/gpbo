#!/usr/bin/env bash
export PYTHONPATH="/home/mark/Dropbox/workspace/gpbo:${PYTHONPATH}"

#implementation plots
if false; then
    cp matplotlibrc gpbo/exps/adaptive/
    (cd gpbo/exps/adaptive ; python2 adaptive.py)
    cp gpbo/exps/adaptive/figs/overheadtime.pdf figures/implementation/
    cp gpbo/exps/adaptive/figs/opcount.pdf figures/implementation/
    cp gpbo/exps/adaptive/figs/useful.pdf figures/implementation/
    cp gpbo/exps/adaptive/figs/updateorder.pdf figures/implementation/
fi
#acquisition search parallel
if false; then
    cp matplotlibrc gpbo/exps/demo/
    (cd gpbo/exps/demo ; python2 eitime.py)
    cp gpbo/exps/demo/figs/timing.pdf figures/implementation/
fi
if $1; then
    cp matplotlibrc gpbo/exps/stopping/hart3
    (cd gpbo/exps/stopping/hart3 ; python2 plot.py)
    cp gpbo/exps/stopping/hart3/figs/stopping_hart3.pdf figures/localstop/

    cp matplotlibrc gpbo/exps/stopping/hart4
    (cd gpbo/exps/stopping/hart4 ; python2 plot.py)
    cp gpbo/exps/stopping/hart4/figs/stopping_hart4.pdf figures/localstop/

    cp matplotlibrc gpbo/exps/stopping/hart6
    (cd gpbo/exps/stopping/hart6 ; python2 plot.py)
    cp gpbo/exps/stopping/hart6/figs/stopping_hart6.pdf figures/localstop/
fi
if $1; then
    cp matplotlibrc gpbo/exps/stopping/branin
    (cd gpbo/exps/stopping/branin ; python2 plot.py)
    cp gpbo/exps/stopping/branin/figs/stopping_branin.pdf figures/localstop/

    cp matplotlibrc gpbo/exps/stopping/camel3
    (cd gpbo/exps/stopping/camel3 ; python2 plot.py)
    cp gpbo/exps/stopping/camel3/figs/stopping_camel3.pdf figures/localstop/

    cp matplotlibrc gpbo/exps/stopping/camel6
    (cd gpbo/exps/stopping/camel6 ; python2 plot.py)
    cp gpbo/exps/stopping/camel6/figs/stopping_camel6.pdf figures/localstop/
fi
