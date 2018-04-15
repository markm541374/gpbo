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
#stopping performance plots
if false; then
    cp matplotlibrc gpbo/exps/stopping/hart3
    (cd gpbo/exps/stopping/hart3 ; python2 plot.py)
    cp gpbo/exps/stopping/hart3/figs/stopping_hart3.pdf figures/localstop/

    cp matplotlibrc gpbo/exps/stopping/hart4
    (cd gpbo/exps/stopping/hart4 ; python2 plot.py)
    cp gpbo/exps/stopping/hart4/figs/stopping_hart4.pdf figures/localstop/

    cp matplotlibrc gpbo/exps/stopping/camel3
    (cd gpbo/exps/stopping/camel3 ; python2 plot.py)
    cp gpbo/exps/stopping/camel3/figs/stopping_camel3.pdf figures/localstop/

    cp matplotlibrc gpbo/exps/stopping/camel6
    (cd gpbo/exps/stopping/camel6 ; python2 plot.py)
    cp gpbo/exps/stopping/camel6/figs/stopping_camel6.pdf figures/localstop/
fi
if false; then
    cp matplotlibrc gpbo/exps/stopping/hart6
    (cd gpbo/exps/stopping/hart6 ; python2 plot.py)
    cp gpbo/exps/stopping/hart6/figs/stopping_hart6.pdf figures/localstop/
fi
#branin separater from teh others because I aslo want the overhead plot
if false; then
    cp matplotlibrc gpbo/exps/stopping/branin
    (cd gpbo/exps/stopping/branin ; python2 plot.py)
    cp gpbo/exps/stopping/branin/figs/stopping_branin.pdf figures/localstop/
    cp gpbo/exps/stopping/branin/figs/overhead_branin.pdf figures/localstop/

fi
if false; then
    cp matplotlibrc gpbo/exps/stopping/branin
    (cd gpbo/exps/stopping/branin ; python2 regretbounds.py)
    cp gpbo/exps/stopping/branin/figs/ends.pdf figures/localstop/

fi
if false; then
    cp matplotlibrc gpbo/exps/stopping/draws
    (cd gpbo/exps/stopping/draws ; python2 exp0.py)
    cp gpbo/exps/stopping/draws/results2d/out11.pdf figures/localstop/
    mv figures/localstop/out11.pdf figures/localstop/d2methods.pdf
fi
if false; then
    cp matplotlibrc gpbo/exps/stopping/gphyp
    (cd gpbo/exps/stopping/gphyp ; python2 plot.py)
    cp gpbo/exps/stopping/gphyp/figs/stopping_GPhyp.pdf figures/localstop/
fi
#theory chap demo figs for gp and bayesopt
if false; then
    cp matplotlibrc gpbo/exps/demo
    (cd gpbo/exps/demo ; python2 gpfigs.py)
    cp gpbo/exps/demo/figs/gpdemo.pdf figures/theory/
fi
if false; then
    cp matplotlibrc gpbo/exps/demo
    (cd gpbo/exps/demo ; python2 bayesoptfig.py)
    cp gpbo/exps/demo/figs/bayesoptdemo.pdf figures/theory/
fi
if false; then
    cp matplotlibrc gpbo/exps/demo
    (cd gpbo/exps/demo ; python2 stoppingdemo.py)
    cp gpbo/exps/demo/figs/stoppingdemo.pdf figures/localstop/
fi
if false; then
    cp matplotlibrc gpbo/exps/stopping/hart4/noisedemo
    (cd gpbo/exps/stopping/hart4/noisedemo ; python2 jitterplot.py)
    cp gpbo/exps/stopping/hart4/noisedemo/results/out11.pdf figures/localstop/
    cp gpbo/exps/stopping/hart4/noisedemo/results/out19.pdf figures/localstop/
    mv figures/localstop/out11.pdf figures/localstop/conditionregret.pdf
    mv figures/localstop/out19.pdf figures/localstop/conditionmagnitude.pdf
fi
if false; then
    cp matplotlibrc gpbo/exps/stopping/hart4/noisedemo
    (cd gpbo/exps/stopping/hart4/noisedemo ; python2 jitterplotthin.py)
    cp gpbo/exps/stopping/hart4/noisedemo/results/out11.pdf figures/localstop/
    cp gpbo/exps/stopping/hart4/noisedemo/results/out19.pdf figures/localstop/
    mv figures/localstop/out11.pdf figures/localstop/conditionregretthin.pdf
    mv figures/localstop/out19.pdf figures/localstop/conditionmagnitudethin.pdf
fi
if false; then
    cp matplotlibrc gpbo/exps/stopping/draws
    (cd gpbo/exps/stopping/draws ; python2 meanregret.py)
    cp gpbo/exps/stopping/draws/results2d/out20.pdf figures/localstop/
    mv figures/localstop/out20.pdf figures/localstop/draw2mean.pdf
fi

if false; then
    cp matplotlibrc gpbo/exps/biasopt/draw2d
    (cd gpbo/exps/biasopt/draw2d ; python2 expicmlF3.py)
    cp gpbo/exps/biasopt/draw2d/F3new/out13.pdf figures/variablefidelity/
    mv figures/variablefidelity/out13.pdf figures/variablefidelity/badev.pdf
    (cd gpbo/exps/biasopt/draw2d ; python2 expicmlF4.py)
    cp gpbo/exps/biasopt/draw2d/F4new/out13.pdf figures/variablefidelity/
    mv figures/variablefidelity/out13.pdf figures/variablefidelity/goodev.pdf
fi
if false; then
    cp matplotlibrc gpbo/exps/biasopt/offsetbranin
    (cd gpbo/exps/biasopt/offsetbranin ; python2 expicmlF5a.py)
    cp gpbo/exps/biasopt/offsetbranin/F5anew/out13.pdf figures/variablefidelity/
    mv figures/variablefidelity/out13.pdf figures/variablefidelity/branin.pdf
    cp gpbo/exps/biasopt/offsetbranin/F5anew/out12.pdf figures/variablefidelity/
    mv figures/variablefidelity/out12.pdf figures/variablefidelity/branincost.pdf
    (cd gpbo/exps/biasopt/offsetbranin ; python2 expicmlF1.py)
    cp gpbo/exps/biasopt/offsetbranin/F1new/out11.pdf figures/variablefidelity/
    mv figures/variablefidelity/out11.pdf figures/variablefidelity/argpost.pdf
    (cd gpbo/exps/biasopt/offsetbranin ; python2 expicmlF2.py)
    cp gpbo/exps/biasopt/offsetbranin/icmlF2/out11.pdf figures/variablefidelity/
    cp gpbo/exps/biasopt/offsetbranin/icmlF2/out15.pdf figures/variablefidelity/
    mv figures/variablefidelity/out11.pdf figures/variablefidelity/ofbcost0.pdf
    mv figures/variablefidelity/out15.pdf figures/variablefidelity/ofbcost1.pdf

fi
if false; then
    cp matplotlibrc gpbo/exps/biasopt/hartmann3
    (cd gpbo/exps/biasopt/hartmann3 ; python2 expicmlF5b.py)
    cp gpbo/exps/biasopt/hartmann3/F5bnew/out13.pdf figures/variablefidelity/
    mv figures/variablefidelity/out13.pdf figures/variablefidelity/hart3.pdf
    cp gpbo/exps/biasopt/hartmann3/F5bnew/out12.pdf figures/variablefidelity/
    mv figures/variablefidelity/out12.pdf figures/variablefidelity/hart3cost.pdf
fi
if false; then
    cp matplotlibrc gpbo/exps/biasopt/hartmann6
    (cd gpbo/exps/biasopt/hartmann6 ; python2 expicmlF5c.py)
    cp gpbo/exps/biasopt/hartmann6/F5cnew/out13.pdf figures/variablefidelity/
    mv figures/variablefidelity/out13.pdf figures/variablefidelity/hart6.pdf
    cp gpbo/exps/biasopt/hartmann6/F5cnew/out12.pdf figures/variablefidelity/
    mv figures/variablefidelity/out12.pdf figures/variablefidelity/hart6cost.pdf
fi
if $1; then
    cp matplotlibrc gpbo/exps/biasopt/mnist
    (cd gpbo/exps/biasopt/mnist ; python2 expicmlF6.py)
    cp gpbo/exps/biasopt/mnist/icmlF6/out13.pdf figures/variablefidelity/mnist.pdf
    cp gpbo/exps/biasopt/mnist/icmlF6/out12.pdf figures/variablefidelity/mnistev.pdf
fi
if false; then
    cp matplotlibrc gpbo/exps/biasopt/powseries
    (cd gpbo/exps/biasopt/powseries ; python2 expicmlF7.py)
    cp gpbo/exps/biasopt/powseries/norm/out6.pdf figures/variablefidelity/powerfit.pdf
    cp gpbo/exps/biasopt/powseries/norm/out5.pdf figures/variablefidelity/powerfitev.pdf
fi
if false; then
    cp matplotlibrc gpbo/exps/predictive/overhead/fixedPES
    (cd gpbo/exps/predictive/overhead/fixedPES ; python2 plots.py)
    cp gpbo/exps/predictive/overhead/fixedPES/figs/iterpes.pdf figures/predictive/
    cp gpbo/exps/predictive/overhead/fixedPES/figs/evcostpes.pdf figures/predictive/
    cp gpbo/exps/predictive/overhead/fixedPES/figs/aqcostpes.pdf figures/predictive/
fi
if false; then
    cp matplotlibrc gpbo/exps/predictive/overhead/fixedEI
    (cd gpbo/exps/predictive/overhead/fixedEI ; python2 plots.py)
    cp gpbo/exps/predictive/overhead/fixedEI/figs/iterei.pdf figures/predictive/
    cp gpbo/exps/predictive/overhead/fixedEI/figs/evcostei.pdf figures/predictive/
    cp gpbo/exps/predictive/overhead/fixedEI/figs/aqcostei.pdf figures/predictive/
fi
if false; then
    cp matplotlibrc gpbo/exps/predictive/
    (cd gpbo/exps/predictive ; python2 plotoverhead.py)
    cp gpbo/exps/predictive/plotfigs/overheadcum.pdf figures/predictive
    cp gpbo/exps/predictive/plotfigs/overheadsingle.pdf figures/predictive
fi
if false; then
    cp matplotlibrc gpbo/exps/predictive/prediction
    (cd gpbo/exps/predictive/prediction ; python2 predplots.py)
    cp gpbo/exps/predictive/prediction/figs/margpredictions.pdf figures/predictive
fi
if false; then
    cp matplotlibrc gpbo/exps/predictive/performance
    (cd gpbo/exps/predictive/performance ; python2 plots.py)
    cp gpbo/exps/predictive/performance/figs/evcostei_300.pdf figures/predictive/perfmodelshort.pdf
    cp gpbo/exps/predictive/performance/figs/evcostei_800.pdf figures/predictive/perfmodelmedium.pdf
    cp gpbo/exps/predictive/performance/figs/evcostei_1500.pdf figures/predictive/perfmodellong.pdf
fi
if false; then
    cp matplotlibrc gpbo/exps/predictive/prediction
    (cd gpbo/exps/predictive/prediction ; python2 plotcontour.py)
    cp gpbo/exps/predictive/prediction/figs/Rmeancontour.pdf figures/predictive/Rmeancontour.pdf
    cp gpbo/exps/predictive/prediction/figs/obsvarcontour.pdf figures/predictive/obsvarcontour.pdf
    cp gpbo/exps/predictive/prediction/figs/Estepscontour.pdf figures/predictive/Estepscontour.pdf
    cp gpbo/exps/predictive/prediction/figs/overheadcontour.pdf figures/predictive/overheadcontour.pdf
fi
if false; then
    cp matplotlibrc gpbo/exps/predictive/prediction
    (cd gpbo/exps/predictive/prediction ; python2 predplots.py)
    cp gpbo/exps/predictive/prediction/figs/margpredictions.pdf figures/predictive/margpredictions.pdf
fi
if false; then
    cp matplotlibrc gpbo/exps/biasopt/denmark
    (cd gpbo/exps/biasopt/denmark ; python2 plotfigs.py)
    cp gpbo/exps/biasopt/denmark/results/out5.pdf figures/variablefidelity/denmarkev.pdf
    cp gpbo/exps/biasopt/denmark/results/out6.pdf figures/variablefidelity/denmarkfc.pdf
fi
if false; then
    cp matplotlibrc gpbo/exps/biasopt/casp
    (cd gpbo/exps/biasopt/casp ; python2 plotfigs.py)
    cp gpbo/exps/biasopt/casp/results/out5.pdf figures/variablefidelity/caspev.pdf
    cp gpbo/exps/biasopt/casp/results/out6.pdf figures/variablefidelity/caspfc.pdf
fi
if false; then
    cp matplotlibrc gpbo/exps/biasopt/cciso
    (cd gpbo/exps/biasopt/cciso ; python2 plotfigs.py)
    cp gpbo/exps/biasopt/cciso/results/out12.pdf figures/variablefidelity/ccisoev.pdf
    cp gpbo/exps/biasopt/cciso/results/out13.pdf figures/variablefidelity/ccisofc.pdf
fi
if false; then
    cp matplotlibrc gpbo/exps/biasopt/ccard
    (cd gpbo/exps/biasopt/ccard ; python2 plotfigs.py)
    cp gpbo/exps/biasopt/ccard/results/out12.pdf figures/variablefidelity/ccardev.pdf
    cp gpbo/exps/biasopt/ccard/results/out13.pdf figures/variablefidelity/ccardfc.pdf
fi
