#!/usr/bin/env bash

(cd gpbo/core;
echo "rebuilding gpbo.core..."
for FILE in *.pyx; do
    NAME=$(echo $FILE | cut -f 1 -d '.');
    A=$(( $(date +%s -r $FILE) > $(date +%s -r $"$NAME.c") ))
    #A=1
    if $1; then #[ "$A" -eq 1 ]; then
        echo "rebuilding $NAME"
        $HOME/.pyenv/shims/cython  $FILE
        gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing -I$HOME/.pyenv/versions/2.7.13/include/python2.7 $"-I$HOME/.pyenv/versions/2.7.13/lib/python2.7/site-packages/numpy/core/include" -o $"$NAME.so" $"$NAME.c"
    else
        echo "$NAME unchanged"
    fi
    done

)
echo "done"