#!/usr/bin/env bash

(cd gpbo/core;
echo "rebuilding gpbo.core..."
for FILE in *.pyx; do
    NAME=$(echo $FILE | cut -f 1 -d '.');
    if (( $(date +%s -r $FILE) > $(date +%s -r $"$NAME.c") ))
    then
        echo "rebuilding $NAME"
        cython $FILE
        gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing -I/usr/include/python2.7 $"-I$HOME/.local/lib/python2.7/site-packages/numpy/core/include" -o $"$NAME.so" $"$NAME.c"

    else
        echo "$NAME unchanged"
    fi
    done

)
echo "done"