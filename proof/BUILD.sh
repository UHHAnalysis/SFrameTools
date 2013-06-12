# $Id: BUILD.sh,v 1.1 2012/05/25 08:39:23 rkogler Exp $

if [ "$1" = "clean" ]; then
    make distclean
    exit 0
fi

if [ "x$ROOTPROOFLITE" != "x" ]; then
    echo "Running on PROOF-Lite, skipping build"
    exit 0
fi

make default
