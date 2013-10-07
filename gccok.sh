#!/bin/bash

# script to test whether gcc version is >= 46X. This is used in the SFrameTools / SFrameAnalysis Makefiles
# to warn users about old compilers.

# You can also execute this script directly. If the gcc version is Ok (>= 4.6.X), it should just print "yes" (note that this is the condition checked in the Makefile,
# so don't modify the output in that case). If it is not Ok, it will print some details about the detected version and path.

majorgcc=`gcc -dM -E - < /dev/null | grep __GNUC__ | cut -f3 -d" "`
minorgcc=`gcc -dM -E - < /dev/null | grep __GNUC_MINOR__ | cut -f3 -d" "`


# we require gcc >= 46X. If this is the case, write "yes" to stdout and exit. If not, output some more
# info for the user explaining what is going on ...

if [ \( "$majorgcc" -gt 4 \) -o \( "$majorgcc" -eq 4 -a "$minorgcc" -ge 6 \) ]; then
  echo "yes"
  exit 0
else
  echo "The compiler is too old: The detected gcc version is ${majorgcc}.${minorgcc} (path to gcc: `which gcc`), but required is >= 4.6.X";
  exit 1;
fi

