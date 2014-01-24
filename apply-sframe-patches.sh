#!/bin/bash

cd $SFRAME_DIR || { echo "ERROR: SFRAME is not setup up"; exit 1; }

function patches_already_applied() {
   grep BOOSTDIR Makefile.common > /dev/null && return 1;
   return 0;
}

function apply_patches() {
  for f in SFrameTools/sframe-patches/*.patch; do
    echo -n Applying $f ...
    patch -s -f -p0 -i $f || { echo "error"; exit 1; }
    echo "done"
  done
}

patches_already_applied
if [ $? -eq 1 ]; then
echo "It looks like the patches have been applied already, not patching again.";
else
apply_patches;
fi


# also write default BOOSTDIR in fullsetup.sh, if noit there already:
grep 'export BOOSTDIR=' fullsetup.sh > /dev/null
if [ $? -eq 1 ]; then
   echo -e '\nexport BOOSTDIR=`cd $CMSSW_BASE; scram tool tag boost include' >> fullsetup.sh
   echo "wrote BOOSTDIR to fullsetup.sh"
fi


