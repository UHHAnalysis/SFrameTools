#!/bin/bash

cd $SFRAME_DIR || { echo "ERROR: SFRAME is not setup up"; exit 1; }

svnversion | grep M > /dev/null

if [ $? -eq 0 ]; then
echo "ERROR: svnversion says SFrame is already modified. Patches have NOT been applied; please apply patches manually.";
exit 1;
fi

for f in SFrameTools/sframe-patches/*.patch; do
echo -n Applying $f ...
patch -s -f -p0 -i $f || { echo "error"; exit 1; }
echo "done"
done

