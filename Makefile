dummy := $(shell ./apply-sframe-patches.sh > /dev/null)

# Package information
LIBRARY = SFrameTools
OBJDIR  = obj
DEPDIR  = $(OBJDIR)/dep
SRCDIR  = src
INCDIR  = include

# configure FastJet
INCLUDES += -I$(FASTJETDIR)/../include

USERCXXFLAGS := -g

# note: this is just to make it compile cleanly. If you want to perform pdf studies, you have to
# install lhapdf and provide the actual path here!
INCLUDES += -I/afs/cern.ch/sw/lcg/external/MCGenerators/lhapdf/5.8.8/x86_64-slc5-gcc46-opt/include

# Include the generic compilation rules
include $(SFRAME_DIR)/Makefile.common
