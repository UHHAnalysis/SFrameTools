# Package information
LIBRARY = SFrameTools
OBJDIR  = obj
DEPDIR  = $(OBJDIR)/dep
SRCDIR  = src
INCDIR  = include

# configure FastJet
FASTJETDIR = /afs/naf.desy.de/user/p/peiffer/fastjet-3.0.2/install

INCLUDES += -I$(FASTJETDIR)/include

# Include the generic compilation rules
include $(SFRAME_DIR)/Makefile.common
