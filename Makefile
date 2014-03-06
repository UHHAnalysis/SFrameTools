dummy := $(shell ./apply-sframe-patches.sh > /dev/null)

# Package information
LIBRARY = SFrameTools
OBJDIR  = obj
DEPDIR  = $(OBJDIR)/dep
SRCDIR  = src
INCDIR  = include

USERLDFLAGS += $(shell root-config --libs) -lGenVector

# Include definitions
include Makefile.defs
