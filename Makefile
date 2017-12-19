UNAME    = $(shell uname)
EXE      = ntupleMerger
 
VPATH  = .:./
vpath %.h ./

CSUF   = cc
HSUF   = h
DICTC  = Dict.$(CSUF)
DICTH  = $(patsubst %.$(CSUF),%.h,$(DICTC))

SRCS   = Event.cc NtupleMerger.cc
OBJS   = $(patsubst %.$(CSUF), %.o, $(SRCS))

LDFLAGS  = -g
SOFLAGS  = -shared 
CXXFLAGS = -I./  

CXX       = g++
CXXFLAGS += -g -std=c++11
##-Wall -Wno-deprecated -std=c++11


HDRS_DICT = Hit.h Cluster.h Cbc.h Stub.h Track.h TelescopeEvent.h Event.h LinkDef.h

bin: ntupleMerger
all: 
	gmake cint 
	gmake bin 
cint: $(DICTC) 

$(DICTC): $(HDRS_DICT)
	@echo "Generating dictionary $(DICTC) with rootcling ..."	
	rootcling -f $@ -rmf NtupleMerger_xr.rootmap -c $^ 

ntupleMerger:  $(OBJS) Dict.o
	$(CXX) $(CXXFLAGS) `root-config --cflags` $(LDFLAGS) $^ -o $@ $(LIBS) `root-config --libs`

## Create object files
%.o : %.$(CSUF)
	$(CXX) $(CXXFLAGS) `root-config --cflags` -o $@ -c $<

## makedepend
depend: $(SRCS:.$(CSUF)=.$(CSUF).dep)
	@cat $(notdir $^) > Makefile.dep
	@-rm -f $(notdir $^) $(patsubst %,%.bak,$(notdir $^))
%.dep:
	@touch $(notdir $@)
	rmkdepend -f$(notdir $@) -- $(CXXFLAGS) `root-config --cflags` -- $*
include Makefile.dep

## Clean 
.PHONY   : clean 
clean : 
	@-rm $(OBJS) $(EXE) Dict_rdict.pcm $(DICTC) *.o
