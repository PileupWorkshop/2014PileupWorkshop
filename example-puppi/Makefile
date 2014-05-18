
# Makefile generated automatically by ../scripts/mkcxx.pl '-f' '-i' '-I../Framework' '-l' '-L../Framework -lEventMixer -lz' '-r'
# run 'make make' to update it if you add new files

CXX = g++  # for macs - otherwise get c++ = clang
CXXFLAGS = -Wall -g -O2

# also arrange for fortran support
FC = gfortran
FFLAGS = -Wall -O2

FJCONFIG = fastjet-config
INCLUDE += `$(FJCONFIG) --cxxflags`
LIBRARIES  += `$(FJCONFIG) --libs --plugins`

INCLUDE += `root-config --cflags`
LIBRARIES  += `root-config --libs`
LIBRARIES  += -lJetCleanser -lGenericSubtractor -lSoftKiller
INCLUDE += $(LCLINCLUDE)

COMMONSRC = puppiContainer.cc
F77SRC = 
COMMONOBJ = puppiContainer.o

PROGSRC = example02.cc
PROGOBJ = example02.o

INCLUDE += -I../Framework
LIBRARIES += -L../Framework -lEventMixer -lz


all:  example02 


example02: example02.o  $(COMMONOBJ)
	$(CXX) $(LDFLAGS) -o $@ $@.o $(COMMONOBJ) $(LIBRARIES)


make:
	../scripts/mkcxx.pl '-f' '-i' '-I../Framework' '-l' '-L../Framework -lEventMixer -lz' '-r'

clean:
	rm -vf $(COMMONOBJ) $(PROGOBJ)

realclean: clean
	rm -vf  example02 

.cc.o:         $<
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@
.cpp.o:         $<
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@
.C.o:         $<
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@
.f.o:         $<
	$(FC) $(FFLAGS) -c $< -o $@
.f90.o:         $<
	$(FC) $(FFLAGS) -c $< -o $@


depend:
	makedepend  $(LCLINCLUDE) -Y --   -- $(COMMONSRC) $(PROGSRC)
# DO NOT DELETE

puppiContainer.o: puppiContainer.hh
example02.o: puppiContainer.hh
