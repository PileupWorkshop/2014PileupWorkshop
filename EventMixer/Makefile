
# Makefile generated automatically by ../scripts/mkcxx.pl '-a' 'EventMixer' '-f'
# run 'make make' to update it if you add new files

CXX = g++  # for macs - otherwise get c++ = clang
CXXFLAGS = -Wall -g -O2

# also arrange for fortran support
FC = gfortran
FFLAGS = -Wall -O2

FJCONFIG = fastjet-config
INCLUDE += `$(FJCONFIG) --cxxflags`
LIBRARIES  += `$(FJCONFIG) --libs --plugins`
INCLUDE += $(LCLINCLUDE)

COMMONSRC = CmdLine.cc EventSource.cc ParticleIDMethods.cc ParticleName.cc Version.cc
F77SRC = 
COMMONOBJ = CmdLine.o EventSource.o ParticleIDMethods.o ParticleName.o Version.o

PROGSRC = 
PROGOBJ = 

INCLUDE += 
LIBRARIES += 


all:  libEventMixer.a


libEventMixer.a: $(COMMONOBJ)
	ar cru libEventMixer.a $(COMMONOBJ)
	ranlib libEventMixer.a


make:
	../scripts/mkcxx.pl '-a' 'EventMixer' '-f'

clean:
	rm -vf $(COMMONOBJ) $(PROGOBJ)

realclean: clean
	rm -vf  libEventMixer.a

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

CmdLine.o: CmdLine.hh
EventSource.o: EventSource.hh FastIStringStream.hh PU14.hh
EventSource.o: HepPID/ParticleIDMethods.hh
ParticleIDMethods.o: HepPID/ParticleIDMethods.hh HepPID/ParticleName.hh
ParticleName.o: HepPID/ParticleName.hh HepPID/ParticleIDMethods.hh
ParticleName.o: HepPID/Version.hh
Version.o: HepPID/Version.hh
