
# Makefile generated automatically by ../scripts/mkcxx.pl '-a' 'EventMixer' '-f' '-g' 'c++'
# run 'make make' to update it if you add new files

CXX = g++  # for macs - otherwise get c++ = clang
CXXFLAGS = -Wall -g -O2

# also arrange for fortran support
FC = gfortran
FFLAGS = -Wall -O2
CXX = c++

FJCONFIG = fastjet-config
INCLUDE += `$(FJCONFIG) --cxxflags`
LIBRARIES  += `$(FJCONFIG) --libs --plugins`
INCLUDE += $(LCLINCLUDE)

COMMONSRC = CmdLine.cc EventMixer.cc EventSource.cc ParticleIDMethods.cc ParticleName.cc PU14.cc Version.cc zfstream.cc
F77SRC = 
COMMONOBJ = CmdLine.o EventMixer.o EventSource.o ParticleIDMethods.o ParticleName.o PU14.o Version.o zfstream.o

PROGSRC = 
PROGOBJ = 

INCLUDE += 
LIBRARIES += 


all:  libEventMixer.a


libEventMixer.a: $(COMMONOBJ)
	ar cru libEventMixer.a $(COMMONOBJ)
	ranlib libEventMixer.a


make:
	../scripts/mkcxx.pl '-a' 'EventMixer' '-f' '-g' 'c++'

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
EventMixer.o: EventMixer.hh CmdLine.hh EventSource.hh
EventSource.o: EventSource.hh FastIStringStream.hh PU14.hh
EventSource.o: HepPID/ParticleIDMethods.hh zfstream.h
ParticleIDMethods.o: HepPID/ParticleIDMethods.hh HepPID/ParticleName.hh
ParticleName.o: HepPID/ParticleName.hh HepPID/ParticleIDMethods.hh
ParticleName.o: HepPID/Version.hh
PU14.o: PU14.hh HepPID/ParticleIDMethods.hh
Version.o: HepPID/Version.hh
zfstream.o: zfstream.h
