2014PileupWorkshop
==================

This repository, currently under development, contains shared
software for the [CERN Pileup
Workshop](https://indico.cern.ch/event/306155/) in May 2014.


Useful reading: see the file [READING.md](READING.md) 
---------------

Project suggestions
===================

- are we better off with massless or massive jet inputs?

- comparing new methods on identical samples

- comparing event-wide v. local PU suppression methods

- your own...

Software
========

The current directory structure is as follows:

- Framework/ provides the basis of a library for reading events, mixing hard
  and pileup events, etc.

- example/ will contain a handful of examples

- scripts/ contains some helper scripts (e.g. for automated Makefile
  generation). 

- sample-events/ contains small event files for basic testing

- [/afs/cern.ch/user/p/puws2014/public/events](file:///afs/cern.ch/user/p/puws2014/public/events) contains large event
  files, also accessible through the web at
  [http://cern.ch/puws2014/events/](http://cern.ch/puws2014/events/).  File names should be fairly self-explanatory. 


Getting started
---------------

Download the software from github if you haven't already done so:

    git clone https://github.com/PileupWorkshop/2014PileupWorkshop.git

If fastjet-config is not in your path, place a .fastjet file containing the path
to the top-level fastjet installation directory (i.e. without the /bin/ part) in each of the
directories where you'll run ./mkmk (see below). Then:

In the Framework/ directory run
  
    ./mkmk
    make

In the example/ directory run

    ./mkmk
    make

First try out a very simple few-line test program that illustrates
some basic features of reading events

    ./test -hard ../sample-events/dummy.pu14

Next try a program that mixes hard and pileup events and prints out a
couple of characteristics

    ./example01 -hard ../sample-events/lhc14-pythia8-4C-dijet50-nev20.pu14.gz \
                -pileup ../sample-events/lhc14-pythia8-4C-minbias-nev100.pu14.gz \
                -npu 20 -nev 2

This adds a fixed number of pileup events (20) to each hard event
(sorry, no Poisson fluctuations yet!).

The last example that's currently in place, takes hard and pileup
events, runs area-based pileup subtraction, and generates some
histograms statistics on the quality of the subtraction. 

    ./example02 -hard ../sample-events/lhc14-pythia8-4C-dijet50-nev20.pu14.gz \
                -pileup ../sample-events/lhc14-pythia8-4C-minbias-nev100.pu14.gz \
                -npu 5 -nev 20 > output.dat

and look inside output.dat to see what's there. To actually get sensible
answers, you'll need more than 20 events ("-nev 20"). For that, run with the
big event samples (see above). You can plot results with

    gnuplot example02.gp

which will produce a file example02.ps.

Look inside the [example02.cc](example/example02.cc)
program to see some of the options. The matching of full jets and hard
jets is performed with a simple geometrical method for now.


FastJet-contrib:
----------------

We have provided a branch with development versions of a few tools 

  svn co https://fastjet.hepforge.org/svn/contrib/branches/1.012-alpha-PUWS14.1

or, create a file in your contrib directory, called contribs.local,
with the following contents

  # addition of a SafeAreaSubtractor (including "rho_m") and SafeNpCSubtractor
  GenericSubtractor                tags/2.0.0-alpha-PUWS14.1
  
  # SoftKiller
  SoftKiller                       tags/1.0.0-alpha-PUWS14.1
  
  # ModifiedMassDropTagged and SoftDrop
  RecursiveTools                   tags/1.0-alpha-PUWS14.1
  
  # other relevant contribs included in stable releases
  JetCleanser                      tags/1.0.0
  ConstituentSubtractor            tags/1.0.0


Installation (see README in fastjet-contrib for more details):  

  ./scripts/update-contribs.sh
  ./configure    # assumes fastjet-config is in your path

  make
  make check  # optional
  make install

Then add the relevant -lContribXXX -lContribYYY to your Makefile.

Note that the "make install" step would install fastjet-contrib in
your FastJet directory which can interfere with prior installations of
fastjet-contrib.

Shared-library alternative:
^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you're using shared FastJet libraries, you can try out 

  make fragile-shared
  make check  # optional
  make fragile-shared-install

Then add the relevant -lfastjetcontribfragile to your Makefile
