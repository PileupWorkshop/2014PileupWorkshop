2014PileupWorkshop
==================

This repository, currently under development, contains shared
software for the [CERN Pileup
Workshop](https://indico.cern.ch/event/306155/) in May 2014.


Useful reading
--------------

In a while we'll be adding a list of useful papers here.

Software
========

The current directory structure is as follows:

- EventMixer/ provides the basis of a library for reading events, mixing hard
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

If fastjet-config is not in your path, place a .fastjet file containing the path
to the top-level fastjet installation directory (i.e. without the /bin/ part) in each of the
directories where you'll run ./mkmk. Then:

In the EventMixer/ directory run
  
    ./mkmk
    make

In the example/ directory run

    ./mkmk
    make
    ./test -hard ../sample-events/dummy.pu14



