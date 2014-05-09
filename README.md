2014PileupWorkshop
==================

This repository, currently under development, will contain shared
software for the [CERN Pileup
Workshop](https://indico.cern.ch/event/306155/) in May 2014.


Useful reading
--------------

In a while we'll be adding a list of useful papers here.

Software
========

The current development plan is as follows:

- EventMixer/ will provide a library for reading events, mixing hard
  and pileup events, etc.

- example/ will contain a handful of examples

- scripts/ contains some helper scripts (e.g. for automated Makefile
  generation). 

- sample-events/ contains small event files for basic testing

- /afs/cern.ch/user/p/puws2014/public/events contains large event
  files. Names should be fairly self-explanatory. (We should add http
  access to the files too.)


Getting started
---------------

- in the EventMixer/ directory run
  
    ./mkmk
    make

in the example/ directory run

    ./mkmk
    make
    ./test -hard ../sample-events/dummy.pu14



