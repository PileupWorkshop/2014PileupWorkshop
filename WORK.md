Work program for comparing various subtraction methods in identical setups
==========================================================================

Ideally, each method should be studied by at least two people/groups

Code and results should be uploaded to a subdirectory of the directory
"Comparisons/" on github


Signal samples
--------------
- dijets (pt > 20, pt > 100, pt > 500), with UE off and massless particles
- both full and CHS

The samples are in the usual location
([/afs/cern.ch/user/p/puws2014/public/events](file:///afs/cern.ch/user/p/puws2014/public/events)
or
[http://cern.ch/puws2014/events/](http://cern.ch/puws2014/events/)). 

    lhc14-pythia8-4C-dijetsel20-noUE-nevsel1e5.pu14.gz
    lhc14-pythia8-4C-dijetsel100-noUE-nevsel1e5.pu14.gz
    lhc14-pythia8-4C-dijetsel500-noUE-nevsel1e5.pu14.gz

Each event in those files has at least one jet with |y|<2.5 and above the
pt [GeV] indicated after the "dijetsel" tag. They are intended to be
used respectively for the analyses with 20, 50 and 100 GeV pt cuts. 

The files have been produced with a generation cut at 80% of jet pt
selection cut. For validation purposes, there are files (on afs)
labelled .res that indicate the average number of jets with pt>20 GeV,
|y|<2.5 in each sample.


List of subtractors
-------------------
- Safearea
- SafeNpC
- ConstituentSubtractor
- SoftKiller
- Puppi
- Cleansing
- corrJVF
- ......

Each one including trimming (Rsub=0.3, fcut=0.05) where possible

Observables
-----------
- pt
- mass
- angularity/width/girth with alpha = 1

Pileup levels: 30, 60, 100, 140
-------------------------------

particles and jet selection
---------------------------
  particles: |y| < 4
  jet definition: antikt R=0.4

  jets: take the two hardest jets, then apply a selection of |y|<2.5,
  and pt > 20 (or 100, or 500), and study the impact of pileup on any
  jets that pass that selection. Pileup jets are matched to the hard
  jets with a deltaR = 0.3 criterion.

  A separate study counts the jets above 20 GeV with |y| < 2.5. That
  count is only in events that have at one jet from the selection
  described in the preceding paragraph.

  Note that in example03.cc prior to revision 107 (git rev-list --count HEAD)
  or hash d1b6590c2f2758d765c3... the number of jets was counted for
  all events.

How to compare (quality measures)
---------------------------------
- offset v. dispersion (use trimmed jet as ref when using trimming). 
  There' a
template offset-v-dispersion.gp gnuplot macro in the example/ directory for quick
plotting
- number of jets above 20 GeV as a function of npu

File format for results
-----------------------
Use the following file format for each subtractor/observable/sample (write out
the "+-" to output):

     # comments  
     # npu    jet_ptmin  <DeltaO>     sigma_DeltaO    corr.coeff.     #_of_jets>20_GeV   matching_efficiency    <O_hard>           obs & method names 
       xxx    xxx        xxx +- yyy   xxx +- yyy      xxx             xxx +- yyy         xxx +- yyy             xxx +- yyy       #  e.g. pt_areasub

Template code
-------------

Look at [example/example03.cc](example/example03.cc) to see code that
has the selection, matching and output as illustrated above. Run it
for example with 

    ./example03 -hard ../sample-events/lhc14-pythia8-4C-dijet50-nev20.pu14.gz \
                -pileup ../sample-events/lhc14-pythia8-4C-minbias-nev100.pu14.gz \
                -massless -npu 5 -nev 20 -jet-ptmin 20 > output.dat

Adapt the -jet-ptmin option depending on what you plan to
study. There's also a -chs option for running with CHS type events.

To be considered at a later stage
---------------------------------
boosted W, pt > 500R=1 (trimmed mass)
angularity alpha=1
tau_32, beta=2, 1-pass-kt-axis
tau_21, beta=1, 1-pass-kt-axis, m > 40
