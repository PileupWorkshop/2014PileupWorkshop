Work program for comparing various subtraction methods in identical setups
==========================================================================

Ideally, each method should be studied by at least two people/groups

Code and results should be uploaded to a subdirectory of the directory
"Comparisons/" on github


Signal samples
--------------
- dijets (pt > 20, pt > 100, pt > 500), with UE off and massless particles
- both full and CHS

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
  antikt R=0.4
  particles: |y| < 4
  jets: select 2 hardest > 20, 100, 500 in hard event, and then |y| < 2.5, 
match to full with deltaR = 0.3 criterion

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
     # npu      <DeltaO>     sigma_DeltaO    corr.coeff.     #_of_jets>20_GeV
       xxx       xxx +- yyy  xxx +- yyy      xxx             xxx +- yyy


To be considered at a later stage
---------------------------------
boosted W, pt > 500R=1 (trimmed mass)
angularity alpha=1
tau_32, beta=2, 1-pass-kt-axis
tau_21, beta=1, 1-pass-kt-axis, m > 40
