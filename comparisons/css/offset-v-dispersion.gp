reset 

set term postscript enhanced color font "Helvetica, 15"
#set output 'offset-v-dispersion-full.eps'
set output 'offset-v-dispersion-chs.eps'


set xrange [-5:5]
set yrange [0:5]
set grid
#set ylabel '{/Symbol s}_{{/Symbol D}O}/<O_{hard}>'
#set xlabel '<{/Symbol D}O>/<O_{hard}>'
set ylabel '{/Symbol s}_{{/Symbol D}O} (GeV)'
set xlabel '<{/Symbol D}O> (GeV)'

#merge='<~/scripts/mergeidx.pl -f output.dat '# npu'

set key left top
#set label ' n_{PU} = 30, 60' at 0.018, 0.19 right
#set label 'O = jet mass' at  0.09, 0.19 right 
#set label '<p_{t,hard}> = 50 GeV' at 0.09, 0.18 right 

set title 'Full event, jet mass'

# full
# set title 'Full event, jet mass'
# plot "<cat dijets20-full-npu*.out | grep m_areasub"           u ($3):($6) w lp lc 1 lw 3 ps 2 t 'areasub' ,\
#      "<cat dijets20-full-npu*.out | grep m_constitsub"        u ($3):($6) w lp lc 2 lw 3 ps 2 t 'constitsub' ,\
#      "<cat dijets20-full-npu*.out | grep m_safeareasub"       u ($3):($6) w lp lc 3 lw 3 ps 2 t 'safeareasub'  ,\
#      "<cat dijets20-full-npu*.out | grep m_soft_killer"       u ($3):($6) w lp lc 4 lw 3 ps 2 t 'soft killer'

# chs
set title 'CHS event, jet mass'
plot "<cat dijets20-chs-npu*.out | grep m_areasub"           u ($3):($6) w lp lc 1 lw 3 ps 2 t 'areasub' ,\
     "<cat dijets20-chs-npu*.out | grep m_constitsub"        u ($3):($6) w lp lc 2 lw 3 ps 2 t 'constitsub' ,\
     "<cat dijets20-chs-npu*.out | grep m_safeareasub"       u ($3):($6) w lp lc 3 lw 3 ps 2 t 'safeareasub'  ,\
     "<cat dijets20-chs-npu*.out | grep m_soft_killer"       u ($3):($6) w lp lc 4 lw 3 ps 2 t 'soft killer' ,\
     "<cat dijets20-chs-npu*.out | grep m_puppi"             u ($3):($6) w lp lc 5 lw 3 ps 2 t 'puppi'  ,\
     "<cat dijets20-chs-npu*.out | grep m_safenpcsub"        u ($3):($6) w lp lc 6 lw 3 ps 2 t 'safenpcsub' ,\
     "<cat dijets20-chs-npu*.out | grep m_linear_cleansing"  u ($3):($6) w lp lc 7 lw 3 ps 2 t 'linear cleansing'
