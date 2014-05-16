reset 

set term postscript enhanced color font "Helvetica, 15"
set output 'offset-v-dispersion.eps'


set xrange [-0.1:0.1]
set yrange [0:0.2]
set grid
set ylabel '{/Symbol s}_{{/Symbol D}O}/<O_{hard}>'
set xlabel '<{/Symbol D}O>/<O_{hard}>'

merge='<~/scripts/mergeidx.pl -f output.dat '# npu'

set key left top
#set label ' n_{PU} = 15, 55, 95, 135' at 0.018, 0.19 right
set label 'O = jetpt' at  0.09, 0.19 right 
set label '<p_{t,hard}> = 50 GeV' at 0.09, 0.18 right 

plot '<~/scripts/mergeidx.pl -f output.dat "# npu"'  u ($2/50):($3/50) w lp lw 3 ps 2 t 'full areasub'
