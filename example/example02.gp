# gnuplot file

set term postscript enhanced 

filename='example02.ps'
set output filename

set xlabel 'rapidity'
set ylabel 'average offset (<{/Symbol D} p_t>) [GeV]'
plot '<../scripts/mergeidx.pl -f output.dat offset_v_rapidity' u 2:4 w histeps

set ylabel 'dispersion (<{/Symbol s}_{{/Symbol D}p_t}>) [GeV]'
plot '<../scripts/mergeidx.pl -f output.dat offset_v_rapidity' u 2:5 w histeps

# close the file
set output
print "File ".filename." produced with the plots"
