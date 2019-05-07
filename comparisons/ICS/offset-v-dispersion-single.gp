# This file was copied from comparisons/review/offset-v-dispersion-single.gp, version 1.0.0 and it was further modified


reset 

set term postscript enhanced color font "Helvetica, 13" dl 1 size 17cm,12cm
#set output 'offset-v-dispersion-single.eps'
set colors classic

set yrange [0:*]
set grid

# set key left top
set key top right

unset key 
set bmargin at screen 0.13
set tmargin at screen 0.92
set xtics 1

label(i,tag)='{/=13 '.tag.'}'

# line styles
set style line 1 lc 1             dt 1         lw 3 pt 6  ps 1   # area-median
set style line 2 lc 7             dt 1         lw 3 pt 4  ps 1 # SK+0
set style line 3 lc rgb "#00aa00" dt 1     lw 3 pt 9  ps 1.4 # PUPPI
set style line 4 lc 3             dt 1 lw 3 pt 11 ps 1.4 # Jet-by-jet CS
set style line 5 lc 4             dt 1     lw 3 pt 68  ps 1.2   # Event-wide CS
set style line 6 lc 5             dt 1     lw 3 pt 2  ps 1.3   # ICS

#do for [radius in "0.4 1.0"]{
do for [radius in "0.4"]{
#do for [massOption in "massless withMass"]{
do for [massOption in "massless"]{

#dir="res_previous3"
dir="res"

fn(pt,npu)=sprintf('%s/sub-dijetsel%s-noUE-nevsel-npu%d-radius'.radius.'-'.massOption.'.res',dir,pt,npu)
fn2(sample,npu)=sprintf('%s/sub-%s-npu%d-radius'.radius.'-'.massOption.'.res',dir,sample,npu)
all_v_npu(pt)='< cat '.fn(pt,30).' '.fn(pt,60).' '.fn(pt,100).' '.fn(pt,140)
all_v_npu2(sample)='< cat '.fn2(sample,30).' '.fn2(sample,60).' '.fn2(sample,100).' '.fn2(sample,140)
all_v_pt(npu)='< cat '.fn('20',npu).' '.fn('50',npu).' '.fn('100',npu).' '.fn('500',npu)

do for [obs in "pt m width"]{
    #------------------------------------------------------------------------
    set output dir.'/offset-v-dispersion-single-'.obs.'-radius'.radius.'-'.massOption.'.eps'
    # NPU dependence, different pt in each panels
    set multiplot


    if (obs eq "pt"){
        set ylabel '{/Symbol s}_{{/Symbol D}p_T} [GeV]' offset 1.5
        set xlabel '<{/Symbol D}p_T> [GeV]'
        set xrange [-1.65:1.2]
        set yrange [0:11.6]
	if (radius eq "1.0"){
	     set xrange [-5:3]
             set yrange [0:24]
	}
        set arrow 1 nohead from 0.0,graph 0.0 to 0.0, graph 1.0 lt 1 dt 1 lc 7 lw 1 back
    }
    if (obs eq "width"){
        set ylabel '{/Symbol s}_{{/Symbol D}width}' offset 1.5
        set xlabel '<{/Symbol D}width>'
        set xrange [-0.03:0.03]
        set yrange [0:0.4]
	if (radius eq "1.0"){
	     set xrange [-0.06:0.06]
             set yrange [0:0.8]
	}
        set arrow 1 nohead from 0.0,graph 0.0 to 0.0, graph 1.0 lt 1 dt 1 lc 7 lw 1 back
    } 
    if (obs eq "m"){
        set ylabel '{/Symbol s}_{{/Symbol D}m} [GeV]' offset 0.0
        set xlabel '<{/Symbol D}m> [GeV]'
        set xrange [-1.1:1.9]
        set yrange [0:6.8]
	if (radius eq "1.0"){
	     set xrange [-2.4:5]
             set yrange [0:14]
	}
        set arrow 1 nohead from 0.0,graph 0.0 to 0.0, graph 1.0 lt 1 dt 1 lc 7 lw 1 back
    }
    set format y "%g"

    unset key

    #set style rectangle fs solid 1.0 fc "#ffffff" back noclip noborder
    set style line 9 lc "#ff0000"

    pts = "20 50 100 500"

    do for [ipt=1:4]{
        pt=word(pts,ipt)
        set lmargin at screen 0.22*ipt-0.12
        set rmargin at screen 0.22*ipt+0.10
        set label 1 'p_T>'.sprintf("%g",pt+0.0).' GeV' center at graph 0.5,1.05


plot             all_v_npu(pt)." | grep ".obs."_area"     u ($3):($6) w lp ls 1 t label(ipt,'Area Subtraction'),\
             all_v_npu(pt)." | grep ".obs."_constit"  u ($3):($6) w lp ls 4 t label(ipt,'Jet-by-jet CS'),\
             all_v_npu(pt)." | grep ".obs."_CSevent"  u ($3):($6) w lp ls 5 t label(ipt,'Event-wide CS'),\
             all_v_npu(pt)." | grep ".obs."_ICS4"  u ($3):($6) w lp ls 6 t label(ipt,'ICS')
        
        set format y ""
        unset ylabel
        if (ipt==1){
            set key at graph 0.42,0.975 width -3.2 spacing 1.05 box lw 1.5 samplen 3.5
            set object 1 rectangle from graph -0.89,0.773 to graph 0.1, 0.975 back  noclip fs solid border rgb "#ffffff" fc "#ffffff"
            set object 2 rectangle from graph -0.01,0.776 to graph 0.01,0.970 front noclip fs solid border rgb "#ffffff" fc "#ffffff"
        } else {
            unset key
            unset object 1
            unset object 2
        }
    }
    unset multiplot
    set yrange [0:*]


    #------------------------------------------------------------------------
    set output dir.'/offset-v-dispersion-single-'.obs.'-radius'.radius.'-'.massOption.'2.eps'
    # NPU dependence, different pt in each panels
    set multiplot
        
    if (obs eq "pt"){
        set ylabel '{/Symbol s}_{{/Symbol D}p_T} [GeV]' offset 0.3
        set xlabel '<{/Symbol D}p_T> [GeV]'
        set xrange [-2.999:1.5]
        set yrange [0:8.5]
	if (radius eq "1.0"){
	     set xrange [-5:3]
             set yrange [0:24]
	}
        set arrow 1 nohead from 0.0,graph 0.0 to 0.0, graph 1.0 lt 1 dt 1 lc 7 lw 1 back
    }
    if (obs eq "width"){
        set ylabel '{/Symbol s}_{{/Symbol D}width}' offset 1.5
        set xlabel '<{/Symbol D}width>'
        set xrange [-0.005:0.01]
        set yrange [0:0.2]
	if (radius eq "1.0"){
	     set xrange [-0.06:0.06]
             set yrange [0:0.8]
	}
        set arrow 1 nohead from 0.0,graph 0.0 to 0.0, graph 1.0 lt 1 dt 1 lc 7 lw 1 back
    } 
    if (obs eq "m"){
        set ylabel '{/Symbol s}_{{/Symbol D}m} [GeV]' offset 0.5
        set xlabel '<{/Symbol D}m> [GeV]'
        set xrange [-0.4:1.05]
        set yrange [0:4.25]
	if (radius eq "1.0"){
	     set xrange [-2.4:5]
             set yrange [0:14]
	}
        set arrow 1 nohead from 0.0,graph 0.0 to 0.0, graph 1 lt 1 dt 1 lc 7 lw 1 back
    }
    set format y "%g"

    unset key

    #set style rectangle fs solid 1.0 fc "#ffffff" back noclip noborder
    set style line 9 lc "#ff0000"

    pts = "20 50 100 500"
    do for [ipt=1:4]{
        pt=word(pts,ipt)
        set lmargin at screen 0.22*ipt-0.12
        set rmargin at screen 0.22*ipt+0.10
        set label 1 'p_T>'.sprintf("%g",pt+0.0).' GeV' center at graph 0.5,1.05

#plot         all_v_npu(pt)." | grep ".obs."_sk45z02"     u ($3):($6) w lp ls 2 t label(ipt,'SK+Zeroing'),\

#plot         all_v_npu(pt)." | grep ".obs."_sk50"     u ($3):($6) w lp ls 2 t label(ipt,'SK'),\
#	     all_v_npu(pt)." | grep ".obs."_puppi"    u ($3):($6) w lp ls 3 t label(ipt,'PUPPI'),\
#             all_v_npu(pt)." | grep ".obs."_CSSK1"  u ($3):($6) w lp ls 5 t label(ipt,'CS+SK'),\
#             all_v_npu(pt)." | grep ".obs."_ICS4"  u ($3):($6) w lp ls 6 t label(ipt,'ICS')

plot         all_v_npu(pt)." | grep ".obs."_ICS5"     u ($3):($6) w lp ls 2 t label(ipt,'ICS5'),\
	     all_v_npu(pt)." | grep ".obs."_ICS6"    u ($3):($6) w lp ls 3 t label(ipt,'ICS6'),\
	     all_v_npu(pt)." | grep ".obs."_CSSK1"    u ($3):($6) w lp ls 4 t label(ipt,'CSSK1'),\
	     all_v_npu(pt)." | grep ".obs."_CSSK2"    u ($3):($6) w lp ls 1 t label(ipt,'CSSK2'),\
           all_v_npu(pt)." | grep ".obs."_CSSK3"  u ($3):($6) w lp ls 5 t label(ipt,'CSSK3'),\
             all_v_npu(pt)." | grep ".obs."_ICS4"  u ($3):($6) w lp ls 6 t label(ipt,'ICS4')

        set format y ""
        unset ylabel
	if (obs eq "width"){
        if (ipt==3){
            set key at graph 0.39,0.975 width -4. spacing 1.05 box lw 1.5 samplen 3.5
            set object 1 rectangle from graph -0.89,0.773 to graph 0.1, 0.975 back  noclip fs solid border rgb "#ffffff" fc "#ffffff"
            set object 2 rectangle from graph -0.01,0.773 to graph 0.01,0.973 front noclip fs solid border rgb "#ffffff" fc "#ffffff"
        } else {
            unset key
            unset object 1
            unset object 2
        }
}
else{
        if (ipt==1){
            set key at graph 0.01,0.975 width 0.9 spacing 1.05 box lw 1.5 samplen 3.5
            set object 1 rectangle from graph -0.89,0.773 to graph 0.1, 0.975 back  noclip fs solid border rgb "#ffffff" fc "#ffffff"
#            set object 2 rectangle from graph -0.01,0.775 to graph 0.01,0.971 front noclip fs solid border rgb "#ffffff" fc "#ffffff"
        } else {
            unset key
            unset object 1
            unset object 2
        }
}

    }
    unset multiplot
    set yrange [0:*]


    if (0) { # commented-out block for other samples
    #------------------------------------------------------------------------
    set output dir.'/offset-v-dispersion-single-'.obs.'-radius'.radius.'-'.massOption.'-allSamples.eps'
    # NPU dependence, all samples
    set multiplot


    if (obs eq "pt"){
        set ylabel '{/Symbol s}_{{/Symbol D}p_T} [GeV]' offset 1.5
	set xlabel '<{/Symbol D}p_T> [GeV]'
        set xrange [-2.5:1.5]
        set yrange [0:12]
        if (radius eq "1.0"){
             set xrange [-5:3]
             set yrange [0:24]
        }
        set arrow 1 nohead from 0.0,graph 0.0 to 0.0, graph 1.0 lt 1 dt 1 lc 7 lw 1 back
    }
    if (obs eq "width"){
    set ylabel '{/Symbol s}_{{/Symbol D}width}' offset 1.5
        set xlabel '<{/Symbol D}width>'
        set xrange [-0.03:0.03]
        set yrange [0:0.4]
        if (radius eq "1.0"){
             set xrange [-0.06:0.06]
             set yrange [0:0.8]
	     }
        set arrow 1 nohead from 0.0,graph 0.0 to 0.0, graph 1.0 lt 1 dt 1 lc 7 lw 1 back
    }
    if (obs eq "m"){
        set ylabel '{/Symbol s}_{{/Symbol D}m} [GeV]' offset 0.0
        set xlabel '<{/Symbol D}m> [GeV]'
        set xrange [-1.2:2.8]
        set yrange [0:7]
	if (radius eq "1.0"){
             set xrange [-2.4:5]
             set yrange [0:14]
        }
        set arrow 1 nohead from 0.0,graph 0.0 to 0.0, graph 0.67 lt 1 dt 1 lc 7 lw 1 back
    }
    set format y "%g"

    unset key
    #set style rectangle fs solid 1.0 fc "#ffffff" back noclip noborder
    set style line 9 lc "#ff0000"



    samples = "WW500-noUE-nev Zprime500-noUE-nev dijetsel100-noUE-nevsel dijetsel500-noUE-nevsel"
    sampleNames = "WW Zprime dijets100 dijets500"
    do for [isample=1:4]{
        sample=word(samples,isample)
        sampleName=word(sampleNames,isample)
        set lmargin at screen 0.22*isample-0.12
        set rmargin at screen 0.22*isample+0.10
        set label 1 sprintf("%s",sampleName) center at graph 0.5,1.05
#        plot \
#all_v_npu2(sample)." | grep ".obs."_ICS1"     u ($3):($6) w lp ls 1 t label(isample,'0.2 0.1 false'),\
#             all_v_npu2(sample)." | grep ".obs."_ICS2" u ($3):($6) w lp ls 2 t label(isample,'0.15 0.25 false'),\
#             all_v_npu2(sample)." | grep ".obs."_ICS3"     u ($3):($6) w lp ls 3 t label(isample,'0.15 0.15 false'),\
#             all_v_npu2(sample)." | grep ".obs."_ICS4"  u ($3):($6) w lp ls 4 t label(isample,'0.15 0.15 true'),\
#             all_v_npu2(sample)." | grep ".obs."_ICS5"    u ($3):($6) w lp ls 5 t label(isample,'0.1 0.15 false'),\
#             all_v_npu2(sample)." | grep ".obs."_ICS6"     u ($3):($6) w lp ls 6 t label(isample,'0.1 0.2 false')
plot             all_v_npu2(sample)." | grep ".obs."_area"     u ($3):($6) w lp ls 1 t label(isample,'Area Subtraction'),\
             all_v_npu2(sample)." | grep ".obs."_sk45z02"     u ($3):($6) w lp ls 2 t label(isample,'SK + Zeroing'),\
	     all_v_npu2(sample)." | grep ".obs."_puppi"    u ($3):($6) w lp ls 3 t label(isample,'PUPPI'),\
             all_v_npu2(sample)." | grep ".obs."_constit"  u ($3):($6) w lp ls 4 t label(isample,'Jet-by-jet CS'),\
             all_v_npu2(sample)." | grep ".obs."_CSevent"  u ($3):($6) w lp ls 5 t label(isample,'Event-wide CS'),\
             all_v_npu2(sample)." | grep ".obs."_ICS3"  u ($3):($6) w lp ls 6 t label(isample,'ICS')
        
        set format y ""
        unset ylabel
        if (isample==1){
            set key at graph 0.42,0.975 width -4 spacing 1.05 box lw 1.5 samplen 3.5
            set object 1 rectangle from graph -0.89,0.672 to graph 0.1, 0.975 back  noclip fs solid border rgb "#ffffff" fc "#ffffff"
            set object 2 rectangle from graph -0.01,0.674 to graph 0.01,0.973 front noclip fs solid border rgb "#ffffff" fc "#ffffff"
        } else {
            unset key
            unset object 1
            unset object 2
        }
    }
    unset multiplot
    set yrange [0:*]


}   # ending brace for comment
}
}
}

set out
