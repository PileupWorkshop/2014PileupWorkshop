reset 

set term postscript enhanced color font "Helvetica, 13" dl 1 size 17cm,12cm
set output 'offset-v-dispersion-single.eps'
set colors classic

set yrange [0:*]
set grid

# set key left top
set key top right

fn(pt,npu)=sprintf('res/sub-dijets%s-npu%d-novor.res',pt,npu)
all_v_npu(pt)='< cat '.fn(pt,30).' '.fn(pt,60).' '.fn(pt,100).' '.fn(pt,140)
all_v_pt(npu)='< cat '.fn('020',npu).' '.fn('050',npu).' '.fn('100',npu).' '.fn('500',npu)

unset key 
set bmargin at screen 0.13
set tmargin at screen 0.92
set xtics 1

label(i,tag)='{/=13 '.tag.'}'

# line styles
set style line 1 lc 4             dt (1,3)     lw 4 pt 2  ps 1   # Cleansing
set style line 2 lc 3             dt (1,3,5,3) lw 4 pt 13 ps 1.2 # Filter
set style line 3 lc 7             dt (5,3)     lw 4 pt 5  ps 1   # SK
set style line 4 lc 7             dt 1         lw 4 pt 4  ps 0.9 # SK+0
set style line 5 lc rgb "#00aa00" dt (3,3)     lw 4 pt 9  ps 1.2 # PUPPI
set style line 6 lc 1             dt 1         lw 4 pt 6  ps 1   # area-median

do for [obs in "pt m"]{
    set output 'offset-v-dispersion-single-'.obs.'.eps'
    #------------------------------------------------------------------------
    # NPU dependence, different pt in each panels
    set multiplot

    if (obs eq "pt"){
        set ylabel '{/Symbol s}_{{/Symbol D}p_t} [GeV]' offset 1.5
        set xlabel '<{/Symbol D}p_t> [GeV]'
        set xrange [-2.5:1.5]
        set yrange [0:12]
        set arrow 1 nohead from 0.0,graph 0.0 to 0.0, graph 1.0 lt 1 dt 1 lc 7 lw 1 back
    } else {
        set ylabel '{/Symbol s}_{{/Symbol D}m} [GeV]' offset 0.0
        set xlabel '<{/Symbol D}m> [GeV]'
        set xrange [-1.2:2.8]
        set yrange [0:7]
        set arrow 1 nohead from 0.0,graph 0.0 to 0.0, graph 0.67 lt 1 dt 1 lc 7 lw 1 back
    }
    set format y "%g"

    unset key

    #set style rectangle fs solid 1.0 fc "#ffffff" back noclip noborder
    set style line 9 lc "#ff0000"

    pts = "020 050 100 500"
    do for [ipt=1:4]{
        pt=word(pts,ipt)
        set lmargin at screen 0.22*ipt-0.12
        set rmargin at screen 0.22*ipt+0.10
        set label 1 'p_t>'.sprintf("%g",pt+0.0).' GeV' center at graph 0.5,1.05
        plot all_v_npu(pt)." | grep ".obs."_clns"     u ($3):($6) w lp ls 1 t label(ipt,'linear cleansing'),\
             all_v_npu(pt)." | grep ".obs."_filt0203" u ($3):($6) w lp ls 2 t label(ipt,'Filter'),\
             all_v_npu(pt)." | grep ".obs."_sk50"     u ($3):($6) w lp ls 3 t label(ipt,'soft killer (SK)'),\
             all_v_npu(pt)." | grep ".obs."_sk45z02"  u ($3):($6) w lp ls 4 t label(ipt,'SK+Zeroing'),\
             all_v_npu(pt)." | grep ".obs."_puppi"    u ($3):($6) w lp ls 5 t label(ipt,'PUPPI'),\
             all_v_npu(pt)." | grep ".obs."_area"     u ($3):($6) w lp ls 6 t label(ipt,'area-median')
        
        set format y ""
        unset ylabel
        if (ipt==1){
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

    #------------------------------------------------------------------------
    # pt dependence, different NPU in each panels
    set output 'offset-v-dispersion-single-'.obs.'-vpt.eps'
    set multiplot
    if (obs eq "pt"){
        set ylabel '{/Symbol s}_{{/Symbol D}p_t} [GeV]' offset 1.5
        set xlabel '<{/Symbol D}p_t> [GeV]'
        set xrange [-2.5:1.5]
        set yrange [0:13]
    } else {
        set ylabel '{/Symbol s}_{{/Symbol D}m} [GeV]' offset 0.0
        set xlabel '<{/Symbol D}m> [GeV]'
        set xrange [-1.2:2.8]
        set yrange [0:7]
    }
    set format y "%g"

    npus="30 60 100 140"
    do for [inpu=1:4]{
        npu=word(npus,inpu)+0.0
        set lmargin at screen 0.22*inpu-0.12
        set rmargin at screen 0.22*inpu+0.10
        set label 1 'N_{PU}='.sprintf("%g",npu) center at graph 0.5,1.05
        plot all_v_pt(npu)." | grep ".obs."_clns"     u ($3):($6) w lp ls 1 t label(inpu,'linear cleansing'),\
             all_v_pt(npu)." | grep ".obs."_filt0203" u ($3):($6) w lp ls 2 t label(inpu,'Filter'),\
             all_v_pt(npu)." | grep ".obs."_sk50"     u ($3):($6) w lp ls 3 t label(inpu,'soft killer (SK)'),\
             all_v_pt(npu)." | grep ".obs."_sk45z02"  u ($3):($6) w lp ls 4 t label(inpu,'SK+Zeroing'),\
             all_v_pt(npu)." | grep ".obs."_puppi"    u ($3):($6) w lp ls 5 t label(inpu,'PUPPI'),\
             all_v_pt(npu)." | grep ".obs."_area"     u ($3):($6) w lp ls 6 t label(inpu,'area-median')

        set format y ""
        unset ylabel
        if (inpu==1){
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

}

set out
