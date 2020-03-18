reset
cd 'C:\Users\gil-c\Documents\Heinz_Research\Bjorken_RTA_2\Project1\output'
datafile = "shear_Ingles_alt.dat"
hbarc = 0.1973

f(x) = a * x**(-4./3.)
fit f(x) datafile u 1:2 via a 

set term wxt 0
set ylabel "energy density [GeV/fm^3]"
set xlabel "proper time [fm/c]"
#set logscale y
plot datafile u 1:3 ps 0.2 lc rgb 'red' title 'Strickland', 'e_exact.dat' u 1:2 ps 0.2 lc rgb 'green' title 'McNelis', 'e_vs_tau_200.dat' u 2:3 ps 0.2 lc rgb 'black' title 'Gojko', datafile u 1:2 ps 0.2 lc rgb 'blue' title 'Ingles'

set terminal pdf
set output 'e_evolution.pdf'
set logscale
plot datafile u 1:2 ps 0.2 lc rgb 'blue' title 'Ingles', datafile u 1:3 ps 0.2 lc rgb 'red' title 'Strickland', 'e_exact.dat' u 1:2 ps 0.2 lc rgb 'green' title 'McNelis', 'e_vs_tau_200.dat' u 2:3 ps 0.2 lc rgb 'black' title 'Gojko'

set term wxt 1
set ylabel "longitude pressure [GeV/fm^3]"
unset logscale y
plot datafile u 1:6 ps 0.2 lc rgb 'blue' title 'Ingles', datafile u 1:7 ps 0.2 lc rgb 'red' title 'Strickland'

set term pdf 
set output 'pl_evolution.pdf'
plot datafile u 1:6 ps 0.2 lc rgb 'blue' title 'Ingles', datafile u 1:7 ps 0.2 lc rgb 'red' title 'Strickland'

set term wxt 2 
set ylabel 'tansverse pressure [GeV/fm^3]'
plot datafile u 1:4 ps 0.2 lc rgb 'blue' title 'Ingles', datafile u 1:5 ps 0.2 lc rgb 'red' title 'Strickland'

set terminal pdf
set output 'pt_evolution.pdf'
plot datafile u 1:4 ps 0.2 lc rgb 'blue' title 'Ingles', datafile u 1:5 ps 0.2 lc rgb 'red' title 'Strickland'

set term wxt 3
plot x
