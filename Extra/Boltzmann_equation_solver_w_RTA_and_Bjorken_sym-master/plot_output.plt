reset
cd 'C:\Users\gil-c\Documents\Heinz_Research\Bjorken_RTA_2\Extra\Boltzmann_equation_solver_w_RTA_and_Bjorken_sym-master'

set ylabel "beta [GeV^{-1}]"
set xlabel "proper time [fm/c]"

FILES = system("dir /b *.dat")
#set logscale y
plot for [data in FILES] data u 2:3 ps 0.2 title data

set term pdf 
set output 'e_evolution_gojko.pdf'
#set logscale y
plot for [data in FILES] data u 2:3 ps 0.2 notitle