reset 

cd  'C:\Users\gil-c\Documents\Heinz_Research\Bjorken_RTA_2\Extra\rta_resum-master\VS_Version\VS_Version\results'
 
set term wxt 0
plot 'pibar_exact.dat' ps 0.2

set term pdf
set output 'pibar_exact.pdf'
plot 'pibar_exact.dat' ps 0.2 

set term wxt 1
plot 'pibar_derivative_exact.dat' ps 0.2

set term pdf
set output 'pibar_derivative_exact.pdf'
plot 'pibar_derivative_exact.dat' ps 0.2

set term wxt 2
plot 'pibar_2nd_derivative_exact.dat' ps 0.2

set term pdf
set output 'pibar_2nd_derivative_exact.pdf'
plot 'pibar_2nd_derivative_exact.dat' ps 0.2

set term wxt 3
plot 'pibar_3rd_derivative_exact.dat' ps 0.2

set term pdf
set output 'pibar_3rd_derivative_exact.pdf'
plot 'pibar_3rd_derivative_exact.dat' ps 0.2

set term wxt 4
plot 'pibar_4th_derivative_exact.dat' ps 0.2

set term pdf
set output 'pibar_4th_derivative_exact.pdf'
plot 'pibar_4th_derivative_exact.dat' ps 0.2

set term wxt 5 
plot 'tau_r_exact.dat' ps 0.2

set term pdf
set output 'tau_r_exact.pdf'
plot 'tau_r_exact.dat' ps 0.2

set term wxt 6
plot 'z_exact.dat' ps 0.2

set term pdf 
set output 'z_exact.pdf'
plot 'z_exact.dat' ps 0.2

set term wxt 7
set logscale
plot 'T_exact.dat'ps 0.2, 'T_aniso.dat' ps 0.2, 'T_ideal.dat' ps 0.2, 

set term pdf
set output 'T_evol.pdf'
set logscale
plot 'T_exact.dat' ps 0.2, 'T_aniso.dat' ps 0.2, 'T_ideal.dat' ps 0.2,  

set term wxt 8
set logscale 
plot 'e_exact.dat' ps 0.2, 'e_aniso.dat' ps 0.2, 'e_ideal.dat' ps 0.2, 'e_vs_tau_200.dat' u 2:3 ps 0.2

set term pdf
set output 'e_evol.pdf'
plot 'e_exact.dat' ps 0.2, 'e_aniso.dat' ps 0.2, 'e_ideal.dat' ps 0.2, 'e_vs_tau_200.dat' u 2:3 ps 0.2

set term wxt 9
plot x
reset