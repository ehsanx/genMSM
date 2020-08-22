clear
clear matrix
clear mata
set more off
forvalues f=1/1000 {
do "https://raw.githubusercontent.com/ehsanx/genMSM/master/genmsm.do"
mata: outputx = msm(newx = `f', subjects=2500, tpoints=10)
svmat double outputx, name(outputx)
renvars outputx1-outputx19 \ id tpoint tpoint2 T0 IT0 chk y ym a am1 l lm1 am1L pA_t T maxT pL psi newx
outsheet using C:/Users/X/MSMsim/data/stata`f'.csv, comma replace
clear
clear matrix
clear mata
}