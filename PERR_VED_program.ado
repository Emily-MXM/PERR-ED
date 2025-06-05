** 04 Jun 2025
** Simulation Settings for recurrent events with Event dependency (ED)
** PERR method with Conditional Frailty (CF) Model
** Varying ED effect
** Data generation using thinning method

clear
set more off

capture program drop PERR_VED_Program
program define PERR_VED_Program, rclass

version 15.1
 syntax [, n(integer 600) var(real 0.5) beta0(real 0.0) hr(real 2.0) ka(real 0.8)   ///
		   phi(real -1.0) maxE(integer 200)  ]

drop _all
set obs `n'

gen id = _n

*follow-up time in years
gen tau = runiform(1,3)
gen entry = 0.001

*time constant covariates
gen z = rbinomial(1,0.5)
gen x = rbinomial(1,0.5)

gen exit = tau

*time to receive treatment
*hazard of treatment: Weibull disctribution with kaC=0.5
local kaC = 0.5
local c0 = -2.0
gen lamC = exp(`c0'+0.5*x+0.5*z)*rgamma(1/0.1,0.1)
*fixed treatment time
gen ct = (-ln(runiform())/lamC)^(1/`kaC')
replace ct = exit if ct>exit
replace ct = entry if ct<entry

*hazard of recurrent events: Weibull disctribution with ka>0
local beta = log(`hr')
local phi2 = 0.5

gen w = rgamma(1/`var',`var')
gen sup = cond(`ka'>1,tau,0.001)
gen lam = w*exp(`beta0'+0.5*x+0.5*z)

gen t0 = entry
gen t_star = entry
gen event0 = 0
gen eventN0 = 0
gen eventN0_effect = 0

qui{
forvalues i = 1(1)`maxE'{

local j = `i'-1
gen v = runiform()
gen rate = `ka'*sup^(`ka'-1)*lam*exp(cond(`beta'>0,`beta',0)+cond(`phi'>0, `phi',0)*eventN`j')

gen u = rexponential(1/rate)  /* scale = 1/rate  */
replace t_star = t_star + u
gen t`i' = t`j'

gen Et= eventN`j'_effect
replace t`i' = t_star if v <= `ka'*(t_star)^(`ka'-1)*lam*exp(`beta'*(ct<exit)*(t_star>ct)+Et)/rate

gen event`i' = (t`i' > t`j' & t`i'<exit & t`i' != ct)
gen eventN`i' = eventN`j'+event`i'*1
gen eventN`i'_effect = eventN`j'_effect+event`i'*`phi'*exp(-`phi2'*t`i')

replace t`i' = exit if t`i'>exit 
drop u v Et rate
	}
}
replace t`maxE' = exit

drop lam* w tau sup eventN* t_star
drop event*

cap drop trt
gen trt = 0
replace trt = 1 if ct<exit

reshape long t ,i(id) j(temp)
drop temp 

duplicates drop
drop if t<=entry

drop if t==ct & trt==1

gen event = 1
replace event = 0 if  t==exit 
replace event = 0 if  t==entry 

by id (t),sort: gen n = _n
by id (t),sort: gen N = _N

gen ep = 1
replace ep = 2 if n==N
expand ep
drop n N
sort id t ep
by id t ep, sort: gen n = _n
by id t ep, sort: replace t = ct if n==2
drop n ep

drop if t>exit
duplicates drop

gen start = entry
gen end = t
by id (t), sort: replace start = end[_n-1] if _n>1
drop t

********************************************
** matched dataset for PERR method

sort id start end
by id: gen n = _n

reshape wide start end event, i(id) j(n)

gen cut = ct
replace cut = exit if cut>=exit

gen entry0 = 0
stset cut, failure(trt) enter(entry0) 

gen inclusive = (trt==0)
*match on z only:
qui sttocc_rev, match(z) n(1) nodots
gen index = ct
by _set (_case), sort: replace index = index[_N]

by _set, sort: gen Nset = _N
drop if Nset==1

drop Nset _set _time _case  cut inclusive entry0  

reshape long start end event, i(id) j(temp)
drop temp
drop if event==.

expand 2 if trt==0

sort id start end
by id start: gen n = _n
replace end = index if n==1 & end>index & start<index & trt==0
replace start = index if n==2 & end>index & start<index & trt==0
drop n
duplicates drop

replace event = 0 if end==index & trt==0

cap drop n
by id (start),sort: gen n = _n

count if n==1
return scalar pop=r(N)

count if event==1
return scalar totalE=r(N)
drop n

*PERR analysis from here:	
gen period = cond(end<=index,0,1)
gen pt = trt*period

*PERR with AG model
stset end, failure(event) exit(t .) id(id) enter(start) 
stcox  pt trt period , nohr cluster(id)

matrix table_AG = r(table)
return scalar coef_AG = table_AG[1,1]
return scalar se_AG = table_AG[2,1]
return scalar pvalue_AG = table_AG[4,1]
return scalar cp_AG = cond(`beta'>=table_AG[5,1]&`beta'<=table_AG[6,1],1,0)

*PERR from CF model

stgen event_order = nfailures()
replace event_order = event_order+1

by id, sort: gen n=_n
by id, sort: egen event_orderN=sum(event)
replace event_orderN=event_orderN+1

qui su event_orderN if n==1,det
local last = int(r(p95))

drop n event_orderN

gen r=runiform(1,_N)

local seednum = _N+int(r[1])
drop r

set matsize 800
strmcure pt trt period , shared(id) strata(event_order) lastpool(`last') seed(`seednum')

matrix table_CF = r(table)
return scalar coef_CF = table_CF[1,1]
return scalar frailty_CF = table_CF[1,4]
return scalar se_CF = table_CF[2,1]
return scalar pvalue_CF = table_CF[4,1]
return scalar cp_CF = cond(`beta'>=table_CF[5,1]&`beta'<=table_CF[6,1],1,0)

eret clear

end

