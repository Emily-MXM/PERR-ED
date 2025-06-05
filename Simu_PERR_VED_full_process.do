** PERR methods for recurrent events with Event dependency (ED)
** Varying ED effect
** 500 replicates

* Full 3*2*2*2 simulation scenarios:
* ka = 0.8 or 1.0 or 1.2
* phi = -1.0 or 1.0
* HR = beta = 0.5 or 2.0
* n = 600 or 1000

clear 
set more off
version 15.1

cd "C:\Users\gmsmxm\NUS Dropbox\Xiangmei Ma\PERR2024\EventDependency\Github"

capture log close
log using SimResults_PERR_VED,replace text

run sttocc_rev.ado
run PERR_VED_program.ado

global reps = 500

disp c(current_time) "  " c(current_date)

foreach n in 600 1000 {
	foreach hr in 0.5 2.0 {
	
		***************************** ka = 0.8 *****************************
		local ka = 0.8
		*********************************************		
		simulate,reps($reps) seed(20250226) dots: ///
			PERR_VED_Program, n(`n') hr(`hr') beta0(0.5) ka(`ka') phi(-1.0)
		su

		duplicates report

		foreach x in  AG  CF  {
			gen hr_`x' = exp(coef_`x')
			qui su hr_`x'
			di "R.bias(HR)_`x'  = " round((r(mean)/`hr'-1)*100,0.1) "%"

			gen mse_`x'=(coef_`x'-log(`hr'))^2
			qui su mse_`x', det
			di "RMSE_`x'  = " round(sqrt(r(mean)), 0.001)
			
			qui su cp_`x'
			di "CP_`x'  = " r(mean)*100 "%"
			
			}
			
		disp c(current_time) "  " c(current_date)

		********************************************
		simulate,reps($reps) seed(20250226) dots: ///
			PERR_VED_Program, n(`n') hr(`hr') beta0(-1.5) ka(`ka') phi(1.0)
		su

		duplicates report

		foreach x in  AG  CF  {
			gen hr_`x' = exp(coef_`x')
			qui su hr_`x'
			di "R.bias(HR)_`x'  = " round((r(mean)/`hr'-1)*100,0.1) "%"

			gen mse_`x'=(coef_`x'-log(`hr'))^2
			qui su mse_`x', det
			di "RMSE_`x'  = " round(sqrt(r(mean)), 0.001)
			
			qui su cp_`x'
			di "CP_`x'  = " r(mean)*100 "%"
			
			}
			
		disp c(current_time) "  " c(current_date)

		***************************** ka = 1.0 *****************************
		local ka = 1.0
		*********************************************
		simulate,reps($reps) seed(20250226) dots: ///
			PERR_VED_Program, n(`n') hr(`hr') beta0(0.4) ka(`ka') phi(-1.0)
		su

		duplicates report

		foreach x in  AG  CF  {
			gen hr_`x' = exp(coef_`x')
			su hr_`x'
			di "R.bias(HR)_`x'  = " round((r(mean)/`hr'-1)*100,0.1) "%"

			gen mse_`x'=(coef_`x'-log(`hr'))^2
			qui su mse_`x', det
			di "RMSE_`x'  = " round(sqrt(r(mean)), 0.001)
			
			qui su cp_`x'
			di "CP_`x'  = " r(mean)*100 "%"
			
			}
			
		disp c(current_time) "  " c(current_date)

		********************************************
		simulate,reps($reps) seed(20250226) dots: ///
			PERR_VED_Program, n(`n') hr(`hr') beta0(-1.6) ka(`ka') phi(1.0)
		su

		duplicates report

		foreach x in  AG  CF  {
			gen hr_`x' = exp(coef_`x')
			su hr_`x'
			di "R.bias(HR)_`x'  = " round((r(mean)/`hr'-1)*100,0.1) "%"

			gen mse_`x'=(coef_`x'-log(`hr'))^2
			qui su mse_`x', det
			di "RMSE_`x'  = " round(sqrt(r(mean)), 0.001)
			
			qui su cp_`x'
			di "CP_`x'  = " r(mean)*100 "%"
			
			}
			
		disp c(current_time) "  " c(current_date)

		***************************** ka = 1.2 *****************************
		local ka = 1.2
		*********************************************
		simulate,reps($reps) seed(20250226) dots: ///
			PERR_VED_Program, n(`n') hr(`hr') beta0(0.3) ka(`ka') phi(-1.0)
		su

		duplicates report

		foreach x in  AG  CF  {
			gen hr_`x' = exp(coef_`x')
			su hr_`x'
			di "R.bias(HR)_`x'  = " round((r(mean)/`hr'-1)*100,0.1) "%"

			gen mse_`x'=(coef_`x'-log(`hr'))^2
			qui su mse_`x', det
			di "RMSE_`x'  = " round(sqrt(r(mean)), 0.001)
			
			qui su cp_`x'
			di "CP_`x'  = " r(mean)*100 "%"
			
			}
			
		disp c(current_time) "  " c(current_date)

		********************************************
		simulate,reps($reps) seed(20250226) dots: ///
			PERR_VED_Program, n(`n') hr(`hr') beta0(-1.7) ka(`ka') phi(1.0)
		su

		duplicates report

		foreach x in  AG  CF  {
			gen hr_`x' = exp(coef_`x')
			su hr_`x'
			di "R.bias(HR)_`x' = " round((r(mean)/`hr'-1)*100,0.1) "%"

			gen mse_`x'=(coef_`x'-log(`hr'))^2
			qui su mse_`x', det
			di "RMSE_`x' = " round(sqrt(r(mean)), 0.001)
			
			qui su cp_`x'
			di "CP_`x' = " r(mean)*100 "%"
			
			}
			
		disp c(current_time) "  " c(current_date)
		
		*********************************************
	}
}
*the end!
*********************************************
log close
