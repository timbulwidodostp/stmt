*! version 1.4.4 10Mar2012

/*
History
PL 12Feb2012: added failure option (1-S(t)).
PR 09sep2011: rmst now works with factor variables.
PL 08sep2011: fixed bug in centiles for uncured in cure models / PR rmst update
PL 15aug2011: if timevar variable = 0 then survival is predicted to be 1 rather than missing.
PR 26may2011: undocumented option nliter() to set option iterate() in predictnl for rmst.
		If not set, predictnl's default iterate(100) is used. For simulations, it may
		be good to take nliter(0) or nliter(some small integer) to avoid long run times.
PR 09apr2011: hidden option power() of _rmst subroutine not implemented as option for predict
PR 25feb2011: add rmst and abc options for mean survival time and area between curves
PR 25feb2011: add stdp option to various predict options
PL 29NOV2010: fixed bug when using meansurv option with a null model.
PL 07Nov2010: corrected bug that caused error if tempory folder had spaces.
PL 01Nov2010: introduced factor variables.
PL 02sep2010: added Therese Anderssons additions for prediction after fitting a cure model
PL 15mar2010: added reverse option to spline generation
PL 01feb2010: change to Stata 11.
PL 12nov2009: added tvc option
PL 21sep2009: ensure predictions work with new rcsbaselineoff option in stpm2.
PL 08Sep2009: modification of msurvpop to work in Stata 11.
PL 12Mar2009: changed to using Newton-Raphson method for estimating centiles.
PL 11Mar2009: fixed problem with strings>244 characters for long varlists.
PL 11dec2008: changed to using rcsgen for spline functions.
PL 20aug2008: Added hdiff1 and hdiff2 options for difference in survival curves.
PL 20aug2008: Added hdiff1 and hdiff2 options for difference in survival curves.
PL 09aug2008: Fixed bug for some predictions involving time-dependent effects.
*/	

program stpm2_pred
	version 11.1
	syntax newvarname [if] [in], [Survival Hazard XB XBNOBaseline DXB DZDY HRNumerator(string) HRDenominator(string) MEANSurv ///
									CENtile(string) CUMHazard CUMOdds NORmal MARTingale DEViance DENSity AT(string) ZEROs FAILure ///
									noOFFset SDIFF1(string) SDIFF2(string) HDIFF1(string) HDIFF2(string) tvc(varname) ///
									CI LEVel(real `c(level)') TIMEvar(varname) STDP PER(real 1)  ///
									CENTOL(real 0.0001) Cure UNCured STARTUNC(real -1) ///
									RMst RSDst TMAx(string) TMIn(string) n(int 1000) POWer(real 1) ABc hr0(string) first(int 1) ///
									NLiter(string) CENTITER(int 100)]
	marksample touse, novarlist
	local newvarname `varlist'
	qui count if `touse'
	if r(N)==0 {
		error 2000          /* no observations */
	}
	
/* Check Options */
/* First check rmst options */
	if "`tmax'" != "" {
		if "`rmst'`rsdst'`abc'" == "" {
			display as error "tmax() valid only with rmst, rsdst or abc"
			exit 198
		}
		confirm number `tmax'
		if `tmax' <= 0 {
			display as error "tmax() must be positive"
			exit 198
		}
	  if "`rsdst'" != "" local rmst rmst
	}
	else local tmax 0

	if "`tmin'" != "" {
		if "`rmst'`rsdst'`abc'" == "" {
			display as error "tmin() valid only with rmst, rsdst or abc"
			exit 198
		}
		confirm number `tmin'
		if `tmin' < 0 {
			display as error "tmin() may not be negative"
			exit 198
		}
	}
	else local tmin 0
	if "`hr0'" != "" {
		if "`rmst'`abc'" == "" {
			display as error "hr0() valid only with rmst or abc"
			exit 198
		}
		cap confirm var `hr0'
		if c(rc) {
			cap confirm number `hr0'
			if c(rc) {
				display as error "hr0() must be a numeric constant or variable"
				exit 198
			}
			if `hr0' <= 0 {
				display as error "hr0() must exceed 0"
				exit 198
			}
		}
	}
	
	if "`hrdenominator'" != "" & "`hrnumerator'" == "" {
		display as error "You must specifiy the hrnumerator option if you specifiy the hrdenominator option"
		exit 198
	}

	if "`sdiff2'" != "" & "`sdiff1'" == "" {
		display as error "You must specifiy the sdiff1 option if you specifiy the sdiff2 option"
		exit 198
	}

	if "`hdiff2'" != "" & "`hdiff1'" == "" {
		display as error "You must specifiy the hdiff1 option if you specifiy the hdiff2 option"
		exit 198
	}
	
	local hratiotmp = substr("`hrnumerator'",1,1)
	local sdifftmp = substr("`sdiff1'",1,1)
	local hdifftmp = substr("`hdiff1'",1,1)
	if wordcount(`"`survival' `hazard' `failure' `meansurv' `hratiotmp' `sdifftmp' `hdifftmp' `centile' `xb' `xbnobaseline' `dxb' `dzdy' `martingale' `deviance' `cumhazard' `cumodds' `normal' `density' `tvc' `cure' `rmst'"') > 1 {
		display as error "You have specified more than one option for predict"
		exit 198
	}
	if wordcount(`"`survival' `hazard' `failure' `meansurv' `hrnumerator' `sdiff1'  `hdifftmp' `centile' `xb' `xbnobaseline' `dxb' `dzdy' `martingale' `deviance' `cumhazard' `cumodds' `normal' `density' `tvc' `cure'`rmst'`abc'"') == 0 {
		display as error "You must specify one of the predict options"
		exit 198
	}
	
	if `per' != 1 & "`hazard'" == "" & "`hdiff1'" == "" {
		display as error "You can only use the per() option in combinaton with the hazard or hdiff1()/hdiff2() options."
		exit 198		
	}

	if "`stdp'" != "" & "`ci'" != "" {
		display as error "You can not specify both the ci and stdp options."
		exit 19
	}
	
	if "`stdp'" != "" & ///
		wordcount(`"`xb' `dxb' `xbnobaseline' `rmst' `abc' `hrnumerator' `hdiff1' `sdiff1'"') == 0 {
		display as error "The stdp option cannot be used with this predict option."
		exit 198
	}

	if "`ci'" != "" & ///
		wordcount(`"`survival' `hazard' `failure' `hrnumerator' `sdiff1' `hdiff1' `centile' `xb' `dxb' `xbnobaseline' `tvc' `meansurv' `cure' `rmst' `abc'"') == 0 {
		display as error "The ci option can not be used with this predict option."
		exit 198
	}
	
	if "`zeros'" != "" & "`meansurv'" != "" {
		display as error "You can not specify the zero option with the meansurv option."
		exit 198
	}

	if "`zeros'" != "" & "`tvc'" != "" {
		display as error "You can not specify the zero option with the tvc option."
		exit 198
	}

	if "`zeros'" != "" & ("`hrnumerator'" != "" | "`hdiff1'" != "" | "`sdiff1'" != "") {
		display as error "You can not specify the zero option with the hrnumerator, hdiff or sdiff options."
		exit 198
	}

	
	if "`at'" != "" & "`hrnumerator'" != "" {
		display as error "You can not use the at option with the hrnumerator option"
		exit 198
	}

	if "`at'" != "" & "`sdiff1'" != "" {
		display as error "You can not use the at option with the sdiff1 and sdiff2 options"
		exit 198
	}
	
	if "`at'" != "" & "`hdiff1'" != "" {
		display as error "You can not use the at option with the hdiff1 and hdiff2 options"
		exit 198
	}

	if ("`uncured'" != "" | "`cure'" != "") &  "`e(cure)'" == "" {
		display as error "You can only use the cure and uncured options after fitting a cure model"
		exit 198
	}
	
	if "`uncured'" != "" & "`survival'" == "" & "`hazard'" == "" & "`centile'" == "" {
		display as error "You must specifiy the survival, hazard or centile option if you specify uncured"
		exit 198
	}

	if `startunc' != -1 & ("`uncured'" == "" | "`centile'" == "") {
		display as error "You must specifiy the uncured and the centile option if you set a starting value"
		exit 198
	}
	
	if "`meansurv'" != "" & ("`ci'" != "") {
		display as error "You can not use the ci option with the meansurv option"
		exit 198
	}
	
/* call Stmeancurve if meansurv option specified */
	if "`meansurv'" != "" {
		Stmeancurve `newvarname' if `touse', timevar(`timevar') at(`at')  `offset'
		exit
	}

/* call _rmst (formerly stpm2_mttf) for SE and CI if rmst option specified */
	if "`rmst'" != "" {
		if "`at'" != "" local At at(`at')
		cap confirm var `newvarname'
		if c(rc) != 0 {
			 _rmst if `touse', tmin(`tmin') tmax(`tmax') generate(`newvarname') ///
			 `At' `zeros' n(`n') power(`power') `rsdst'			
			 if "`ci'`stdp'" == "" {
				exit
			}
		}
		cap drop `newvarname'
		if "`ci'" != "" {
			local predictnl_opts , ci(`newvarname'_lci `newvarname'_uci)
		}
		else if "`stdp'" != "" {
			local predictnl_opts , se(`newvarname'_se)
		}
		if "`nliter'"!="" local nli iterate(`nliter')
		else local nli
		qui predictnl `newvarname' = predict(rmst tmin(`tmin') tmax(`tmax') ///
		 `At' `zeros' n(`n') power(`power') `rsdst') `nli' `predictnl_opts'
		exit
	}

	if "`abc'" != "" {
		// if "`at'" != "" local At at(`at')
		if "`hrdenominator'" != "" local hrd hrdenominator(`hrdenominator')
		if `first' local wt wt
		else local wt
		if "`tmin'" != "" local tmin tmin(`tmin')
		if "`hr0'" != "" local hr0 hr0(`hr0')

		cap confirm var `newvarname', exact
		if c(rc) != 0 { // `newvarname' does not exist
			qui _abc `newvarname' if `touse', `hr0' hrnumerator(`hrnumerator') `hrd' ///
			 tmax(`tmax') `tmin' `At' `zeros' n(`n') `wt'
			cap drop _weight
			if "`ci'`stdp'" == "" {
				*cap erase _predict_weights.dta
				exit
			}
		}
		cap drop `newvarname'
		if "`ci'" != "" {
			local predictnl_opts ci(`newvarname'_lci `newvarname'_uci)
		}
		else if "`stdp'" != "" {
			local predictnl_opts se(`newvarname'_se)
		}
		cap predictnl double `newvarname' = predict(abc `hr0' hrnumerator(`hrnumerator') `hrd' ///
		 tmax(`tmax') `tmin' `At' `zeros' n(`n') first(0)) if `touse', `predictnl_opts' level(`level')
		local rc = c(rc)
		// Tidy up weight variable and temporary file created by _abc
		cap drop _weight
		cap erase _predict_weights.dta
		if `rc' {
			noi di as err "could not estimate SE or CI - predictnl failed to converge"
			exit 498
		}
		else exit
	}

	

/* calculate midt for centile option */
	summ _t, meanonly
	local midt = (r(max) - r(min))/2
	
/* calculate startt for centile option if uncured option is specified*/

	if `startunc' == -1 {
		summ _t, meanonly
	 	local startt = (r(max) - r(min))/8
	}
	else {
		local startt = `startunc'
	}
	
/* store time-dependent covariates and main varlist */
	local etvc `e(tvc)'
	local main_varlist `e(varlist)'
		
/* dydx option of old version of stpm */
	if "`dzdy'" != "" {
		local dxb dxb
	}
/* generate ocons for use when orthogonalising splines */
	tempvar ocons
	gen `ocons' = 1
	
/* Use _t if option timevar not specified */
	tempvar t lnt 
	if "`timevar'" == "" {
		qui gen double `t' = _t if `touse'
		qui gen double `lnt' = ln(_t) if `touse'
	}
	else {
		qui gen double `t' = `timevar' if `touse'
		qui gen double `lnt' = ln(`timevar') if `touse'
	}
	
/* Check to see if nonconstant option used */
	if "`e(noconstant)'" == "" {
		tempvar cons
		qui gen `cons' = 1 if `touse'
	}	

/* Preserve data for out of sample prediction  */	
	tempfile newvars 
	preserve

/* Calculate new spline terms if timevar option specified */
	if "`timevar'" != "" & "`e(rcsbaseoff)'" == "" {
		capture drop _rcs* _d_rcs*
		if "`e(orthog)'" != "" {
			tempname rmatrix
			matrix `rmatrix' = e(R_bh)
			local rmatrixopt rmatrix(`rmatrix')
		}
		qui rcsgen `lnt' if `touse', knots(`e(ln_bhknots)') gen(_rcs) dgen(_d_rcs) `e(reverse)' `rmatrixopt'
	}
	
/* calculate new spline terms if timevar option or hrnumerator option is specified */

	if "`timevar'" != "" | "`hrnumerator'" != "" | "`sdiff1'" != "" | "`hdiff1'" != "" {
		foreach tvcvar in `e(tvc)' {
			if (("`hrnumerator'" != "" | "`sdiff1'" != "" | "`hdiff1'" != "") & "`timevar'" == "") | "`e(rcsbaseoff)'" != "" {
				capture drop _rcs_`tvcvar'* _d_rcs_`tvcvar'*
			}
			if "`e(orthog)'" != "" {
				tempname rmatrix_`tvcvar'
				matrix `rmatrix_`tvcvar'' = e(R_`tvcvar')
				local rmatrixopt rmatrix(`rmatrix_`tvcvar'')
			}
			qui rcsgen `lnt' if `touse',  gen(_rcs_`tvcvar') knots(`e(ln_tvcknots_`tvcvar')') dgen(_d_rcs_`tvcvar') `e(reverse)' `rmatrixopt'
			if "`hrnumerator'" == "" & "`sdiff1'"  == "" & "`hdiff1'" == "" {
				forvalues i = 1/`e(df_`tvcvar')'{
					qui replace _rcs_`tvcvar'`i' = _rcs_`tvcvar'`i'*`tvcvar' if `touse'
					qui replace _d_rcs_`tvcvar'`i' = _d_rcs_`tvcvar'`i'*`tvcvar' if `touse'
				}
			}

		}
	}	
	
/* zeros */
	if "`zeros'" != "" {
		local tmptvc `e(tvc)'
		foreach var in `e(varlist)' {
			_ms_parse_parts `var'
			if `"`: list posof `"`r(name)'"' in at'"' == "0" { 
				qui replace `r(name)' = 0 if `touse'
				if `"`: list posof `"`r(name)'"' in tmptvc'"' != "0" { 
				forvalues i = 1/`e(df_`r(name)')' {
						qui replace _rcs_`r(name)'`i' = 0 if `touse'
						qui replace _d_rcs_`r(name)'`i' = 0 if `touse'
					}
				}
			}
		}
	}

/* Out of sample predictions using at() */
	if "`at'" != "" {
		tokenize `at'
		while "`1'"!="" {
			fvunab tmpfv: `1'
			local 1 `tmpfv'
			_ms_parse_parts `1'
			if "`r(type)'"!="variable" {
				display as error "level indicators of factor" /*
								*/ " variables may not be individually set" /*
								*/ " with the at() option; set one value" /*
								*/ " for the entire factor variable"
				exit 198
			}
			cap confirm var `2'
			if _rc {
				cap confirm num `2'
				if _rc {
					di as err "invalid at(... `1' `2' ...)"
					exit 198
				}
			}
			qui replace `1' = `2' if `touse'
			if `"`: list posof `"`1'"' in etvc'"' != "0" {
				local tvcvar `1'
				if "`e(orthog)'" != "" {
					tempname rmatrix_`tvcvar'
					matrix `rmatrix_`tvcvar'' = e(R_`tvcvar')
					local rmatrixopt rmatrix(`rmatrix_`tvcvar'')
				}
				capture drop _rcs_`tvcvar'* _d_rcs_`tvcvar'*
				qui rcsgen `lnt' if `touse', knots(`e(ln_tvcknots_`tvcvar')') gen(_rcs_`tvcvar') dgen(_d_rcs_`tvcvar') `e(reverse)' `rmatrixopt'
				forvalues i = 1/`e(df_`tvcvar')'{
					qui replace _rcs_`tvcvar'`i' = _rcs_`tvcvar'`i'*`tvcvar' if `touse'
					qui replace _d_rcs_`tvcvar'`i' = _d_rcs_`tvcvar'`i'*`tvcvar' if `touse'
				}
			}
			mac shift 2
		}
	}

/* Add offset term if exists unless no offset option is specified */
	if "`e(offset1)'" !=  "" & /* !! PR */ "`offset'" != "nooffset" {
		local addoff "+ `e(offset1)'" 
	}

/* check ci and stdp options */
	if "`ci'" != "" & "`stdp'" != "" {
		display as error "Only one of the ci and se options can be specified"
		exit 198
	}
	
/* Deviance and Martingale Residuals */
	if "`deviance'" != "" | "`martingale'" != "" {
		tempvar cH res
		qui predict `cH' if `touse', cumhazard timevar(`t')  `offset'
		gen double `res' = _d - `cH' if `touse'
		if "`deviance'" != "" {
			gen double `newvarname' = sign(`res')*sqrt( -2*(`res' + _d*(ln(_d -`res')))) if `touse'
        }
        else rename `res' `newvarname'
	}
	
/* Failure (1-S(t)) */
	else if "`failure'" != "" {
		qui predict `newvarname' if `touse', s timevar(`t')  `offset' ci
		qui replace `newvarname' = 1 - `newvarname' if `touse'
		tempvar tmpSt
		qui gen double `tmpSt' = 1 - `newvarname'_uci if `touse'
		qui replace `newvarname'_uci = 1 - `newvarname'_lci if `touse'
		qui replace `newvarname'_lci = `tmpSt' if `touse'
	}
	
/* Cumulative Hazard */
	else if "`cumhazard'" != "" {
		tempvar S
		predict `S' if `touse', s timevar(`t')  `offset'
		gen double `newvarname' = -ln(`S') if `touse'
	}

/* Cumulative Odds */
	else if "`cumodds'" != "" {
		tempvar S
		predict `S' if `touse', s timevar(`t') /* !! pr */ `offset'
		gen double `newvarname' = (1 -`S')/`S' if `touse'
	}
	
/* Standard Normal Deviate */
	else if "`normal'" != "" {
		tempvar S
		predict `S' if `touse', s timevar(`t') /* !! pr */ `offset'
		gen double `newvarname' = -invnormal(`S') if `touse'
	}
	
/* density */
	else if "`density'" != "" {
		tempvar S h
		predict  `S' if `touse', s timevar(`t') /* !! pr */ `offset'
		predict  `h' if `touse', h timevar(`t') /* !! pr */ `offset'
		gen double `newvarname' = `S'*`h' if `touse'
	}	
	
/* linear predictor */	
	else if "`xb'" != "" {
		if "`ci'" != "" {
			local prednlopt ci(`newvarname'_lci `newvarname'_uci)
		}
		else if "`stdp'" != "" {
			local prednlopt se(`newvarname'_se)
		}
		qui predictnl double `newvarname' = xb(xb) `addoff' if `touse', `prednlopt' level(`level')
	}
			
/* derivative of linear predictor */	
	else if "`dxb'" != "" {
		if "`ci'" != "" {
			local prednlopt ci(`newvarname'_lci `newvarname'_uci)
		}
		else if "`stdp'" != "" {
			local prednlopt se(`newvarname'_se)
		}
		qui predictnl double `newvarname' = xb(dxb) if `touse', `prednlopt' level(`level')
	}
	
/* linear predictor exluding spline terms */
	else if "`xbnobaseline'" != "" {
		if "`ci'" != "" {
			local prednlopt ci(`newvarname'_lci `newvarname'_uci)
		}
		else if "`stdp'" != "" {
			local prednlopt se(`newvarname'_se)
		}
/* commented out for now - ignores the constant (gamma_0) - may be needed later
		if "`e(noconstant)'" == "" {	
			local xbnobhpred [xb][_cons]
		}	
*/
		foreach var in `e(varlist)' {
			if "`xbnobhpred'" == "" {
				local xbnobhpred [xb][`var']*`var'
			}
			else {
				local xbnobhpred `xbnobhpred' + [xb][`var']*`var'
			}
			if `"`: list posof `"`var'"' in etvc'"' != "0" {
				forvalues i = 1/`e(df_`var')' {
					local xbnobhpred `xbnobhpred' + [xb][_rcs_`var'`i']*_rcs_`var'`i'
				}
			}
		}
*		if "`e(noconstant)'" != "" {	
*			local xbnobhpred = subinstr("`xbnobhpred'","+","",1) 
*		}
		predictnl double `newvarname' = `xbnobhpred' if `touse', `prednlopt' level(`level')
	}
	
/* tvc option */
	else if "`tvc'" != "" {
		if "`ci'" != "" {
			local prednlopt ci(`newvarname'_lci `newvarname'_uci)
		}
		else if "`stdp'" != "" {
			local prednlopt se(`newvarname'_se)
		}
		local tvcpred [xb][`tvc']*`tvc'
		if `"`: list posof `"`tvc'"' in etvc'"' != "0" {
			forvalues i = 1/`e(df_`tvc')' {
				local tvcpred `tvcpred' + [xb][_rcs_`tvc'`i']*_rcs_`tvc'`i'
			}
		}
*		predictnl double `newvarname' = (`tvcpred')/`tvc' if `touse', `prednlopt' level(`level')
		predictnl double `newvarname' = (`tvcpred') if `touse', `prednlopt' level(`level')
	}
		
/* Survival Function */
	else if "`survival'" != "" & "`uncured'" == "" {
		tempvar sxb 
		if "`ci'" != "" {
			tempvar sxb_lci sxb_uci
			local prednlopt ci(`sxb_lci' `sxb_uci')
		}
		if "`e(scale)'" != "theta" {
			qui predictnl double `sxb' = xb(xb) `addoff' if `touse', `prednlopt' level(`level') 
		}
/* predict on ln(-ln S(t)) scale for theta */
		else if "`e(scale)'" == "theta" {
			qui predictnl double `sxb' = ln(ln(exp(xb(ln_theta))*exp(xb(xb)`addoff')+1)/exp(xb(ln_theta))) if `touse', `prednlopt'  level(`level') 		
		}
/* Transform back to survival scale */
		if "`e(scale)'" == "hazard" {
			qui gen double `newvarname' = exp(-exp(`sxb')) if `touse'
			if "`ci'" != "" {
				qui gen `newvarname'_lci = exp(-exp(`sxb_uci'))  if `touse'
				qui gen `newvarname'_uci =  exp(-exp(`sxb_lci')) if `touse'
			}
		}
		else if "`e(scale)'" == "odds" {
			qui gen double `newvarname' = (1 +exp(`sxb'))^(-1) if `touse'
			if "`ci'" != "" {
				qui gen `newvarname'_lci = (1 +exp(`sxb_uci'))^(-1) if `touse'
				qui gen `newvarname'_uci = (1 +exp(`sxb_lci'))^(-1) if `touse'
			}
		}
		else if "`e(scale)'" == "normal" {
			qui gen double `newvarname' = normal(-`sxb') if `touse'
			if "`ci'" != "" {
				qui gen `newvarname'_lci = normal(-`sxb_uci') if `touse'
				qui gen `newvarname'_uci = normal(-`sxb_lci') if `touse' 
			}
		}		
		else if "`e(scale)'" == "theta" {
			qui gen double `newvarname' = exp(-exp(`sxb')) if `touse'
			if "`ci'" != "" {
				qui gen `newvarname'_lci = exp(-exp(`sxb_lci')) if `touse'
				qui gen `newvarname'_uci = exp(-exp(`sxb_uci')) if `touse' 
			}
		}
		qui replace `newvarname' = 1 if `t' == 0 & `touse'
		if "`ci'" != "" {
			qui replace `newvarname'_lci = 1 if `t' == 0 & `touse'
			qui replace `newvarname'_uci = 1 if `t' == 0 & `touse'
		}
	}
	

/* Hazard Function */

	else if "`hazard'" != "" & "`uncured'" == "" {
		tempvar lnh 
		if "`ci'" != "" {
			tempvar lnh_lci lnh_uci
			local prednlopt ci(`lnh_lci' `lnh_uci')
		}
		if "`e(scale)'" == "hazard" {
			qui predictnl double `lnh' = -ln(`t') + ln(xb(dxb)) + xb(xb) `addoff'  if `touse', `prednlopt' level(`level') 
		}
		if "`e(scale)'" == "odds" {
			qui predictnl double `lnh' = -ln(`t') + ln(xb(dxb)) + (xb(xb)`addoff')  -ln(1+exp(xb(xb)`addoff'))   if `touse', `prednlopt' level(`level') 
		}		
		if "`e(scale)'" == "normal" {
			qui predictnl double `lnh' = -ln(`t') + ln(xb(dxb)) + ln(normalden(xb(xb)`addoff')) - ln(normal(-(xb(xb)`addoff')))   if `touse', `prednlopt' level(`level') 
		}		
		if "`e(scale)'" == "theta" {
			qui predictnl double `lnh' = -ln(`t') + ln(xb(dxb)) + xb(`xb') - ln(exp(xb(ln_theta))*exp(xb(xb)`addoff') + 1)  if `touse', `prednlopt' level(`level') 
		}		

/* Transform back to hazard scale */
		qui gen double `newvarname' = exp(`lnh')*`per' if `touse'
		if "`ci'" != "" {
			qui gen `newvarname'_lci = exp(`lnh_lci')*`per'  if `touse'
			qui gen `newvarname'_uci =  exp(`lnh_uci')*`per' if `touse'
		}
	}
	
/* Predict Hazard Ratio */
	else if "`hrnumerator'" != "" {
		tempvar lhr
		if `"`ci'"' != "" {
			tempvar lhr_lci lhr_uci
			local predictnl_opts ci(`lhr_lci' `lhr_uci')
		}
		else if "`stdp'" != "" {
			tempvar lhr_se
			local predictnl_opts se(`lhr_se')
		}		
		
		forvalues i=1/`e(dfbase)' {
			local dxb1 `dxb1' [xb][_rcs`i']*_d_rcs`i' 
			local dxb0 `dxb0' [xb][_rcs`i']*_d_rcs`i'
			local xb1_plus `xb1_plus' [xb][_rcs`i']*_rcs`i'
			local xb0_plus `xb0_plus' [xb][_rcs`i']*_rcs`i'
			if `i' != `e(dfbase)' {
				local dxb0 `dxb0' + 
				local dxb1 `dxb1' + 
				local xb1_plus `xb1_plus' +
				local xb0_plus `xb0_plus' +
			}
		}
		
/* use Parse_list to select appropriate values of factor variables */
		Parse_list, listname(hrnumerator) parselist(`hrnumerator')
		tokenize `r(retlist)'
		while "`1'"!="" {
			cap confirm var `2'
			if _rc {
				if "`2'" == "." {
					local 2 `1'
				}
				else {
					cap confirm num `2'
					if _rc {
						di in red "invalid hrnumerator(... `1' `2' ...)"
						exit 198
					}
				}
			}
			if "`xb10'" != "" & "`2'" != "0" & `: list posof `"`1'"' in main_varlist' != 0{
				local xb10 `xb10' +
			}
			if "`xb1_plus'" != "" & "`2'" != "0" &  `: list posof `"`1'"' in main_varlist' != 0{
				local xb1_plus `xb1_plus' +
			}
			if "`2'" != "0" & `: list posof `"`1'"' in main_varlist' != 0 {
				local xb10 `xb10' [xb][`1']*`2' 
				local xb1_plus `xb1_plus' [xb][`1']*`2' 
			}
			if `"`: list posof `"`1'"' in etvc'"' != "0" & "`2'" != "0" {
				if "`e(rcsbaseoff)'" == "" {
					local dxb1 `dxb1' +
				}
				local xb10 `xb10' +
				local xb1_plus `xb1_plus' +

				forvalues i=1/`e(df_`1')' {
					local dxb1 `dxb1' [xb][_rcs_`1'`i']*_d_rcs_`1'`i'*`2' 
					local xb10 `xb10' [xb][_rcs_`1'`i']*_rcs_`1'`i'*`2'  
					local xb1_plus `xb1_plus' [xb][_rcs_`1'`i']*_rcs_`1'`i'*`2'  
					if `i' != `e(df_`1')' {
						local dxb1 `dxb1' +
						local xb10 `xb10' +
						local xb1_plus `xb1_plus' +
					}
				}
			}
			mac shift 2
		}			
			
		if "`hrdenominator'" != "" {

/* use Parse_list to select appropriate values of factor variables */
		Parse_list, listname(hrdenominator) parselist(`hrdenominator')
		tokenize `r(retlist)'
			while "`1'"!="" {
				cap confirm var `2'
				if _rc {
					if "`2'" == "." {
						local 2 `1'
					}
					else {
						cap confirm num `2'
						if _rc {
							di as err "invalid hrdenominator(... `1' `2' ...)"
							exit 198
						}
					}
				}
				if "`2'" != "0" & `: list posof `"`1'"' in main_varlist' != 0 {
					local xb10 `xb10' - [xb][`1']*`2'
					if "`e(rcsbaseoff)'" == "" & `: list posof `"`1'"' in main_varlist' != 0 {
						local xb0_plus `xb0_plus' + [xb][`1']*`2' 
					}
					else if `: list posof `"`1'"' in main_varlist' != 0 {
						local xb0_plus `xb0_plus' [xb][`1']*`2' 
					}
				}
				if `"`: list posof `"`1'"' in etvc'"' != "0" & "`2'" != "0" {
					if "`e(rcsbaseoff)'" == "" {
						local dxb0 `dxb0' +
					}
					local xb0_plus `xb0_plus' + 
					local xb10 `xb10' - 
					forvalues i=1/`e(df_`1')' {
						local dxb0 `dxb0' [xb][_rcs_`1'`i']*_d_rcs_`1'`i'*`2'
						local xb10 `xb10' [xb][_rcs_`1'`i']*_rcs_`1'`i'*`2'
						local xb0_plus `xb0_plus' [xb][_rcs_`1'`i']*_rcs_`1'`i'*`2'
						if `i' != `e(df_`1')' {
							local dxb0 `dxb0' +
							local xb10 `xb10' -
							local xb0_plus `xb0_plus' +
						}
					}
				}
				mac shift 2
			}
		}
		if "`e(noconstant)'" == "" {
			local xb0_plus `xb0_plus' + [xb][_cons]
			local xb1_plus `xb1_plus' + [xb][_cons]
		}
		
		if "`e(scale)'" =="hazard" {
			qui predictnl double `lhr' = ln(`dxb1') - ln(`dxb0') + `xb10' if `touse', `predictnl_opts' level(`level')
		}
		else if "`e(scale)'" =="odds" {
			qui predictnl double `lhr' =  	ln(`dxb1') - ln(`dxb0') + `xb10' - ///
											ln(1+exp(`xb1_plus')) + ln(1+exp(`xb0_plus')) ///
											if `touse', `predictnl_opts' level(`level')
		}
		else if "`e(scale)'" =="normal" {
			qui predictnl double `lhr' =  	ln(`dxb1') - ln(`dxb0') + ///
											ln(normalden(`xb1_plus')) - ln(normalden(`xb0_plus')) - ///
											ln(normal(-(`xb1_plus'))) + ln(normal(-(`xb0_plus'))) ///
											if `touse', `predictnl_opts' level(`level')
		}
		else if "`e(scale)'" =="theta" {
			qui predictnl double `lhr' =  	ln(`dxb1') - ln(`dxb0') + `xb10' ///
											-ln(exp(xb(ln_theta))*exp(`xb1_plus') + 1) + ln(exp(xb(ln_theta))*exp(`xb0_plus') + 1) ///
											if `touse', `predictnl_opts' level(`level')
		}
		qui gen double `newvarname' = exp(`lhr') if `touse'
		if `"`ci'"' != "" {
			qui gen double `newvarname'_lci=exp(`lhr_lci')  if `touse'
			qui gen double `newvarname'_uci=exp(`lhr_uci')  if `touse'
		}
		else if "`stdp'" != "" {
			qui gen double `newvarname'_se = `lhr_se' * `newvarname'
		}
	}

/* Predict Difference in Hazard Functions */
	else if "`hdiff1'" != "" {
		if `"`ci'"' != "" {
			local predictnl_opts "ci(`newvarname'_lci `newvarname'_uci)"
		}
		else if "`stdp'" != "" {
			local predictnl_opts se(`newvarname'_se)
		}
		
		forvalues i=1/`e(dfbase)' {
			local dxb1 `dxb1' [xb][_rcs`i']*_d_rcs`i' 
			local dxb0 `dxb0' [xb][_rcs`i']*_d_rcs`i'
			local xb1_plus `xb1_plus' [xb][_rcs`i']*_rcs`i'
			local xb0_plus `xb0_plus' [xb][_rcs`i']*_rcs`i'
			if `i' != `e(dfbase)' {
				local dxb0 `dxb0' + 
				local dxb1 `dxb1' + 
				local xb1_plus `xb1_plus' +
				local xb0_plus `xb0_plus' +
			}
		}
/* use Parse_list to select appropriate values of factor variables */
		Parse_list, listname(hdiff1) parselist(`hdiff1')
		tokenize `r(retlist)'
		while "`1'"!="" {
			cap confirm var `2'
			if _rc {
				if "`2'" == "." {
					local 2 `1'
				}
				else {
					cap confirm num `2'
					if _rc {
						di as err "invalid hdiff1(... `1' `2' ...)"
						exit 198
					}
				}
			}
			if "`xb1_plus'" != "" & "`2'" != "0" & `: list posof `"`1'"' in main_varlist' != 0 {
				local xb1_plus `xb1_plus' +
			}
			if "`2'" != "0" & `: list posof `"`1'"' in main_varlist'!=0 {
				local xb1_plus `xb1_plus' [xb][`1']*`2' 
			}
			if `"`: list posof `"`1'"' in etvc'"' != "0" & "`2'" != "0" {
				if "`e(rcsbaseoff)'"  == "" {
					local dxb1 `dxb1' +
				}
				local xb1_plus `xb1_plus' +
				local xb10 `xb10' +

				forvalues i=1/`e(df_`1')' {
					local dxb1 `dxb1' [xb][_rcs_`1'`i']*_d_rcs_`1'`i'*`2' 
					local xb1_plus `xb1_plus' [xb][_rcs_`1'`i']*_rcs_`1'`i'*`2'  
					if `i' != `e(df_`1')' {
						local dxb1 `dxb1' +
						local xb1_plus `xb1_plus' +
					}
				}
			}
			mac shift 2
		}			
		if "`hdiff2'" != "" {
/* use Parse_list to select appropriate values of factor variables */
		Parse_list, listname(hdiff2) parselist(`hdiff2')
		tokenize `r(retlist)'
			while "`1'"!="" {
				cap confirm var `2'
				if _rc {
					if "`2'" == "." {
						local 2 `1'
					}
					else {
						cap confirm num `2'
						if _rc {
							di as err "invalid hdiff2(... `1' `2' ...)"
							exit 198
						}
					}
				}
				if "`2'" != "0" {
					if "`e(rcsbaseoff)'" == "" & `: list posof `"`1'"' in main_varlist' != 0 {
						local xb0_plus `xb0_plus' + [xb][`1']*`2' 
					}
					else if `: list posof `"`1'"' in main_varlist' != 0 {
						local xb0_plus `xb0_plus' [xb][`1']*`2' 
					}
				}
				if `"`: list posof `"`1'"' in etvc'"' != "0" & "`2'" != "0" {
					if "`e(rcsbaseoff)'" == "" {
						local dxb0 `dxb0' +
					}
					local xb0_plus `xb0_plus' + 
					forvalues i=1/`e(df_`1')' {
						local dxb0 `dxb0' [xb][_rcs_`1'`i']*_d_rcs_`1'`i'*`2'
						local xb0_plus `xb0_plus' [xb][_rcs_`1'`i']*_rcs_`1'`i'*`2'
						if `i' != `e(df_`1')' {
							local dxb0 `dxb0' +
							local xb0_plus `xb0_plus' +
						}
					}
				}
				mac shift 2
			}
		}
		if "`e(noconstant)'" == "" {
			local xb0_plus `xb0_plus' + [xb][_cons]
			local xb1_plus `xb1_plus' + [xb][_cons]
		}
		if "`e(scale)'" =="hazard" {
			qui predictnl double `newvarname' = (1/`t' * (`dxb1')*exp(`xb1_plus') - 1/`t' * (`dxb0')*exp(`xb0_plus'))*`per' ///
												if `touse', `predictnl_opts' level(`level')
		}
		else if "`e(scale)'" =="odds" {
			qui predictnl double `newvarname' =  (1/`t' *(`dxb1')*exp(`xb1_plus')/((1 + exp(`xb1_plus'))) - ///
												1/`t' *(`dxb0')*exp(`xb0_plus')/((1 + exp(`xb0_plus'))))*`per' ///
												if `touse', `predictnl_opts' level(`level')
		}
		else if "`e(scale)'" =="normal" {
				qui predictnl double `newvarname' = (1/`t' *(`dxb1')*normalden(`xb1_plus')/normal(-(`xb1_plus')) - /// 
													1/`t' *(`dxb0')*normalden(`xb0_plus')/normal(-(`xb0_plus')))*`per' ///
													if `touse', `predictnl_opts' level(`level')
		}
		else if "`e(scale)'" =="theta" {
			qui predictnl double `newvarname' = (1/`t' *((`dxb1')*exp(`xb1_plus'))/((exp([ln_theta][_cons])*exp(`xb1_plus') + 1)) - ///
												1/`t' *((`dxb0')*exp(`xb0_plus'))/((exp([ln_theta][_cons])*exp(`xb0_plus') + 1)))*`per' ///
												if `touse', `predictnl_opts' level(`level')
		}
	}

/* Predict Difference in Survival Curves */
	else if "`sdiff1'" != "" {
		if `"`ci'"' != "" {
			local predictnl_opts "ci(`newvarname'_lci `newvarname'_uci)"
		}
		else if "`stdp'" != "" {
			local predictnl_opts se(`newvarname'_se)
		}

		forvalues i=1/`e(dfbase)' {
			local xb1_plus `xb1_plus' [xb][_rcs`i']*_rcs`i'
			local xb0_plus `xb0_plus' [xb][_rcs`i']*_rcs`i'
			if `i' != `e(dfbase)' {
				local xb1_plus `xb1_plus' +
				local xb0_plus `xb0_plus' +
			}
		}

/* use Parse_list to select appropriate values of factor variables */
		Parse_list, listname(sdiff1) parselist(`sdiff1')
		tokenize `r(retlist)'
		while "`1'"!="" {
			cap confirm var `2'
			if _rc {
				if "`2'" == "." {
					local 2 `1'
				}
				else {
					cap confirm num `2'
					if _rc {
						di as err "invalid sdiff1(... `1' `2' ...)"
						exit 198
					}
				}
			}
			if "`xb1_plus'" != "" & "`2'" != "0" & `: list posof `"`1'"' in main_varlist' != 0 {
				local xb1_plus `xb1_plus' +
			}
			if "`2'" != "0" & `: list posof `"`1'"' in main_varlist' != 0 {
				local xb1_plus `xb1_plus' [xb][`1']*`2' 
			}
			if `"`: list posof `"`1'"' in etvc'"' != "0" & "`2'" != "0" {
				local xb1_plus `xb1_plus' +

				forvalues i=1/`e(df_`1')' {
					local xb1_plus `xb1_plus' [xb][_rcs_`1'`i']*_rcs_`1'`i'*`2'  
					if `i' != `e(df_`1')' {
						local xb1_plus `xb1_plus' +
					}
				}
			}
			mac shift 2
		}			

		if "`sdiff2'" != "" {
/* use Parse_list to select appropriate values of factor variables */
		Parse_list, listname(sdiff2) parselist(`sdiff2')
		tokenize `r(retlist)'
			while "`1'"!="" {
				cap confirm var `2'
				if _rc {
					if "`2'" == "." {
						local 2 `1'
					}
					else {
						cap confirm num `2'
						if _rc {
							di as err "invalid sdiff2(... `1' `2' ...)"
							exit 198
						}
					}
				}
				if "`2'" != "0" & `: list posof `"`1'"' in main_varlist' != 0 {
					local xb0_plus `xb0_plus' + [xb][`1']*`2' 
				}
				if `"`: list posof `"`1'"' in etvc'"' != "0" & "`2'" != "0" {
					local xb0_plus `xb0_plus' + 
					forvalues i=1/`e(df_`1')' {
						local xb0_plus `xb0_plus' [xb][_rcs_`1'`i']*_rcs_`1'`i'*`2'
						if `i' != `e(df_`1')' {
							local xb0_plus `xb0_plus' +
						}
					}
				}
				mac shift 2
			}
		}
		if "`e(noconstant)'" == "" {
			local xb0_plus `xb0_plus' + [xb][_cons]
			local xb1_plus `xb1_plus' + [xb][_cons]
		}

		if "`e(scale)'" =="hazard" {
			qui predictnl double `newvarname' = exp(-exp(`xb1_plus')) - exp(-exp(`xb0_plus')) if `touse', `predictnl_opts' level(`level')
		}
		else if "`e(scale)'" =="odds" {
			qui predictnl double `newvarname' =  	1/(exp(`xb1_plus')+1) - 1/(exp(`xb0_plus')+1) if `touse', `predictnl_opts' level(`level')
		}
		else if "`e(scale)'" =="normal" {
			qui predictnl double `newvarname' =  	normal(-(`xb1_plus')) - normal(-(`xb0_plus')) if `touse', `predictnl_opts' level(`level')
		}
		else if "`e(scale)'" =="theta" {
			qui predictnl double `newvarname' =  	(exp([ln_theta][_cons])*exp(`xb1_plus') + 1)^(-1/exp([ln_theta][_cons])) ///
											-(exp([ln_theta][_cons])*exp(`xb0_plus') + 1)^(-1/exp([ln_theta][_cons])) ///
											if `touse', `predictnl_opts' level(`level')
		}
	}

/* estimate cure, survival of uncured or hazard of uncured */
	else if "`cure'" != "" | "`uncured'" != "" {
		local xblist [xb][_cons]
		local rcslist
		local drcslist
		tempvar temp
		if "`ci'" != "" {
			local prednlopt ci(`temp'_lci `temp'_uci)
		}
		foreach var in `e(varlist)' {
			local xblist `xblist' + [xb][`var']*`var'
		}
		
		if "`cure'" != "" {		/*if cure is specified this is what we want to estimate*/
			qui predictnl double `temp' = `xblist' if `touse', `prednlopt' level(`level')
			qui gen double `newvarname' = exp(-exp(`temp'))	if `touse'	/*we model on log(-log) scale*/
			if "`ci'" != "" {		
				qui gen double `newvarname'_lci = exp(-exp(`temp'_uci))	if `touse'
				qui gen double `newvarname'_uci = exp(-exp(`temp'_lci)) if `touse'
			}
		}
		
		else {		/*continue, estimate survival or hazard of uncured or Predicted survival time among uncured for a given centile*/
			forvalues i = 1/`e(dfbase)' {	
				if "`rcslist'" == "" local rcslist [xb][_rcs`i']*_rcs`i'			/*create a list of the sum of all spline variables*/
				else local rcslist `rcslist' + [xb][_rcs`i']*_rcs`i'
			}
			foreach var in `e(tvc)' {
				forvalues i = 1/`e(df_`var')' {
					local rcslist `rcslist' + [xb][_rcs_`var'`i']*_rcs_`var'`i'
				}
			}
			forvalues i = 1/`e(dfbase)' {
				if "`drcslist'" == "" local drcslist [xb][_rcs`i']*_d_rcs`i'		/*need derivatives of rcslist to calculate hazard and for centile*/
				else local drcslist `drcslist' + [xb][_rcs`i']*_d_rcs`i'
			}
			foreach var in `e(tvc)' {
				forvalues i = 1/`e(df_`var')' {
						local drcslist `drcslist' + [xb][_rcs_`var'`i']*_d_rcs_`var'`i'
				}
			}		
			local pi exp(-exp(`xblist')) 		/*we need cure for estimation of survival and hazard*/
			local exprcs exp(`rcslist') 		/*we need exp of the sum of all spline variables for estimation of survival and hazard*/
	/*predicted survival of uncured*/		
			if "`survival'" != "" & "`uncured'" != "" {			
				tokenize `e(boundary_knots)'
				local lastknot = `2' 
				qui predictnl double `temp' = ln(-(ln(`pi'^(`exprcs') - `pi') - ln(1 - `pi'))) if `touse', `prednlopt' level(`level')
				qui gen double `newvarname' = exp(-exp(`temp')) if `touse'
				qui replace `newvarname' = 0 if `newvarname' == . & `t'>=`lastknot' & `touse'
				if "`ci'" != "" {
					qui gen double `newvarname'_lci = exp(-exp(`temp'_uci)) if `touse'
					qui gen double `newvarname'_uci = exp(-exp(`temp'_lci)) if `touse'
				}
			}
	/*predicted hazard of uncured*/		
			else if "`hazard'" != "" & "`uncured'" != "" {	  
				qui predictnl double `temp' = ln(-ln(`pi')*((`drcslist')/`t')*`exprcs'*`pi'^(`exprcs'))- ln(`pi'^(`exprcs') - `pi')  if `touse', `prednlopt' level(`level') 
				qui gen double `newvarname' = exp(`temp') if `touse'
				if "`ci'" != "" {
					qui gen double `newvarname'_lci = exp(`temp'_lci) if `touse'
					qui gen double `newvarname'_uci = exp(`temp'_uci) if `touse'
				}
			}
	/* Predicted survival time among uncured for a given centile */
			else if "`centile'" != "" & "`uncured'" != ""{
				tempvar centilevar	
				gen `centilevar' = 1 - `centile'/100 if `touse'	
				/* initial values */	
				tempvar nr_time nr_time_old nr_surv nr_haz maxerr 
				qui gen double `nr_time' = `startt'
				qui gen double `nr_time_old' = `nr_time' if `touse'
				/* loop */
				local done 0
				local itercount 0
				while !`done' {
					local itercount = `itercount' + 1
					if `itercount'>=`centiter' {
						continue, break
					}
					qui predict `nr_surv' if `touse', s uncured timevar(`nr_time_old')
					qui predict `nr_haz' if `touse', h uncured timevar(`nr_time_old')
					qui replace `nr_time' = exp(ln(`nr_time_old') - (`nr_surv' - `centilevar')/(-`nr_haz'*`nr_surv'*`nr_time_old')) if `touse'
					qui gen double `maxerr' = abs(`nr_time' - `nr_time_old') if `touse'
					summ `maxerr' if `touse', meanonly
					if r(max)<`centol' {
						local done 1
					}
					else {
						drop `nr_surv' `nr_haz' `maxerr'
						qui replace `nr_time_old' = `nr_time' if `touse'
					}
				}

				qui gen double `newvarname' = `nr_time' if `touse'
				
				if "`ci'" != "" {
					tempvar ln_nr_time lnln_s lnln_s_se h tp_se s_tp
					drop _rcs*
					drop _d_rcs*
					qui gen double `ln_nr_time' = ln(`nr_time') if `touse'
					if "`e(rcsbaseoff)'" == "" {
						if "`e(orthog)'" != "" {
							tempname rmatrix
							matrix `rmatrix' = e(R_bh)
							local rmatrixopt rmatrix(`rmatrix')
						}
						qui rcsgen `ln_nr_time'  if `touse', knots(`e(ln_bhknots)') gen(_rcs) dgen(_d_rcs) `e(reverse)' `rmatrixopt'	
					}

					foreach tvcvar in `e(tvc)' {
						if "`e(orthog)'" != "" {
							tempname rmatrix_`tvcvar'
							matrix `rmatrix_`tvcvar'' = e(R_`tvcvar')
							local rmatrixopt rmatrix(`rmatrix_`tvcvar'')
						}
						if `e(df_`tvcvar')' == 1 {
							qui rcsgen `ln_nr_time' if `touse',  gen(_rcs_`tvcvar') dgen(_d_rcs_`tvcvar') `e(reverse)' `rmatrixopt'
						}
						else if `e(df_`tvcvar')' != 1 {
							qui rcsgen `ln_nr_time' if `touse', knots(`e(ln_tvcknots_`tvcvar')') gen(_rcs_`tvcvar') dgen(_d_rcs_`tvcvar') `e(reverse)' `rmatrixopt'
						}

					
						forvalues k = 1/`e(df_`tvcvar')' {
							qui replace _rcs_`tvcvar'`k' = _rcs_`tvcvar'`k' * `tvcvar' if `touse'
							qui replace _d_rcs_`tvcvar'`k' = _d_rcs_`tvcvar'`k' * `tvcvar' if `touse'
						}
					}	
					local rcslisttp
					forvalues i = 1/`e(dfbase)' {	
						if "`rcslisttp'" == "" local rcslisttp [xb][_rcs`i']*_rcs`i'			
						else local rcslisttp `rcslisttp' + [xb][_rcs`i']*_rcs`i'
					}
					foreach var in `e(tvc)' {
						forvalues i = 1/`e(df_`var')' {
							local rcslisttp `rcslisttp' + [xb][_rcs_`var'`i']*_rcs_`var'`i'
						}
					}
					
					local exprcstp exp(`rcslisttp') 		
					qui predictnl double `lnln_s' = ln(-(ln(`pi'^(`exprcstp') - `pi') - ln(1 - `pi'))) if `touse', se(`lnln_s_se')
					predict  `s_tp' if `touse', s uncured `offset' 
					qui replace `s_tp' = ln(`s_tp') if `touse'
					predict `h' if `touse', h uncured timevar(`nr_time')
					qui gen double `tp_se' = abs(`s_tp'*`lnln_s_se'/`h') if `touse'
					qui gen double `newvarname'_lci = `newvarname' - invnormal(1-0.5*(1-`level'/100))*`tp_se' if `touse'
					qui gen double `newvarname'_uci = `newvarname' + invnormal(1-0.5*(1-`level'/100))*`tp_se' if `touse'
				}	
			}		
		}		
	}
	
	
/* Predicted survival time for a given centile */
/* Estimated using Newton-Raphson alogorithm (updated from ridder method in version 1.2.2) */
	else if "`centile'" != "" & "`uncured'" == "" {
/* transform to appropriate scale */
		tempvar transcent centilevar
		
		gen `centilevar' = 1 - `centile'/100 if `touse'
		if "`e(scale)'" == "hazard" {
			qui gen double `transcent' = ln(-ln(`centilevar')) if `touse'
		}
		else if "`e(scale)'" == "odds" {
			qui gen double `transcent' = ln(1/`centilevar'-1) if `touse'
		}
		else if "`e(scale)'" == "normal" {
			qui gen double `transcent' = -invnorm(`centilevar') if `touse'
		}
		else if "`e(scale)'" == "theta" {
			qui gen double `transcent' = ln((`centilevar'^(-exp([ln_theta][_cons]))-1)/exp([ln_theta][_cons])) if `touse'
		}
		
/* initial values */	
		
		tempvar tmpxb nr_time nr_time_old nr_xb nr_dxb maxerr 
		qui gen double `nr_time' = `midt'
		qui gen double `nr_time_old' = `nr_time' if `touse'
		
/* loop */
		local done 0
		while !`done' {
			qui predict `nr_xb' if `touse', xb timevar(`nr_time_old') `offset'
			qui predict `nr_dxb' if `touse', dxb timevar(`nr_time_old')
			qui replace `nr_time' = exp(ln(`nr_time_old') - (`nr_xb' - `transcent')/`nr_dxb') if `touse'
			qui gen double `maxerr' = abs(`nr_time' - `nr_time_old') if `touse'
			summ `maxerr' if `touse', meanonly
			if r(max)<`centol' {
				local done 1
			}
			else {
				drop `nr_xb' `nr_dxb' `maxerr'
				qui replace `nr_time_old' = `nr_time' if `touse'
			}
		}
		
		qui gen double `newvarname' = `nr_time' if `touse'


		if "`ci'" != "" {
			tempvar lnln_s lnln_s_se h tp_se
			drop _rcs*
			drop _d_rcs*
			tempvar ln_nr_time
			qui gen double `ln_nr_time' = ln(`nr_time') if `touse'
			if "`e(rcsbaseoff)'" == "" {
				if "`e(orthog)'" != "" {
					tempname rmatrix
					matrix `rmatrix' = e(R_bh)
					local rmatrix rmatrix(`rmatrix')
				}
				qui rcsgen `ln_nr_time'  if `touse', knots(`e(ln_bhknots)') gen(_rcs) dgen(_d_rcs) `e(reverse)' `rmatrixopt'
			
				unab rcsterms :_rcs*		
				unab drcsterms :_d_rcs*		
			}

			foreach tvcvar in `e(tvc)' {
				if "`e(orthog)'" != "" {
					tempname rmatrix_`tvcvar'
					matrix `rmatrix_`tvcvar'' = e(R_`tvcvar')
					local rmatrixopt rmatrix(`rmatrix_`tvcvar'')
				}
				if `e(df_`tvcvar')' == 1 {
					qui rcsgen `ln_nr_time' if `touse',   gen(_rcs_`tvcvar') dgen(_d_rcs_`tvcvar') `e(reverse)' `rmatrixopt'
				}
				else if `e(df_`tvcvar')' != 1 {
					qui rcsgen `ln_nr_time' if `touse', knots(`e(ln_tvcknots_`tvcvar')') gen(_rcs_`tvcvar') dgen(_d_rcs_`tvcvar') `e(reverse)' `rmatrixopt' 
				}

				unab rcsterms_`tvcvar' : _rcs_`tvcvar'*
				unab drcsterms_`tvcvar' : _d_rcs_`tvcvar'*		

/*
				if "`e(orthog)'" != "" {		
					mata st_store(.,tokens(st_global("e(rcsterms_`tvcvar')")),"`touse'",(st_data(.,(tokens(st_global("e(rcsterms_`tvcvar')") + " `ocons'" )),"`touse'")*luinv(st_matrix("e(R_`tvcvar')")))[,1..`e(df_`tvcvar')'])							
					mata st_store(.,tokens(st_global("e(drcsterms_`tvcvar')")),"`touse'",(st_data(.,(tokens(st_global("e(drcsterms_`tvcvar')"))),"`touse'")*luinv(st_matrix("e(R_`tvcvar')"))[1..`e(df_`tvcvar')',1..`e(df_`tvcvar')'])[,1..`e(df_`tvcvar')'])
				}
*/				
				forvalues k = 1/`e(df_`tvcvar')' {
					qui replace _rcs_`tvcvar'`k' = _rcs_`tvcvar'`k' * `tvcvar' if `touse'
					qui replace _d_rcs_`tvcvar'`k' = _d_rcs_`tvcvar'`k' * `tvcvar' if `touse'
				}
			}	

			if "`e(scale)'" == "hazard" {
				qui predictnl double `lnln_s' = xb(xb)`addoff' if `touse', se(`lnln_s_se') 
			}
			else if "`e(scale)'" == "odds" {
				qui predictnl double `lnln_s' = ln(ln(1+exp(xb(xb)`addoff'))) if `touse', se(`lnln_s_se')
			}
			
			else if "`e(scale)'" == "normal" {
				qui predictnl double `lnln_s' = ln(-ln(normal(xb(xb)`addoff'))) if `touse', se(`lnln_s_se') 
			}
			
			else if "`e(scale)'" == "theta" {
				qui predictnl double `lnln_s' = ln(ln(exp(xb(ln_theta))*exp(xb(xb)`addoff')+1)/exp(xb(ln_theta))) if `touse', se(`lnln_s_se')
			}
			
			if "`e(scale)'" == "hazard" {
				qui predictnl double `h'=(1/(`newvarname'))*(xb(dxb))*exp(xb(xb)`addoff') if `touse' 
			}
			else if "`e(scale)'" == "odds" {
				qui predictnl double `h'=1/(`newvarname')*(xb(dxb))*exp(xb(xb)`addoff')/(1+exp(xb(xb)`addoff')) if `touse'
			}
			else if "`e(scale)'" == "normal" {
				qui predictnl double `h'=1/(`newvarname')*(xb(dxb))*normalden(xb(xb)`addoff')/(normal(-(xb(xb)`addoff'))) if `touse' 
			}
			else if "`e(scale)'" == "theta" {
				qui predictnl double `h'=1/(`newvarname')*(xb(dxb))*exp(xb(xb)`addoff')/(exp(xb(ln_theta))*exp(xb(xb)`addoff') + 1) if `touse' 
			}
			tempvar s_tp
			predict  `s_tp' if `touse', s `offset'
			qui replace `s_tp' = ln(`s_tp') if `touse'
			qui gen double `tp_se' = -`s_tp'*`lnln_s_se'/`h' if `touse'

			qui gen double `newvarname'_lci = `newvarname' - invnormal(1-0.5*(1-`level'/100))*`tp_se' if `touse'
			qui gen double `newvarname'_uci = `newvarname' + invnormal(1-0.5*(1-`level'/100))*`tp_se' if `touse'
		}	
	}
	
/* restore original data and merge in new variables */
	local keep `newvarname'
	if "`ci'" != "" { 
		local keep `keep' `newvarname'_lci `newvarname'_uci
	}
	else if "`stdp'" != "" {
		local keep `keep' `newvarname'_se 
	}
	keep `keep'
	qui save `"`newvars'"'
	restore
	merge 1:1 _n using `"`newvars'"', nogenerate noreport
end




/* meansurv added to stpm2_pred as sub program */
* 10March2009 - added averages for models on odds, probit or theta scales.

program Stmeancurve, sortpreserve 
	version 10.0
	syntax newvarname [if] [in],[TIMEvar(varname) AT(string) noOFFSET] 
	marksample touse, novarlist
	local newvarname `varlist'

	tempvar t lnt touse_time
		
	preserve
	
	/* use timevar option or _t */
	if "`timevar'" == "" {
		qui gen `t' = _t if `touse'
		qui gen double `lnt' = ln(_t) if `touse'
	}
	else {
		qui gen double `t' = `timevar' 
		qui gen double `lnt' = ln(`timevar') 
	}

	/* index which time units are selected */
	gen `touse_time' = `t' != .
	
	/* generate ocons for use when orthogonalising splines */
	tempvar ocons
	gen `ocons' = 1

	/* Calculate new spline terms */
	if "`timevar'" != "" & "`e(rcsbaseoff)'" == "" {
		drop _rcs* _d_rcs*
		qui rcsgen `lnt' if `touse_time', knots(`e(ln_bhknots)') gen(_rcs) `e(reverse)'
		if "`e(orthog)'" != "" {
			mata st_store(.,tokens(st_global("e(rcsterms_base)")),"`touse_time'",(st_data(.,(tokens(st_global("e(rcsterms_base)") + " `ocons'" )),"`touse_time'")*luinv(st_matrix("e(R_bh)")))[,1..`e(dfbase)'])							
		}
	}
	
	if "`e(tvc)'" != "" {
         capture drop _rcs_* 
    }
 
	foreach tvcvar in `e(tvc)' {
		qui rcsgen `lnt' if `touse_time',  gen(_rcs_`tvcvar') knots(`e(ln_tvcknots_`tvcvar')') `e(reverse)'
		if "`e(orthog)'" != "" {
			mata st_store(.,tokens(st_global("e(rcsterms_`tvcvar')")),"`touse_time'",(st_data(.,(tokens(st_global("e(rcsterms_`tvcvar')") + " `ocons'" )),"`touse_time'")*luinv(st_matrix("e(R_`tvcvar')")))[,1..`e(df_`tvcvar')'])							
		}
	}
	
	/* Out of sample predictions using at() */
	if "`at'" != "" {
		tokenize `at'
		while "`1'"!="" {
			unab 1: `1'
			cap confirm var `2'
			if _rc {
				cap confirm num `2'
				if _rc {
					di as err "invalid at(... `1' `2' ...)"
					exit 198
				}
			}
			qui replace `1' = `2' if `touse'
			mac shift 2
		}
	}

	
	foreach tvcvar in `e(tvc)' {
		local rcstvclist `rcstvclist' `e(rcsterms_`tvcvar')'
	}

	if "`e(scale)'" == "theta" {
		local theta = exp([ln_theta][_cons])
	}
	
	qui count if `touse_time'
	local Nt `r(N)'
	
	mata: msurvpop("`newvarname'","`touse'","`touse_time'","`rcstvclist'")

/* restore original data and merge in new variables */
	local keep `newvarname'
	if "`ci'" != "" { 
		local keep `keep' `newvarname'_lci `newvarname'_uci
	}
	else if "`stdp'" != "" {
		local keep `keep' `newvarname'_se 
	}
	keep `keep'
	tempfile newvars
	qui save `"`newvars'"'
	restore
	merge 1:1 _n using `"`newvars'"', nogenerate noreport
end

/* 
program Parse_list converts the options give in hrnum, hrdenom, hdiff1, hdiff2, sdiff1 and sdiff2 
to factor notation
*/
program define Parse_list, rclass
	syntax, listname(string) parselist(string)
	tokenize `parselist'
	while "`1'"!="" {
		fvunab tmpfv: `1'
		local 1 `tmpfv'
		_ms_parse_parts `1'
		if "`r(type)'"!="variable" {
			display as error "level indicators of factor" /*
							*/ " variables may not be individually set" /*
							*/ " with the `listname'() option; set one value" /*
							*/ " for the entire factor variable"
			exit 198
		}
		cap confirm var `2'
		if _rc {
			cap confirm num `2'
			if _rc {
				if "`2'" != "." {
					di as err "invalid `listname'(... `1' `2' ...)"
					exit 198
				}
			}
		}
		local _varlist `_varlist' `1'
		local `1'_value `2'
		mac shift 2
	}

	_ms_extract_varlist `e(varlist)', noomitted
	local varlist_omitted `r(varlist)' 

/* check if any tvc variables  not in varlist_omitted */
	local tvconly
	foreach tvcvar in `e(tvc)' {
		mata: st_local("addtvconly",strofreal(subinword(st_local("varlist_omitted"),st_local("tvcvar"),"")==st_local("varlist_omitted")))
		if `addtvconly' {
			local tvconly `tvconly' `tvcvar'
		}
	}
	foreach var in `varlist_omitted' `tvconly'{
		_ms_parse_parts `var'
		local vartype `r(type)'
		local intmult
		local partadded
		foreach parse_var in `_varlist' {
			/* check parse_var in model */
			/*
			_ms_extract_varlist  `parse_var'
			if "`r(varlist)'" == "" {
				display as error "`parse_var' is not included in the model"
				exit 198
			}
			*/
			
			
			* NOW SEE IF MODEL VARIABLE IS LISTED IN PARSE_VAR
*			mata: st_local("invar",strofreal(strpos("`var'","`parse_var'")>0))
			_ms_parse_parts `var'
			local invar 0
			if "`r(k_names)'" == "" {
				if "`r(name)'" == "`parse_var'" {
					local invar 1
				}
			}
			else {
				forvalues i = 1/`r(k_names)' {
					if "`r(name`i')'" == "`parse_var'" {
						local invar 1
					}
				}
			}
			if `invar' {
				if "`vartype'" == "variable" {
					local retlist `retlist' `var' ``parse_var'_value'
				}
				else if "`vartype'" == "factor" {
					if `r(level)' == ``parse_var'_value' {
						local retlist `retlist' `var' 1
					}
					else {
						_ms_extract_varlist ``parse_var'_value'.`parse_var'
					}
				}
				else if "`vartype'" == "interaction" {
					if strpos("`var'","`parse_var'") >0 {
							_ms_parse_parts `var'
							forvalues i = 1/`r(k_names)' {
								if "`r(name`i')'" == "`parse_var'" {
									if "`r(op`i')'" == "``parse_var'_value'" {
										local intmult `intmult'*1
										local partadded yes
									}
									else if "`r(op`i')'" == "c" {
										local intmult `intmult'*``parse_var'_value'
										local partadded yes
									}
									else {
										local intmult `intmult'*0
									}
								}
							}
					}
					else {
						local intmult `intmult'*0
					}
				}
				else if "`vartype'" == "product" {
						display "products not currently available"
				}
			}
		}
		if "`vartype'" == "interaction" {
			local intmult = 1`intmult'
			if `intmult' != 0 {
				local retlist `retlist' `var' `intmult'
			}
		}
		return local retlist `retlist'
	}
end

program define _rmst, sortpreserve
/*
	Mean or restricted mean or restriced SD of time to failure for an stpm2 model,
	based on covariate patterns or an `at()' specification.

*/

version 11.1
syntax [if] [in], [ at(string) GENerate(string) N(int 1000) ///
		POWer(real 1) RSDst TMAx(real 0) TMIn(real 0) ZEROs ]
local cmdline `e(cmdline)'
if (`tmin' < 0) local tmin 0
if (`tmax' <= 0) local tmax .

// Use original variable names saved by stpm2
local varlist `e(varnames)'

quietly {
	marksample touse
	replace `touse' = 0 if _st == 0
	tempvar t s
	if missing(`tmax') {
		noi di as txt "[tmax() not specified, calculating mean survival to t = +infinity]"
		// Find t corresponding to 99.999999 centile for range of t
		local cmax 99.999999
		predict `t', centile(`cmax')
		sum `t', meanonly
		local tmax = r(max)
		drop `t'
	}
	local n_orig = _N
	if `n' > `n_orig' {
		tempvar mark_orig
		gen byte `mark_orig' = 1
	}
	range `t' `tmin' `tmax' `n'
	tempvar integral
	/* `power' option is undocumented and not implemented for rsdst */
	if `power' != 1 {
		tempvar lnt
		gen `lnt' = cond(`power' == 0, ln(`t') - ln(`tmax'),  ///
		 (`t'^`power' - 1) / `power' - (`tmax'^`power' - 1) / `power')
	}
	gen `integral' = .
	
/* identify covariate patterns and compute mean survival for them */
	if "`zeros'" == "" {
		if "`at'" != "" {
			// Identify covariates `s(omitted)' in `varlist' but not in `at'.
			CheckAt "`at'" "`varlist'"
			local covpatvarlist `s(omitted)'
		}
		else local covpatvarlist `varlist'
	}
	if "`covpatvarlist'" != "" {
		tempvar covpat
		covpat `covpatvarlist' if `touse', generate(`covpat')
		sum `covpat', meanonly
		local ncovpat = r(max)
		noi di as txt _n "Processing " as res `ncovpat' as txt " covariate patterns ..."
	}

	else {
		local ncovpat 1
		local covpat 1
	}

	forvalues i = 1 / `ncovpat' {
		local At
		tokenize `covpatvarlist'
		while "`1'" != "" {
			sum `1' if (`covpat' == `i') & (`touse' == 1), meanonly
			local mm = r(mean)
			local At `At' `1' `mm'
			mac shift
		}

		predict double `s', timevar(`t') at(`At' `at') survival `zeros'

		if `power' != 1 {
			// work with f(t) - f(tmax), so that on power-transformed scale, tmax (effective) = 0.
			replace `s' = 1 - `s' // distribution function
			integ `s' `lnt' if !missing(`lnt')
			replace `integral' = cond(`power' == 0, ln(`tmax'), (`tmax'^`power'-1)/`power') - r(integral) if `covpat' == `i'
		}
		else {
			replace `s' = cond(`t' == 0, 0, 1 - `s')
			integ `s' `t'

			local ex = `tmax' - r(integral)
			if "`rsdst'" != "" {
				replace `s' = `s' * `t'
				integ `s' `t'
				local ex2 = `tmax'^2 - 2 * r(integral)
				replace `integral' = sqrt(`ex2' - `ex'^2) if `covpat' == `i'
			}
			else replace `integral' = `ex' if `covpat' == `i'
		}		
		drop `s'
		if mod(`i', 10) == 0 noi di as txt `i', _cont
	}

	if !missing("`generate'") {
		cap confirm var `generate', exact
		if (c(rc) == 0) replace `generate' = `integral'
		else rename `integral' `generate'
/*
		local txt = cond("`rsdst'" != "", "SD", "mean")
		if (missing("`cmax'")) lab var `generate' "restricted `txt' at time `tmax'"
		else lab var `generate' "`txt' time"
*/
	}

	if `n' > `n_orig' {
		drop if missing(`mark_orig')
		drop `mark_orig'
	}
}
end

* version 1.0.1 PR 26Sep2002
* Based on PR 9-Jan-94. Based on covariate-pattern counter in lpredict.ado.
program define covpat, sortpreserve
version 11
syntax [varlist(default=none)] [if] [in], Generate(string)
confirm new var `generate'
tempvar keep
quietly {
	marksample keep
	sort `keep' `varlist'
	gen long `generate' = .
	by `keep' `varlist': replace `generate'=cond(_n==1 & `keep',1,.)
	replace `generate' = sum(`generate')
	replace `generate'=. if `generate'==0
	label var `generate' "covariate pattern"
}
end


program define CheckAt, sclass
version 11
/*
	Checks whether the `at' list specifies all variables in varlist.
	If not, provides in s(omitted) a list of omitted variables
*/
args at varlist
sret clear
local omitted `varlist'
tokenize `at'
while "`1'"!="" {
	fvunab v: `1'
	local 1 `v'
	_ms_parse_parts `1'
	if "`r(type)'"!="variable" {
		display as error "level indicators of factor" /*
						*/ " variables may not be individually set" /*
						*/ " with the at() option; set one value" /*
						*/ " for the entire factor variable"
		exit 198
	}
	cap confirm var `2'
	if _rc {
		cap confirm num `2'
		if _rc {
			di as err "invalid at(... `1' `2' ...)"
			exit 198
		}
	}
	ChkIn `1' "`varlist'"
	local omitted : list omitted - 1
	macro shift 2
}
sreturn local omitted `omitted'
end


program define ChkIn, sclass
version 9.2
* Returns s(k) = index # of target variable v in varlist, or 0 if not found.
args v varlist
sret clear
local k: list posof "`v'" in varlist
sret local k `k'
if `s(k)' == 0 {
   	di as err "`v' is not a valid covariate"
   	exit 198
}
end


/* mata program to obtain mean survival */
mata:
void msurvpop(string scalar newvar, string scalar touse, string scalar touse_time, | string scalar tvcrcslist) 
{
/* Transfer data from Stata */
	Nt = strtoreal(st_local("Nt"))
	hascons = st_global("e(noconstant)") == ""	

	if (st_global("e(rcsbaseoff)") == "") {
		rcsbase = st_data( ., tokens(st_global("e(rcsterms_base)")), touse_time)
	}
	else {
		rcsbase = J(Nt,0,0)
	}

	if (tvcrcslist != "") {
		rcstvc = st_data( ., tokens(tvcrcslist), touse_time)
	}
	else {
		rcstvc = I(0)
	}
	
	if (st_global("e(varlist)") != "") {
		x = st_data(.,tokens(st_global("e(varlist)")),touse)
	}
	else {
		x = J(1,0,.)	
	}

	tvcvar = tokens(st_global("e(tvc)"))
	ntvc = cols(tvcvar)

/* Beta matrix */		
	beta = st_matrix("e(b)")'[1..cols(rcsbase)+cols(rcstvc) + cols(tokens(st_global("e(varlist)"))) + hascons,1]

	scale = st_global("e(scale)")
	if (scale == "theta") theta = strtoreal(st_local("theta"))

/*	check whether to include offset */
	offset_name = st_global("e(offset1)")
	if (offset_name != "" & st_local("offset") != "nooffset") offset = st_data(.,offset_name,touse)
	else offset = 0


	startstop = J(ntvc,2,.)
	tvcpos = J(ntvc,1,.)
	tmpstart = 1

/* Loop over number of time observations */
	for (i=1;i<=ntvc;i++){
		startstop[i,1] = tmpstart
		tmpntvc = cols(tokens(st_global("e(rcsterms_"+tvcvar[1,i]+")")))
		startstop[i,2] = tmpstart + tmpntvc - 1
		tmpstart = startstop[i,2] + 1

		pos = .	
		maxindex(tvcvar[1,i]:==tokens(st_global("e(varlist)"))	,1,pos,.)
		tvcpos[i,1] = pos
	}

	
	Nx = rows(x)

/* Loop over all selected observations */	
	meansurv = J(Nt,1,0)

	if (hascons) {
		addcons = J(Nt,1,1)
	}
	else {
		addcons = J(Nt,0,0)
	}

	for (i=1;i<=Nx;i++) {
		tmprcs = J(Nt,0,.)
		for (j=1;j<=ntvc;j++) {
			tmprcs = tmprcs, rcstvc[,startstop[j,1]..startstop[j,2]]:*	x[i,tvcpos[j,1]]
		}
		if (scale == "hazard")	{
			meansurv = meansurv + exp(-exp((J(Nt,1,x[i,]),rcsbase, tmprcs,addcons)*beta :+ offset)):/Nx
		}
		else if (scale == "odds")	meansurv = meansurv + ((1 :+ exp((J(Nt,1,x[i,]),rcsbase, tmprcs,addcons)*beta :+ offset)):^(-1)):/Nx
		else if (scale == "normal")	meansurv = meansurv + normal(-(J(Nt,1,x[i,]),rcsbase, tmprcs,addcons)*beta :+ offset):/Nx
		else if (scale == "theta")	meansurv = meansurv + ((theta:*exp((J(Nt,1,x[i,]),rcsbase, tmprcs,addcons)*beta :+ offset) :+ 1):^(-1/theta)):/Nx
	}

	(void) st_addvar("double",newvar)	
	st_store(., newvar, touse_time,meansurv)
}	
end

	
	
