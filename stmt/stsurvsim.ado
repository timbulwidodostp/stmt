*! v 1.0.3 PR 04may2014
program define stsurvsim
// Simulate sample from Royston-Parmar distribution
syntax newvarname [if] [in], [ BEta(string) BKnots(numlist ascending) centol(real 0.001) ///
 Knots(numlist) KNScale(string) MAXiter(int 100) RMATrix(name) Scale(string) ]

if "`knscale'" != "" & "`knscale'" != "log" {
	di as err "invalid knscale(`knscale')"
	exit 198
}

tempname eb R
if "`beta'`bknots'`knots'`rmatrix'`scale'" == "" { // using stpm2 parameters
	// Extract beta from e(b)
	cap matrix `eb' = e(b)
	if c(rc) {
		di as err "coefficient vector not found - no stpm2 model was fitted"
		exit 198
	}
	if "`e(tvc)'" != "" {
		di as err "simulation from models with time-dependent effects not supported"
		exit 198
	}

	if ("`knscale'" == "log") di as txt "[knscale(log) ignored]"

	local names : colnames `eb'
	local i 1
	tokenize `names'
	while "``i''" != "" {
		if substr("``i''", 1, 6) == "_d_rcs" {
			local dcolnames `dcolnames' ``i''
			local dcoleq `dcoleq' dxb
			local bi = [dxb]_b[``i'']
			local b `b' `bi',
		}
		else {
			local coleq `coleq' xb
			local colnames `colnames' ``i''
			local bi = [xb]_b[``i'']
			local b `b' `bi',
		}
		local ++i
	}
	local b = substr("`b'", 1, length("`b'") - 1)

	// Extract knots
	local ln_bhknots `e(ln_bhknots)'
	local bk1 : word 1 of `e(boundary_knots)'
	local bk2 : word 2 of `e(boundary_knots)'
	local midt = (`bk1' + `bk2') / 2

	// Extract R-matrix
	cap confirm matrix e(R_bh)
	if c(rc) == 0 {
		matrix `R' = e(R_bh)
		local rmatrix rmatrix(`R')
	}

	// Extract scale
	local scale `e(scale)'
}
else { // user input values
	if "`scale'" == "" {
		di as err "scale() required"
		exit 198
	}
	if "`scale'" == "h" ///
	 | "`scale'" == "ha" ///
	 | "`scale'" == "haz" ///
	 | "`scale'" == "haza" ///
	 | "`scale'" == "hazar" ///
	 | "`scale'" == "hazard" {
	 	local scale hazard
	 }
	else if "`scale'" == "o" ///
	 | "`scale'" == "od" ///
	 | "`scale'" == "odd" ///
	 | "`scale'" == "odds" {
	 	local scale odds
	 }
	else if "`scale'" == "n" ///
	 | "`scale'" == "no" ///
	 | "`scale'" == "nor" ///
	 | "`scale'" == "norm" ///
	 | "`scale'" == "norma" ///
	 | "`scale'" == "normal" {
	 	local scale normal
	}
	if ("`scale'" != "hazard") & ("`scale'" != "odds") & ("`scale'" != "normal") {
		di as err "scale(`scale') invalid"
		exit 198
	}
	if "`bknots'" != "" & wordcount("`bknots'") != 2 {
		di as err "you must supply exactly 2 boundary knots"
		exit 198
	}
	if "`rmatrix'" != "" {
		confirm matrix `rmatrix'
		local rmatrix rmatrix(`rmatrix')
	}

	// Assemble string comprising log of (left boundary knot, interior knot(s), right boundary knot)
	local nk = wordcount("`knots'")
	local df = `nk' + 1
	local terms = `df' + 1
	if `nk' > 0 {
		local bk1 : word 1 of `bknots'
		if "`knscale'" == "log" {
			local ln_bhknots `bk1'
			local bk1 = exp(`bk1')
		}
		else local ln_bhknots = ln(`bk1')
		forvalues i = 1 / `nk' {
			local lnk : word `i' of `knots'
			if ("`knscale'" != "log") local lnk = ln(`lnk')
			local ln_bhknots `ln_bhknots' `lnk'
		}
		local bk2 : word 2 of `bknots'
		if "`knscale'" == "log" {
			local lnbk `bk2'
			local bk2 = exp(`bk2')
		}
		else local lnbk = ln(`bk2')
		local ln_bhknots `ln_bhknots' `lnbk'
		local midt = (`bk1' + `bk2') / 2
	}
	else {
		sum _t if _d==1
		local midt = (r(min) + r(max)) / 2
	}
	// Construct e(b) matrix
	// Extract covariates from `beta'
	local colnames
	local coleq
	local b
	tokenize `beta', parse(" =")
	local beta
	while "`1'" != "" {
		local vn `1'
		if ("`2'" == "=") macro shift
		confirm num `2'
		if substr("`vn'", 1, 4) == "_rcs" | substr("`vn'", 1, 6) == "_d_rcs" | "`vn'" == "_cons" {
			// Spline var begins with _rcs or _d_rcs
			local beta `beta' `2'
		}
		else {
			// Covariates must exist in the data. Splines don't need to.
			confirm var `vn'
			unab vn : `vn'
			local colnames `colnames' `vn'
			local coleq `coleq' xb
			local b `b' `2',
		}
		macro shift 2
	}

	forvalues i = 1 / `terms' {
		if `i' < `terms' {
			local colnames `colnames' _rcs`i'
			local dcolnames `dcolnames' _d_rcs`i'
			local dcoleq `dcoleq' dxb
		}
		else local colnames `colnames' _cons
		local coleq `coleq' xb
		local bi : word `i' of `beta'
		local b `b' `bi',
	}
	forvalues i = 1 / `df' {
		local bi : word `i' of `beta'
		local b `b' `bi'
		if (`i' < `df') local b `b',
	}
}

// Post vector of regression coefficients
matrix `eb' = `b'
matrix colnames `eb' = `colnames' `dcolnames'
matrix coleq `eb' = `coleq' `dcoleq'

// Predict survival time for given centiles. Estimated using Newton-Raphson algorithm.
// Transform survival function to appropriate scale
marksample touse, novarlist
tempvar transcent
if "`scale'" == "hazard" {
	qui gen double `transcent' = ln(-ln(runiform())) if `touse'
}
else if "`scale'" == "odds" {
	qui gen double `transcent' = ln(1/runiform()-1) if `touse'
}
else if "`scale'" == "normal" {
	qui gen double `transcent' = -invnormal(runiform()) if `touse'
}

// initial values	
tempvar nr_time ln_nr_time_old nr_time_old nr_xb nr_dxb maxerr 
qui gen double `nr_time' = `midt' if `touse'
qui gen double `nr_time_old' = `midt'

// Post beta
tempname esthold
capture _estimates hold `esthold'
local has_est = (c(rc) == 0)
ereturn post `eb'

// loop
local done 0
local i 0
quietly while !`done' {
	// Create spline basis functions
	cap drop _rcs* _d_rcs*
	gen double `ln_nr_time_old' = ln(`nr_time_old')
	rcsgen `ln_nr_time_old' if `touse', knots(`ln_bhknots') gen(_rcs) dgen(_d_rcs) `rmatrix'
	_predict `nr_xb' if `touse', xb equation(xb)
	_predict `nr_dxb' if `touse', xb equation(dxb)
	replace `nr_time' = exp(`ln_nr_time_old' - (`nr_xb' - `transcent')/`nr_dxb') if `touse'
//	gen double `maxerr' = abs(`nr_time' - `nr_time_old') if `touse' // !! change to reldif convergence, 02may2014
	gen double `maxerr' = reldif(`nr_time', `nr_time_old') if `touse'
	sum `maxerr' if `touse', meanonly
	if r(max) < `centol' {
		local done 1
	}
	else {
		drop `nr_xb' `nr_dxb' `maxerr' `ln_nr_time_old'
		qui replace `nr_time_old' = `nr_time' if `touse'
	}
	local ++i
	if (`i' > `maxiter') {
		di as err "convergence not achieved"
		local done 1
	}
}
di as txt "[" `i' " iterations performed]"
qui gen double `varlist' = `nr_time' if `touse'
lab var `varlist' "survival time simulated from RP `scale' model"
if `has_est' _estimates unhold `esthold'
end

program define rcsgen, rclass
* version 1.5.4 13Oct2013
	version 10.0
	syntax  [varlist(default=none)] [if] [in] ///
		, [Gen(string) DGen(string) Knots(numlist) BKnots(numlist max=2) Orthog Percentiles(numlist ascending) RMATrix(name) ///
			DF(int 0)  IF2(string) FW(varname)  REVerse SCAlar(string) NOSecondder NOFirstder]      

	marksample touse
	
// sort knots
	if "`knots'" != "" {
		numlist "`knots'", sort
		local knots `r(numlist)'
	}
    
/* Error checks */
	if "`scalar'"  != "" {
		if "`varlist'" != "" {
			display as error "You can't specify both a varname and the scalar option"
			exit 198
		}
		if "`df'" != "0" {
			display as error "You can't specify the df option with the scalar option"
			exit 198
		} 
		if "`percentiles'" != "" {
			display as error "You can't specify the percentiles option with the scalar option"
			exit 198
		} 
		if "`orthog'" != "" {
			display as error "You can't specify the orthog option with the scalar option"
			exit 198
		}
		if "`fw'" != "" {
			display as error "You can't specify the fw option with the scalar option"
			exit 198
		}
	}

	if "`knots'" != "" & "`percentiles'" != "" {
		display as err "Only one of the knots, df and percentiles options can be used"
		exit 198
	}
        
	if "`knots'" != "" & "`df'" != "0" {
		display as err "Only one of the knots, df and percentiles options can be used"
		exit 198
	}
	
	if "`df'" != "0" & "`percentiles'" != "" {
		display as err "Only one of the knots, df and percentiles options can be used"
		exit 198
	} 

	if "`bknots'" != "" & "`df'" == "0" {
		display as err "Boundary knots can only be defined with the degrees of freedom option"
		exit 198
	} 

	if "`orthog'" != "" & "`rmatrix'" != "" {
		display as error "Only one of the orthog and rmatrix options  can be specified"
		exit 198
	}
	
	if "`gen'" == "" {
		di in red "Must specify name for cubic splines basis"
		exit 198
	}
    
/* percentiles option */             
	if "`percentiles'" != "" {
		if "`fw'" != "" {
			local fw [fw=`fw']
		}
		if "`if2'" != "" {
			local aif & `if2'
		}
		local knots

		foreach ptile in `percentiles' {
			summ `varlist' if `touse' `aif', meanonly
			if `ptile' == 0 {
				local knots `r(min)'
			}
			else if `ptile' == 100 {
				local knots `knots' `r(max)'
			}
			else {
				_pctile `varlist' if `touse' `aif' `fw', p(`ptile')
				local knots `knots' `r(r1)'
			}
		}
	}

/* Find knot locations if df option is used */
	if "`df'" > "1" {
		if "`fw'" != "" {
			local fw [fw=`fw']
		}
		if "`if2'" != "" {
			local aif & `if2'
		}
		if "`bknots'"!="" {
			local lowerknot: word 1 of `bknots'
			local upperknot: word 2 of `bknots'
		}
		else {
			quietly summ `varlist' if `touse' `aif', meanonly
			local lowerknot `r(min)'
			local upperknot `r(max)'
		}
		local dfm1=`df'-1        

		forvalues y= 1/`dfm1' {
			local centile=(100/`df')*`y'
			local centilelist `centilelist' `centile'
		}

		local intknots
		foreach ctile in `centilelist' {
			_pctile `varlist' if `touse' `aif' `fw', p(`ctile')
			local intknots `intknots' `r(r1)'
		}
		if real(word("`intknots'",1))<=`lowerknot' {
			display as err "Lowest internal knot is not greater than lower boundary knot"
			exit 198
		}
		if real(word("`intknots'",`dfm1'))>=`upperknot' {
			display as err "Highest internal knot is not greater than upper boundary knot"
			exit 198
		}
	 	local knots
		local knots  `lowerknot' `intknots' `upperknot'
	}

/*Derive the spline variables in the default way (not backwards)*/
		
	if "`reverse'" == "" & "`nosecondder'"  == "" & "`nofirstder'" == "" {
	/* Start to derive spline variables */
		if "`scalar'" == "" {
			quietly gen double `gen'1 = `varlist' if `touse'
		}
		else {
			scalar `gen'1 = `scalar'
		}			

	/* generate first derivative if dgen option is specified */
		if "`dgen'" != "" {
			if "`scalar'" == "" {
				quietly gen double `dgen'1 = 1 if `touse'
			}
			else {
				scalar `dgen'1 = 1
			}			
		}
		local rcslist `gen'1 
		local drcslist `dgen'1
	
		local nk : word count `knots'
		if "`knots'" == "" {
			local interior  = 0
		}
		else {
			local interior  = `nk' - 2
		}
		local nparams = `interior' + 1
	
		if "`knots'" != "" {
			local i = 1 
			tokenize "`knots'"
			while "``i''" != "" {
				local k`i' ``i''
				local i = `i' + 1
			}
	
			local kmin = `k1'
			local kmax = `k`nk''
	
			forvalues j=2/`nparams' {
				local lambda = (`kmax' - `k`j'')/(`kmax' - `kmin')
				if "`scalar'" == "" {
					quietly gen double `gen'`j' = ((`varlist'-`k`j'')^3)*(`varlist'>`k`j'') - ///
								`lambda'*((`varlist'-`kmin')^3)*(`varlist'>`kmin') - ///
								(1-`lambda')*((`varlist'-`kmax')^3)*(`varlist'>`kmax')  if `touse'
				}
				else {
					scalar `gen'`j' = ((`scalar'-`k`j'')^3)*(`scalar'>`k`j'') - ///
								`lambda'*((`scalar'-`kmin')^3)*(`scalar'>`kmin') - ///
								(1-`lambda')*((`scalar'-`kmax')^3)*(`scalar'>`kmax')  

				}
				
				local rcslist `rcslist' `gen'`j'
	
	/* calculate derivatives */
				if "`dgen'"!="" {
					if "`scalar'" == "" {
						quietly gen double `dgen'`j' = (3*(`varlist'- `k`j'')^2)*(`varlist'>`k`j'') - ///
									`lambda'*(3*(`varlist'-`kmin')^2)*(`varlist'>`kmin') - ///
									(1-`lambda')*(3*(`varlist'-`kmax')^2)*(`varlist'>`kmax') 
					}
					else {
						scalar `dgen'`j' = (3*(`scalar'- `k`j'')^2)*(`scalar'>`k`j'') - ///
									`lambda'*(3*(`scalar'-`kmin')^2)*(`scalar'>`kmin') - ///
									(1-`lambda')*(3*(`scalar'-`kmax')^2)*(`scalar'>`kmax') 
					}
					local drcslist `drcslist' `dgen'`j'
				}       
			}
		}
	}

/*Derive the spline variables in reversed order */		/*ADDED: 2010-02-02 by Therese Andersson*/
		
	else if "`reverse'" != "" & "`nosecondder'" == "" & "`nofirstder'" == "" {
		local rcslist  
		local drcslist

		local nk : word count `knots'
		if "`knots'" == "" {
			local interior  = 0
		}
		else {
			local interior  = `nk' - 2
		}
		local nparams = `interior' + 1

		if "`knots'" != "" {
			local i = 1 
			tokenize "`knots'"
			while "``i''" != "" {
				local k`i' ``i''
				local i = `i' + 1
			}

			local kmin = `k1'
			local kmax = `k`nk''

			forvalues j=1/`interior' {
				local h = `nk'-`j'
				local lambda = (`k`h''-`kmin')/(`kmax' - `kmin')
				if "`scalar'" == "" {
					quietly gen double `gen'`j' = ((`k`h''-`varlist')^3)*(`k`h''>`varlist') - ///
								`lambda'*((`kmax'-`varlist')^3)*(`kmax'>`varlist') - ///
								(1-`lambda')*((`kmin'-`varlist')^3)*(`kmin'>`varlist')  if `touse'
				}
				else {
					scalar `gen'`j' = ((`k`h''-`scalar')^3)*(`k`h''>`scalar') - ///
								`lambda'*((`kmax'-`scalar')^3)*(`kmax'>`scalar') - ///
								(1-`lambda')*((`kmin'-`scalar')^3)*(`kmin'>`scalar') 
				}
				local rcslist `rcslist' `gen'`j'

/* calculate derivatives */
				if "`dgen'"!="" {
					if "`scalar'" == "" {
						quietly gen double `dgen'`j' = (-3*(`k`h''-`varlist')^2)*(`k`h''>`varlist') - ///
									`lambda'*(-3*(`kmax'-`varlist')^2)*(`kmax'>`varlist') - ///
									(1-`lambda')*(-3*(`kmin'-`varlist')^2)*(`kmin'>`varlist')  if `touse'
					}
					else {
						scalar `dgen'`j' = (-3*(`k`h''-`scalar')^2)*(`k`h''>`scalar') - ///
									`lambda'*(-3*(`kmax'-`scalar')^2)*(`kmax'>`scalar') - ///
									(1-`lambda')*(-3*(`kmin'-`scalar')^2)*(`kmin'>`scalar') 

					}
					local drcslist `drcslist' `dgen'`j'
				}       
			}
/* Derive last spline variable */
			if "`scalar'" == "" {
				quietly gen double `gen'`nparams' = `varlist' if `touse'
			}
			else {
				scalar `gen'`nparams' = `scalar' 
			}
			local rcslist `rcslist' `gen'`nparams'

/* generate first derivative if dgen option is specified */
			if "`dgen'" != "" {
				if "`scalar'" == "" {
					quietly gen double `dgen'`nparams' = 1 if `touse'
				}
				else {
					scalar `dgen'`nparams' = 1 
				}
				local drcslist `drcslist' `dgen'`nparams'
			}
		}
	}

//  no second derivative
		
	else if "`nosecondder'" != "" & "`reverse'" == "" & "`nofirstder'" == "" {
  				
				/* Start to derive spline variables */
                if "`scalar'" == "" {
                        quietly gen double `gen'1 = `varlist' if `touse'
                }
                else {
                        scalar `gen'1 = `scalar'
                }                       

        /* generate first derivative if dgen option is specified */
                if "`dgen'" != "" {
                        if "`scalar'" == "" {
                                quietly gen double `dgen'1 = 1 if `touse'
                        }
                        else {
                                scalar `dgen'1 = 1
                        }                       
                }
                local rcslist `gen'1 
                local drcslist `dgen'1

                local nk : word count `knots'
                if "`knots'" == "" {
                        local interior  = 0
                }
                else {
                        local interior  = `nk' - 2
                }
                local nparams = `interior' + 2
				local npar = `interior' + 1

                if "`knots'" != "" {
                        local i = 1 
                        tokenize "`knots'"
                        while "``i''" != "" {
                                local k`i' ``i''
                                local i = `i' + 1
                        }

                        local kmin = `k1'
                        local kmax = `k`nk''

                        forvalues j=2/`npar' {
                                local lambda = (`kmax' - `k`j'')/(`kmax' - `kmin')
                                if "`scalar'" == "" {
                                        quietly gen double `gen'`j' = ((`varlist'-`k`j'')^3)*(`varlist'>`k`j'') - ///
                                                                `lambda'*((`varlist'-`kmin')^3)*(`varlist'>`kmin') - ///
                                                                (1-`lambda')*((`varlist'-`kmax')^3)*(`varlist'>`kmax') if `touse'
                                }
                                else {
                                        scalar `gen'`j' = ((`scalar'-`k`j'')^3)*(`scalar'>`k`j'') - ///
                                                                `lambda'*((`scalar'-`kmin')^3)*(`scalar'>`kmin') - ///
                                                                (1-`lambda')*((`varlist'-`kmax')^3)*(`scalar'>`kmax')  

                                }
                                
                                local rcslist `rcslist' `gen'`j'
        
        /* calculate derivatives */
                                if "`dgen'"!="" {
                                        if "`scalar'" == "" {
                                                quietly gen double `dgen'`j' = (3*(`varlist'- `k`j'')^2)*(`varlist'>`k`j'') - ///
                                                                        `lambda'*(3*(`varlist'-`kmin')^2)*(`varlist'>`kmin') - ///
                                                                        (1-`lambda')*(3*(`varlist'-`kmax')^2)*(`varlist'>`kmax') 
                                        }
                                        else {
                                                scalar `dgen'`j' = (3*(`scalar'- `k`j'')^2)*(`scalar'>`k`j'') - ///
                                                                        `lambda'*(3*(`scalar'-`kmin')^2)*(`scalar'>`kmin') - ///
                                                                        (1-`lambda')*(3*(`scalar'-`kmax')^2)*(`scalar'> `kmax') 
                                        }
                                        local drcslist `drcslist' `dgen'`j'
                                }       
                        }
				
				/* Derive last spline variable */
						local c=(1/(3*(`kmax' - `kmin')))						
                        if "`scalar'" == "" {
                                quietly gen double `gen'`nparams' = (`varlist'-`kmin')^2*(`varlist'>`kmin') - ///
								`c'*((`varlist'-`kmin')^3)*(`varlist'>`kmin') + ///
                                                                 `c'*((`varlist'-`kmax')^3)*(`varlist'>`kmax') if `touse'
                        }
                        else {
                                scalar `gen'`nparams' = (`scalar'-`kmin')^2*(`scalar'>`kmin') - ///
								`c'*((`scalar'-`kmin')^3)*(`scalar'>`kmin') + ///
                                                                 `c'*((`scalar'-`kmax')^3)*(`scalar'>`kmax') 
                        }
                        local rcslist `rcslist' `gen'`nparams'

/* generate first derivative if dgen option is specified */
                        if "`dgen'" != "" {
                                if "`scalar'" == "" {
                                        quietly gen double `dgen'`nparams' = 2*(`varlist'-`kmin')*(`varlist'>`kmin') - ///
								 3*`c'*((`varlist'-`kmin')^2)*(`varlist'>`kmin') + ///
                                                                 3*`c'*((`varlist'-`kmax')^2)*(`varlist'>`kmax') if `touse' 
                                }
                                else {
                                        scalar `dgen'`nparams' = 2*(`scalar'-`kmin')*(`scalar'>`kmin') - ///
								 3*`c'*((`scalar'-`kmin')^2)*(`scalar'>`kmin') + ///
                                                                 3*`c'*((`scalar'-`kmax')^2)*(`scalar'>`kmax') 
                                }
                                local drcslist `drcslist' `dgen'`nparams'
                        }
                }
        }
 
 ****************************!!!*********************			
		
	else if "`nosecondder'" != "" & "`reverse'" != "" & "`nofirstder'" == ""{
                			
                local nk : word count `knots'
                if "`knots'" == "" {
                        local interior  = 0
                }
                else {
                        local interior  = `nk' - 2
                }
                local nparams = `interior' + 2
		local npar = `interior' + 1

                if "`knots'" != "" {
                        local i = 1 
                        tokenize "`knots'"
                        while "``i''" != "" {
                                local k`i' ``i''
                                local i = `i' + 1
                        }

                        local kmin = `k1'
                        local kmax = `k`nk''


	/* Derive first spline variable */
                        local c=(1/(3*(`kmax' - `kmin')))						
                        if "`scalar'" == "" {
                                quietly gen double `gen'1 = (`kmax'-`varlist')^2*(`kmax'>`varlist') - ///
								`c'*((`kmax'-`varlist')^3)*(`kmax'>`varlist') + ///
                                                                 `c'*((`kmin'-`varlist')^3)*(`kmin'>`varlist') if `touse'
                        }
                        else {
                                scalar `gen'1 = (`kmax'-`scalar')^2*(`kmax'>`scalar') - ///
								`c'*((`kmax'-`scalar')^3)*(`kmax'>`scalar') + ///
                                                                 `c'*((`kmin'-`scalar')^3)*(`kmin'>`scalar') 
                        }
                        local rcslist `gen'1

/* generate first derivative if dgen option is specified */
                        if "`dgen'" != "" {
                                if "`scalar'" == "" {
                                        quietly gen double `dgen'1 = -2*(`kmax'-`varlist')*(`kmax'>`varlist') - ///
								 (-3)*`c'*((`kmax'-`varlist')^2)*(`kmax'>`varlist') + ///
                                                                 (-3)*`c'*((`kmin'-`varlist')^2)*(`kmin'>`varlist') if `touse' 
                                }
                                else {
                                        scalar `dgen'1 = -2*(`kmax'-`scalar')*(`kmax'>`scalar') - ///
								 (-3)*`c'*((`kmax'-`scalar')^2)*(`kmax'>`scalar') + ///
                                                                 (-3)*`c'*((`kmin'-`scalar')^2)*(`kmin'>`scalar')
                                }
                                local drcslist `dgen'1
                        }


                        forvalues j=2/`npar' {
                                local h = `nk'-(`j'-1)
                                local lambda = (`k`h''-`kmin')/(`kmax' - `kmin')
                                if "`scalar'" == "" {
                                        quietly gen double `gen'`j' = ((`k`h''-`varlist')^3)*(`k`h''>`varlist') - ///
                                                                `lambda'*((`kmax'-`varlist')^3)*(`kmax'>`varlist') - ///
                                                                (1-`lambda')*((`kmin'-`varlist')^3)*(`kmin'>`varlist') if `touse'
                                }
                                else {
                                        scalar `gen'`j' = ((`k`h''-`scalar')^3)*(`k`h''>`scalar') - ///
                                                                `lambda'*((`kmax'-`scalar')^3)*(`kmax'>`scalar') - ///
                                                                (1-`lambda')*((`kmin'-`scalar')^3)*(`kmin'>`scalar') 
                                }
                                local rcslist `rcslist' `gen'`j'

/* calculate derivatives */
                                if "`dgen'"!="" {
                                        if "`scalar'" == "" {
                                                quietly gen double `dgen'`j' = (-3*(`k`h''-`varlist')^2)*(`k`h''>`varlist') - ///
                                                                        `lambda'*(-3*(`kmax'-`varlist')^2)*(`kmax'>`varlist') - ///
                                                                        (1-`lambda')*(-3*(`kmin'-`varlist')^2)*(`kmin'>`varlist')  if `touse'
                                        }
                                        else {
                                                scalar `dgen'`j' = (-3*(`k`h''-`scalar')^2)*(`k`h''>`scalar') - ///
                                                                        `lambda'*(-3*(`kmax'-`scalar')^2)*(`kmax'>`scalar') - ///
                                                                        (1-`lambda')*(-3*(`kmin'-`scalar')^2)*(`kmin'>` scalar') 

                                        }
                                        local drcslist `drcslist' `dgen'`j'
                                }       
                        }
/* Derive last spline variable */
                        if "`scalar'" == "" {
                                quietly gen double `gen'`nparams' = `varlist' if `touse'
                        }
                        else {
                                scalar `gen'`nparams' = `scalar' 
                        }
                        local rcslist `rcslist' `gen'`nparams'

/* generate first derivative if dgen option is specified */
                        if "`dgen'" != "" {
                                if "`scalar'" == "" {
                                        quietly gen double `dgen'`nparams' = 1 if `touse'
                                }
                                else {
                                        scalar `dgen'`nparams' = 1 
                                }
                                local drcslist `drcslist' `dgen'`nparams'
                        }
                }
        }
 
**************************** !!!! ************************************ 

******** ADDED 2012-01-10 Relax assumption of continuous first and second derivative at the first knot (or the last knot if reverse option is used)

****************************!!!*********************	
		
	else if "`nosecondder'" == "" & "`reverse'" == "" & "`nofirstder'" != "" {
  				
				/* Start to derive spline variables */
                if "`scalar'" == "" {
                        quietly gen double `gen'1 = `varlist' if `touse'
                }
                else {
                        scalar `gen'1 = `scalar'
                }                       

        /* generate first derivative if dgen option is specified */
                if "`dgen'" != "" {
                        if "`scalar'" == "" {
                                quietly gen double `dgen'1 = 1 if `touse'
                        }
                        else {
                                scalar `dgen'1 = 1
                        }                       
                }
                local rcslist `gen'1 
                local drcslist `dgen'1

                local nk : word count `knots'
                if "`knots'" == "" {
                        local interior  = 0
                }
                else {
                        local interior  = `nk' - 2
                }
                local nparams = `interior' + 3
				local npar = `interior' + 1
				local par = `interior' + 2

                if "`knots'" != "" {
                        local i = 1 
                        tokenize "`knots'"
                        while "``i''" != "" {
                                local k`i' ``i''
                                local i = `i' + 1
                        }

                        local kmin = `k1'
                        local kmax = `k`nk''

                        forvalues j=2/`npar' {
                                local lambda = (`kmax' - `k`j'')/(`kmax' - `kmin')
                                if "`scalar'" == "" {
                                        quietly gen double `gen'`j' = ((`varlist'-`k`j'')^3)*(`varlist'>`k`j'') - ///
                                                                `lambda'*((`varlist'-`kmin')^3)*(`varlist'>`kmin') - ///
                                                                (1-`lambda')*((`varlist'-`kmax')^3)*(`varlist'>`kmax') if `touse'
                                }
                                else {
                                        scalar `gen'`j' = ((`scalar'-`k`j'')^3)*(`scalar'>`k`j'') - ///
                                                                `lambda'*((`scalar'-`kmin')^3)*(`scalar'>`kmin') - ///
                                                                (1-`lambda')*((`varlist'-`kmax')^3)*(`scalar'>`kmax')  

                                }
                                
                                local rcslist `rcslist' `gen'`j'
        
        /* calculate derivatives */
                                if "`dgen'"!="" {
                                        if "`scalar'" == "" {
                                                quietly gen double `dgen'`j' = (3*(`varlist'- `k`j'')^2)*(`varlist'>`k`j'') - ///
                                                                        `lambda'*(3*(`varlist'-`kmin')^2)*(`varlist'>`kmin') - ///
                                                                        (1-`lambda')*(3*(`varlist'-`kmax')^2)*(`varlist'>`kmax') 
                                        }
                                        else {
                                                scalar `dgen'`j' = (3*(`scalar'- `k`j'')^2)*(`scalar'>`k`j'') - ///
                                                                        `lambda'*(3*(`scalar'-`kmin')^2)*(`scalar'>`kmin') - ///
                                                                        (1-`lambda')*(3*(`scalar'-`kmax')^2)*(`scalar'> `kmax') 
                                        }
                                        local drcslist `drcslist' `dgen'`j'
                                }       
                        }
				
		/* Derive the first extra spline variable */
						local c=(1/(3*(`kmax' - `kmin')))						
                        if "`scalar'" == "" {
                                quietly gen double `gen'`par' = (`varlist'-`kmin')^2*(`varlist'>`kmin') - ///
								`c'*((`varlist'-`kmin')^3)*(`varlist'>`kmin') + ///
                                                                 `c'*((`varlist'-`kmax')^3)*(`varlist'>`kmax') if `touse'
                        }
                        else {
                                scalar `gen'`par' = (`scalar'-`kmin')^2*(`scalar'>`kmin') - ///
								`c'*((`scalar'-`kmin')^3)*(`scalar'>`kmin') + ///
                                                                 `c'*((`scalar'-`kmax')^3)*(`scalar'>`kmax') 
                        }
                        local rcslist `rcslist' `gen'`par'

/* generate first derivative if dgen option is specified */
                        if "`dgen'" != "" {
                                if "`scalar'" == "" {
                                        quietly gen double `dgen'`par' = 2*(`varlist'-`kmin')*(`varlist'>`kmin') - ///
								 3*`c'*((`varlist'-`kmin')^2)*(`varlist'>`kmin') + ///
                                                                 3*`c'*((`varlist'-`kmax')^2)*(`varlist'>`kmax') if `touse' 
                                }
                                else {
                                        scalar `dgen'`par' = 2*(`scalar'-`kmin')*(`scalar'>`kmin') - ///
								 3*`c'*((`scalar'-`kmin')^2)*(`scalar'>`kmin') + ///
                                                                 3*`c'*((`scalar'-`kmax')^2)*(`scalar'>`kmax') 
                                }
                                local drcslist `drcslist' `dgen'`par'
                        }
		/* Derive the last spline variable */
						if "`scalar'" == "" {
                                quietly gen double `gen'`nparams' = (`varlist'-`kmin')*(`varlist'>`kmin') if `touse'
                        }
                        else {
                                scalar `gen'`nparams' = (`scalar'-`kmin')*(`scalar'>`kmin')  
                        }
                        local rcslist `rcslist' `gen'`nparams'

/* generate first derivative if dgen option is specified */
                        if "`dgen'" != "" {
                                if "`scalar'" == "" {
                                        quietly gen double `dgen'`nparams' = 1*(`varlist'>`kmin')  if `touse' 
                                }
                                else {
                                        scalar `dgen'`nparams' = 1*(`scalar'>`kmin')
                                }
                                local drcslist `drcslist' `dgen'`nparams'
                        }
                }
        }
 
 ****************************!!!*********************			
		
	else if "`nosecondder'" == "" & "`reverse'" != "" & "`nofirstder'" != ""{
                			
                local nk : word count `knots'
                if "`knots'" == "" {
                        local interior  = 0
                }
                else {
                        local interior  = `nk' - 2
                }
                local nparams = `interior' + 3
				local npar = `interior' + 1
				local par = `interior' + 2

                if "`knots'" != "" {
                        local i = 1 
                        tokenize "`knots'"
                        while "``i''" != "" {
                                local k`i' ``i''
                                local i = `i' + 1
                        }

                        local kmin = `k1'
                        local kmax = `k`nk''
	/* Derive first spline variable */
                       if "`scalar'" == "" {
                                quietly gen double `gen'1 = (`kmax'-`varlist')*(`kmax'>`varlist')  if `touse'
                        }
                        else {
                                scalar `gen'1 = (`kmax'-`scalar')*(`kmax'>`scalar')  
                        }
                        local rcslist `gen'1

/* generate first derivative if dgen option is specified */
                        if "`dgen'" != "" {
                                if "`scalar'" == "" {
                                        quietly gen double `dgen'1 = -1*(`kmax'>`varlist') if `touse' 
                                }
                                else {
                                        scalar `dgen'1 = -1
                                }
                                local drcslist `dgen'1
                        }

	/* Derive second spline variable */
                        local c=(1/(3*(`kmax' - `kmin')))						
                        if "`scalar'" == "" {
                                quietly gen double `gen'2 = (`kmax'-`varlist')^2*(`kmax'>`varlist') - ///
								`c'*((`kmax'-`varlist')^3)*(`kmax'>`varlist') + ///
                                                                 `c'*((`kmin'-`varlist')^3)*(`kmin'>`varlist') if `touse'
                        }
                        else {
                                scalar `gen'2 = (`kmax'-`scalar')^2*(`kmax'>`scalar') - ///
								`c'*((`kmax'-`scalar')^3)*(`kmax'>`scalar') + ///
                                                                 `c'*((`kmin'-`scalar')^3)*(`kmin'>`scalar') 
                        }
                        local rcslist `rcslist' `gen'2

/* generate first derivative if dgen option is specified */
                        if "`dgen'" != "" {
                                if "`scalar'" == "" {
                                        quietly gen double `dgen'2 = -2*(`kmax'-`varlist')*(`kmax'>`varlist') - ///
								 (-3)*`c'*((`kmax'-`varlist')^2)*(`kmax'>`varlist') + ///
                                                                 (-3)*`c'*((`kmin'-`varlist')^2)*(`kmin'>`varlist') if `touse' 
                                }
                                else {
                                        scalar `dgen'2 = -2*(`kmax'-`scalar')*(`kmax'>`scalar') - ///
								 (-3)*`c'*((`kmax'-`scalar')^2)*(`kmax'>`scalar') + ///
                                                                 (-3)*`c'*((`kmin'-`scalar')^2)*(`kmin'>`scalar')
                                }
                                local drcslist `drcslist' `dgen'2
                        }


                        forvalues j=3/`par' {
                                local h = `nk'-(`j'-2)
                                local lambda = (`k`h''-`kmin')/(`kmax' - `kmin')
                                if "`scalar'" == "" {
                                        quietly gen double `gen'`j' = ((`k`h''-`varlist')^3)*(`k`h''>`varlist') - ///
                                                                `lambda'*((`kmax'-`varlist')^3)*(`kmax'>`varlist') - ///
                                                                (1-`lambda')*((`kmin'-`varlist')^3)*(`kmin'>`varlist') if `touse'
                                }
                                else {
                                        scalar `gen'`j' = ((`k`h''-`scalar')^3)*(`k`h''>`scalar') - ///
                                                                `lambda'*((`kmax'-`scalar')^3)*(`kmax'>`scalar') - ///
                                                                (1-`lambda')*((`kmin'-`scalar')^3)*(`kmin'>`scalar') 
                                }
                                local rcslist `rcslist' `gen'`j'

/* calculate derivatives */
                                if "`dgen'"!="" {
                                        if "`scalar'" == "" {
                                                quietly gen double `dgen'`j' = (-3*(`k`h''-`varlist')^2)*(`k`h''>`varlist') - ///
                                                                        `lambda'*(-3*(`kmax'-`varlist')^2)*(`kmax'>`varlist') - ///
                                                                        (1-`lambda')*(-3*(`kmin'-`varlist')^2)*(`kmin'>`varlist')  if `touse'
                                        }
                                        else {
                                                scalar `dgen'`j' = (-3*(`k`h''-`scalar')^2)*(`k`h''>`scalar') - ///
                                                                        `lambda'*(-3*(`kmax'-`scalar')^2)*(`kmax'>`scalar') - ///
                                                                        (1-`lambda')*(-3*(`kmin'-`scalar')^2)*(`kmin'>` scalar') 

                                        }
                                        local drcslist `drcslist' `dgen'`j'
                                }       
                        }
/* Derive last spline variable */
                        if "`scalar'" == "" {
                                quietly gen double `gen'`nparams' = `varlist' if `touse'
                        }
                        else {
                                scalar `gen'`nparams' = `scalar' 
                        }
                        local rcslist `rcslist' `gen'`nparams'

/* generate first derivative if dgen option is specified */
                        if "`dgen'" != "" {
                                if "`scalar'" == "" {
                                        quietly gen double `dgen'`nparams' = 1 if `touse'
                                }
                                else {
                                        scalar `dgen'`nparams' = 1 
                                }
                                local drcslist `drcslist' `dgen'`nparams'
                        }
                }
        }
 
**************************** !!!! ************************************  
 
/* orthogonlise */      
	if "`orthog'" != "" {
		tempname R Rinv cons
		mata: orthgs("`rcslist'","`touse'") 
		matrix `Rinv' = inv(`R')
		if "`dgen'" != "" {
			gen `cons' = 1 if `touse'
			mata st_store(.,tokens(st_local("drcslist")), /// 
							"`touse'",st_data(.,tokens(st_local("drcslist")), ///
							"`touse'")*st_matrix("`Rinv'")[1..`nparams',1..`nparams'])
		}
	}
	else if "`rmatrix'" != "" {
		tempname Rinv cons
		matrix `Rinv' = inv(`rmatrix')
		if "`scalar'" == "" {
			gen `cons' = 1 if `touse'
			mata st_store(.,tokens(st_local("rcslist")), ///  
					"`touse'",(st_data(.,   tokens(st_local("rcslist") + " `cons'"), ///
					"`touse'"))*st_matrix("`Rinv'")[,1..`nparams'])	
			if "`dgen'" != "" {
				mata st_store(.,tokens(st_local("drcslist")), ///
					"`touse'",st_data(.,tokens(st_local("drcslist")), ///
					"`touse'")*st_matrix("`Rinv'")[1..`nparams',1..`nparams'])
		
			}
		}
		else {
			tempname scalarmatrix
			matrix `scalarmatrix' = `gen'1
			forvalues i = 2/`nparams'{
				matrix `scalarmatrix' = `scalarmatrix',`gen'`i'
			}
			matrix `scalarmatrix' = `scalarmatrix',1
			mata st_matrix("`scalarmatrix'",st_matrix("`scalarmatrix'")*st_matrix("`Rinv'")[,1..`nparams']) 
			forvalues i = 1/`nparams'{
				scalar `gen'`i' = el(`scalarmatrix',1,`i')
			}
			if "`dgen'" != "" {
				tempname dscalarmatrix
				matrix `dscalarmatrix' = `dgen'1
				forvalues i = 2/`nparams'{
					matrix `dscalarmatrix' = `dscalarmatrix',`dgen'`i'
				}
				mata st_matrix("`dscalarmatrix'",st_matrix("`dscalarmatrix'")*st_matrix("`Rinv'")[1..`nparams',1..`nparams'])
				forvalues i = 1/`nparams'{
					scalar `dgen'`i' = el(`dscalarmatrix',1,`i')
				}
			}
		}
	}
	
	/*subtract the value at the last knot, so that there is an interpretable baseline*/ 		/*ADDED: 2010-04-14 by Therese Andersson*/
	if ("`orthog'" != "" | "`rmatrix'" != "") & "`reverse'" != "" {					/*2010-04-28 TA, also make this work with rmatrix*/

		tempname rcsvaluevector
		local rcsvaluelist
		
		if "`knots'" != "" {
			local i = 1 
			tokenize "`knots'"
			while "``i''" != "" {
				local k`i' ``i''
				local i = `i' + 1
			}

			local kmin = `k1'
			local kmax = `k`nk''

		if "`nosecondder'" == "" & "`nofirstder'" == "" {
			forvalues j=1/`interior' {
				local h = `nk'-`j'
				local lambda = (`k`h''-`kmin')/(`kmax' - `kmin')
				local rcsvalue`j' = ((`k`h''-`kmax')^3)*(`k`h''>`kmax') - ///
												`lambda'*((`kmax'-`kmax')^3)*(`kmax'>`kmax') - ///
												(1-`lambda')*((`kmin'-`kmax')^3)*(`kmin'>`kmax')
				local rcsvaluelist `rcsvaluelist' `rcsvalue`j''
			}
		}
						if "`nosecondder'" != "" & "`nofirstder'" == "" {
							local c=(1/(3*(`kmax' - `kmin')))						
							local rcsvalue1 = (`kmax'-`kmax')^2*(`kmax'>`kmax') - ///
												  `c'*((`kmax'-`kmax')^3)*(`kmax'>`kmax') + ///
                                                  `c'*((`kmin'-`kmax')^3)*(`kmin'>`kmax') 
							local rcsvaluelist `rcsvalue1'
							
							forvalues j=2/`npar' {
									local h = `nk'-(`j'-1)
									local lambda = (`k`h''-`kmin')/(`kmax' - `kmin')
                                    local rcsvalue`j' = ((`k`h''-`kmax')^3)*(`k`h''>`kmax') - ///
                                                         `lambda'*((`kmax'-`kmax')^3)*(`kmax'>`kmax') - ///
                                                         (1-`lambda')*((`kmin'-`kmax')^3)*(`kmin'>`kmax') 
									local rcsvaluelist `rcsvaluelist' `rcsvalue`j''
							}
						}
						
						
						if "`nosecondder'" == "" & "`nofirstder'" != "" {
							local rcsvalue1 = (`kmax'-`kmax')*(`kmax'>`kmax')
							local rcsvaluelist `rcsvalue1'
							
							local c=(1/(3*(`kmax' - `kmin')))
							local rcsvalue2 = (`kmax'-`kmax')^2*(`kmax'>`kmax') - ///
								        `c'*((`kmax'-`kmax')^3)*(`kmax'>`kmax') + ///
                                        `c'*((`kmin'-`kmax')^3)*(`kmin'>`kmax') 
							local rcsvaluelist `rcsvaluelist' `rcsvalue2'
							
							forvalues j=3/`par' {
                                local h = `nk'-(`j'-2)
                                local lambda = (`k`h''-`kmin')/(`kmax' - `kmin')
                                local rcsvalue`j' = ((`k`h''-`kmax')^3)*(`k`h''>`kmax') - ///
                                            `lambda'*((`kmax'-`kmax')^3)*(`kmax'>`kmax') - ///
                                            (1-`lambda')*((`kmin'-`kmax')^3)*(`kmin'>`kmax')
								local rcsvaluelist `rcsvaluelist' `rcsvalue`j''
							}
						}
						
                }
		/* Derive last spline variable */
		local rcsvaluelist `rcsvaluelist' `kmax'
		
		matrix input rcsvaluevector=(`rcsvaluelist' 1)
		matrix rcsvalueorthog=rcsvaluevector*`Rinv'
		
		if "`scalar'" == "" {
			forvalues j=1/`nparams' {
				qui replace `gen'`j' = `gen'`j' - rcsvalueorthog[1,`j']
			}
		}
		else {
			forvalues j=1/`nparams' {
				scalar `gen'`j' = `gen'`j' - rcsvalueorthog[1,`j']
			}
		}
	}

/* report new variables created */        
	if "`scalar'" != "" {
		local type Scalars
	}
	else {
		local type Variables
	}
	if "`dgen'"!="" {
		di in green "`type' `gen'1 to `gen'`nparams' and `dgen'1 to `dgen'`nparams' were created"
	}
	else {
		di in green "`type' `gen'1 to `gen'`nparams' were created"
	}
	if "`knots'" == "" {
		di in green "Warning: Only `gen'1 has been created as you did not specifiy any the knots, df or percentile options"
	}
	
	if "`orthog'" != "" {
		return matrix R = `R'
	}
	return local knots `knots'
end

/* Gram-Schmidt orthogonalization in Mata */        
mata:
void orthgs(string scalar varlist, |string scalar touse) 
{
	x = st_data(.,tokens(varlist),touse)
	meanx = mean(x)
	v = x :- meanx ,J(rows(x),1,1) 
	q = J(rows(v),0,.)
	R = J(cols(v),cols(v),0)
	R[cols(v),] = (meanx,1)
	for (i=1;i<=cols(x);i++){
		r = norm(v[,i])/sqrt(rows(v))
		q = q, (v[,i]:/ r)
		R[i,i] = r
		for (j = i + 1; j<=cols(x); j++){
			r = (q[,i]' * v[,j])/rows(v)
			v[,j] = v[,j] - r*q[,i]
			R[i,j] = r 
		}
	}
	st_store(.,tokens(varlist),touse,q)
	st_local("R",Rname=st_tempname())
	st_matrix(Rname,R)
}
end
	
