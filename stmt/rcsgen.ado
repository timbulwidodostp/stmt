*! version 1.5 24oct2010
* Added various options + use mata for orthogonalisation
* Based on an original program by Chris Nelson
* Chris Nelson 24/APR/2006
* Paul Lambert 28/AUG/2008
* Mark Rutherford 05/DEC/2008
* Patrick royston 16/OCT2009 - allow numlist for knots to be in non-ascending order
* Therese Andersson 01/FEB/2010 - add an option to calculate the spline variables backwards, to use for cure estimation in stpm2
* Patrick royston 24/OCT/2010 - return in `r(rcslist)' (and `r(drcslist)') names of created variables (and derivatives, if specified)
program define rcsgen, rclass
	version 10.0
	syntax  varlist(max=1) [if] [in] ///
		, [Gen(string) DGen(string) Knots(numlist) BKnots(numlist max=2) Orthog Percentiles(numlist ascending) RMATrix(name) ///
			DF(int 0)  IF2(string) FW(varname)  REVerse]      

	marksample touse
	
	// sort knots
	if "`knots'" != "" {
		numlist "`knots'", sort
		local knots `r(numlist)'
	}
    
/* Error checks */	
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

	if "`df'" != "0" {
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
		
	if "`reverse'" == "" {
/* Start to derive spline variables */
		quietly gen double `gen'1 = `varlist' if `touse'

/* generate first derivative if dgen option is specified */
		if "`dgen'" != "" {
			quietly gen double `dgen'1 = 1 if `touse'
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
				quietly gen double `gen'`j' = ((`varlist'-`k`j'')^3)*(`varlist'>`k`j'') - ///
											`lambda'*((`varlist'-`kmin')^3)*(`varlist'>`kmin') - ///
											(1-`lambda')*((`varlist'-`kmax')^3)*(`varlist'>`kmax')  if `touse'
				local rcslist `rcslist' `gen'`j'
	
/* calculate derivatives */
				if "`dgen'"!="" {
					quietly gen double `dgen'`j' = (3*(`varlist'- `k`j'')^2)*(`varlist'>`k`j'') - ///
												`lambda'*(3*(`varlist'-`kmin')^2)*(`varlist'>`kmin') - ///
												(1-`lambda')*(3*(`varlist'-`kmax')^2)*(`varlist'>`kmax')  if `touse'
					local drcslist `drcslist' `dgen'`j'
				}       
			}
		}
	}
/*Derive the spline variables in reversed order */		/*ADDED: 2010-02-02 by Therese Andersson*/
		
	else if "`reverse'" != "" {
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
				quietly gen double `gen'`j' = ((`k`h''-`varlist')^3)*(`k`h''>`varlist') - ///
												`lambda'*((`kmax'-`varlist')^3)*(`kmax'>`varlist') - ///
												(1-`lambda')*((`kmin'-`varlist')^3)*(`kmin'>`varlist')  if `touse'
				local rcslist `rcslist' `gen'`j'

/* calculate derivatives */
				if "`dgen'"!="" {
					quietly gen double `dgen'`j' = (-3*(`k`h''-`varlist')^2)*(`k`h''>`varlist') - ///
													`lambda'*(-3*(`kmax'-`varlist')^2)*(`kmax'>`varlist') - ///
												(1-`lambda')*(-3*(`kmin'-`varlist')^2)*(`kmin'>`varlist')  if `touse'
					local drcslist `drcslist' `dgen'`j'
				}       
			}
/* Derive last spline variable */
			quietly gen double `gen'`nparams' = `varlist' if `touse'
			local rcslist `rcslist' `gen'`nparams'

/* generate first derivative if dgen option is specified */
			if "`dgen'" != "" {
				quietly gen double `dgen'`nparams' = 1 if `touse'
				local drcslist `drcslist' `dgen'`nparams'
			}
		}
	}
 
/* orthogonlise */      
	if "`orthog'" != "" {
		tempname R Rinv cons
		mata: orthgs("`rcslist'","`touse'") 
		if "`dgen'" != "" {
			gen `cons' = 1 if `touse'
			matrix `Rinv' = inv(`R')
			mata st_store(.,tokens(st_local("drcslist")), /// 
							"`touse'",st_data(.,tokens(st_local("drcslist")), ///
							"`touse'")*st_matrix("`Rinv'")[1..`nparams',1..`nparams'])
		}
	}
	else if "`rmatrix'" != "" {
		tempname Rinv cons
		matrix `Rinv' = inv(`rmatrix')
		gen `cons' = 1 if `touse'
		mata st_store(.,tokens(st_local("rcslist")),"`touse'", ///
						(st_data(.,tokens(st_local("rcslist") + " `cons'"), ///
						"`touse'"))*st_matrix("`Rinv'")[,1..`nparams'])
		if "`dgen'" != "" {
			mata st_store(.,tokens(st_local("drcslist")), ///
						"`touse'",st_data(.,tokens(st_local("drcslist")), ///
						"`touse'")*st_matrix("`Rinv'")[1..`nparams',1..`nparams'])
		}
	}
	
	
/* report new variables created */        
	if "`dgen'"!="" {
		di in green "Variables `gen'1 to `gen'`nparams' and `dgen'1 to `dgen'`nparams' were created"
	}
	else {
		di in green "Variables `gen'1 to `gen'`nparams' were created"
	}
	if "`knots'" == "" {
		di in green "Warning: Only `gen'1 has been created as you did not specifiy any the knots, df or percentile options"
	}
	
	if "`orthog'" != "" {
		return matrix R = `R'
	}
	return local knots `knots'
	// PR additions:
	return local rcslist `rcslist'
	if "`dgen'" != ""  return local drcslist `drcslist'
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
