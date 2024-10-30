{smcl}
{* 19oct2011}{...}
{cmd:help stsurvsim}{right: ({browse "http://www.stata-journal.com/article.html?article=up0043":SJ14-2: st0274_1})}
{hline}

{title:Title}

{p2colset 5 18 20 2}{...}
{p2col :{hi:stsurvsim} {hline 2}}Simulate survival times from Royston-Parmar
model{p_end}
{p2colreset}{...}
 

{title:Syntax}

{p 8 17 2}
{cmd:stsurvsim}
{newvar}
{ifin}
[{cmd:,} {it:options}]

{synoptset 22}{...}
{synopthdr}
{synoptline}
{synopt :{opt be:ta(item [item ...])}}regression coefficients for spline basis functions
and optionally for covariates{p_end}
{synopt :{opth bk:nots(numlist)}}two boundary knots on scale of time or log time{p_end}
{synopt :{opth k:nots(numlist)}}interior knots on scale of time or log time{p_end}
{synopt :{opt kns:cale(string)}}determine scale of knots (time or log time){p_end}
{synopt :{opt centol(#)}}convergence criterion{p_end}
{synopt :{opt max:iter(#)}}maximum number of iterations{p_end}
{synopt :{opt rmat:rix(matrix_name)}}supply spline basis orthogonalization matrix{p_end}
{synopt :{opt s:cale(scale_name)}}specify the scale of the Royston-Parmar model{p_end}
{synoptline}
{p2colreset}{...}

{pstd}
Weights are not allowed.  Note that {cmd:stsurvsim} is likely to
require {cmd:stpm2} to be installed -- it can be downloaded from the
Statistical Software Components archive (see {helpb ssc}).


{title:Description}

{pstd}
{cmd:stsurvsim} creates in {it:newvar} a simulated sample of
uncensored survival times from the baseline distribution function of a
Royston-Parmar flexible parametric survival model, implemented by 
{helpb stpm2}.  Covariate effects may optionally be included.

{pstd}
{cmd:stsurvsim} can be used in two different modes.  If any of the
options {opt beta()}, {opt bknots()}, {opt knots()}, or {opt scale()}
are supplied, then they must all be supplied as well as possibly 
{opt rmatrix()}.  If, however, a linear model is fit (no spline terms),
knots and boundary knots are not required, and only {opt beta()} and
{opt scale()} must be supplied.  {cmd:stsurvsim} uses their information
to define the model from which to simulate.  Otherwise, none of these
options are provided, and {cmd:stsurvsim} uses information from the most
recent fit of {cmd:stpm2}.

{pstd}
Important: If the beta coefficients for the spline basis functions
relate to orthogonalized basis functions, you must supply an R-matrix in
the {opt rmatrix()} option (see {it:{help stsurvsim##remarks:Remarks}}
below).  {cmd:stsurvsim} has no way of checking whether an R-matrix is
needed; it is up to you to provide it if necessary.


{title:Options}

{phang}
{cmd:beta(}{it:item} [{it:item} ...]{cmd:)} specifies the
regression coefficients for the spline basis functions and also for any
covariates that may be required.  The syntax of {it:item} is
{it:name}[{cmd:=}]{it:#} (where {it:name} is the name of a covariate),
{cmd:_rcs}{it:#} (where {it:#} is 1, 2, ... corresponding to the spline
basis functions in the model), or {cmd:_cons} (for the intercept term).
For example: {cmd:beta(_rcs1=2.751 _rcs2=0.229 _cons=-2.079)}.  If 
{opt beta()} is not specified, the regression coefficients are picked up
from the most recent fit of {cmd:stpm2}.

{phang}
{opth bknots(numlist)} specifies the two boundary knots for the
spline function on the scale of time or, if the {cmd:knscale(log)}
option is used, log time.  If {cmd:bknots()} is not specified, the
boundary knots are picked up from the most recent fit of {cmd:stpm2}.

{phang}
{opth knots(numlist)} specifies interior knots for the spline
function on the scale of time or, if the {cmd:knscale(log)} option is
used, log time.  If {cmd:knots()} is not specified, the interior knots
are picked up from the most recent fit of {cmd:stpm2}.

{phang}
{opt knscale(string)} determines the scale of the knots supplied
in {opt bknots()} and {opt knots()}.  If {opt knscale()} is not
specified, knots are given in units of time.  If {cmd:knscale(log)} is
specified, knots are given in units of log time.  The default is
{cmd:knscale("")}, meaning knots in units of time.

{phang}
{opt centol(#)} defines the convergence criterion.  This option is
seldom used.  At the heart of {cmd:stsurvsim} is an iterative procedure that
converts a pseudorandom sample from a standard uniform distribution into
times to event with a Royston-Parmar distribution.  Typically, convergence is
rapid; about four to six iterations are needed.  If convergence fails, a
message is displayed.  It may be possible to cure the problem by specifying a
larger value of {it:#} in {opt centol(#)}.  The default is
{cmd:centol(0.001)}.  Thus we suggest trying {cmd:centol(0.01)} in such
instances.  See also {opt maxiter(#)}.

{phang}
{opt maxiter(#)} specifies the maximum number of iterations of the
procedure used to estimate simulated times to event.  This option is seldom
used.  See also {opt centol(#)}.  The default is {cmd:maxiter(100)}.

{phang}
{opt rmatrix(matrix_name)} supplies the matrix used to
orthogonalize the spline basis functions.  If the {opt noorthog} option
of {cmd:stpm2} was not used when fitting the original model,
orthogonalization is done by default.  It is then essential to save the
R-matrix after running {cmd:stpm2} and supply it to {cmd:stsurvsim} as
{opt rmatrix(matrix_name)}.  {cmd:stpm2} stores the R-matrix in
{cmd:e(R_bh)}.

{phang}
{opt scale(scale_name)} specifies the scale on which the
Royston-Parmar model was fit.  Valid {it:scale_name}s are {cmd:hazard},
{cmd:odds}, and {cmd:normal}.  The {cmd:scale(theta)} option of
{cmd:stpm2} is not supported.  If {cmd:scale()} is not specified, the
scale is picked up from the most recent fit of {cmd:stpm2}.


{marker remarks}{...}
{title:Remarks}

{pstd}
If a Royston-Parmar model is fit by {cmd:stpm2}, the spline basis
functions, created in variables called {cmd:_rcs1}, {cmd:_rcs2}, etc.,
are orthogonalized by default.  A linear transformation to do this is
stored by {cmd:stpm2} as an R-matrix in {cmd:e(R_bh)}.  To obtain
correct simulations, you need to store this matrix after running
{cmd:stpm2} and supply it to {cmd:stsurvsim} in the 
{opt rmatrix(matrix_name)} option.  Failure to do this will produce
incorrect simulations.

{pstd}
You can avoid the bother of the R-matrix by specifying the
{cmd:noorthog} option when running {cmd:stpm2}.  The {opt rmatrix()}
option of {cmd:stsurvsim} is then unnecessary.


{title:Examples}

{phang}{cmd:. webuse brcancer}{p_end}
{phang}{cmd:. stset rectime, failure(censrec) scale(365.25)}{p_end}

{phang}{cmd:. stpm2, df(2) scale(hazard) noorthog}{p_end}
{phang}{cmd:. stsurvsim t_sim}{p_end}
{phang}{cmd:. stsurvsim t_sim, knots(1.77) bknots(.197 6.72) beta(_rcs1 = 2.75 _rcs2 = 0.23 _cons = -2.08) scale(hazard)}{p_end}
{phang}{cmd:. stsurvsim t_sim, knots(1.77) bknots(.197 6.72) beta(_rcs1 = 2.75 _rcs2 = 0.23 _cons = -2.08 hormon = -0.37) scale(hazard)}{p_end}

{phang}{cmd:. replace _d = !_d}{p_end}
{phang}{cmd:. stpm2, df(2) scale(hazard) lininit}{p_end}
{phang}{cmd:. stsurvsim t_cens}{p_end}

{phang}{cmd:. generate t = min(t_sim, t_cens)}{p_end}
{phang}{cmd:. generate byte d = cond(t_cens <= t_sim, 0, 1)}{p_end}
{phang}{cmd:. stset t, failure(d)}{p_end}

{phang}{cmd:. stpm2, df(1) scale(odds) noorthog}{p_end}
{phang}{cmd:. stsurvsim t}{p_end}
{phang}{cmd:. stsurvsim t, beta(_rcs1=1.53 _cons=-2.3) scale(odds)}{p_end}


{title:Author}

{pstd}Patrick Royston{p_end}
{pstd}MRC Clinical Trials Unit{p_end}
{pstd}London, UK{p_end}
{pstd}pr@ctu.mrc.ac.uk{p_end}


{title:Also see}

{p 4 14 2}Article:  {it:Stata Journal}, volume 14, number 2: {browse "http://www.stata-journal.com/article.html?article=up0043":st0274_1},{break}
                    {it:Stata Journal}, volume 12, number 4: {browse "http://www.stata-journal.com/article.html?article=st0274":st0274}

{p 7 14 2}Help:  {helpb stpm2} (if installed){p_end}
