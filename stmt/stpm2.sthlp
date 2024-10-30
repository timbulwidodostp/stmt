{smcl}
{* *! version 1.4.4 10Mar2012}{...}
{cmd:help stpm2}{right: ({browse "http://www.stata-journal.com/article.html?article=st0165_1":SJ12-4: st0165_1})}
{hline}

{title:Title}

{p2colset 5 14 16 2}{...}
{p2col :{hi:stpm2} {hline 2}}Flexible parametric survival models{p_end}
{p2colreset}{...}


{title:Syntax}

{p 8 16 2}{cmd:stpm2} [{varlist}] {ifin} [{cmd:,} {it:options}]

{marker options}{...}
{synoptset 22}{...}
{synopthdr}
{synoptline}
{synopt :{opth bhaz:ard(varname)}}invoke relative survival models where
{it:varname} holds the expected mortality rate (hazard) at the time of death
{p_end}
{synopt :{opt bk:nots(knotslist)}}boundary knots for baseline{p_end}
{synopt :{opt bknotstvc(knotslist)}}boundary knots for time-dependent effects{p_end}
{synopt :{opt cure}}fit a cure model{p_end}
{synopt :{opt df(#)}}specify degrees of freedom for baseline hazard function{p_end}
{synopt :{opt dft:vc(df_list)}}specify degrees of freedom for each time-dependent effect{p_end}
{synopt :{opt failconvlininit}}automatically try {cmd:lininit} option if convergence fails{p_end}
{synopt :{opth knots(numlist)}}specify knot locations for baseline hazard{p_end}
{synopt :{opt knotst:vc(knotslist)}}specify knot locations for time-dependent effects{p_end}
{synopt :{opt knscale(scale)}}scale for user-defined knots (default scale is time){p_end}
{synopt :{opt nocons:tant}}suppress constant term{p_end}
{synopt :{opt rcsbaseoff}}do not include baseline spline variables{p_end}
{synopt :{opt noorth:og}}do not use orthogonal transformation of spline variables{p_end}
{synopt :{opt sc:ale(scalename)}}specify scale on which survival model is
 to be fit{p_end}
{synopt :{opth st:ratify(varlist)}}for backward compatibility with {cmd:stpm}{p_end}
{synopt :{cmdab:th:eta(est}|{it:#}{cmd:)}}for backward compatibility with
{cmd:stpm}{p_end}
{synopt :{opth tvc(varlist)}}{it:varlist} of time-varying effects{p_end}
{synopt :{opt alleq}}report all equations{p_end}
{synopt :{opt ef:orm}}report exponentiate coefficients{p_end}
{synopt :{opt keepc:ons}}do not drop constraints used in ml routine{p_end}
{synopt :{opt l:evel(#)}}set confidence level; default is {cmd:level(95)}{p_end}
{synopt :{opt showc:ons}}list constraints in output{p_end}
{synopt :{opt const:heta(#)}}constrain value of theta when using Aranda-Ordaz family of link functions{p_end}
{synopt :{opt initt:heta(#)}}initial value of theta (default 1: log cumulative-odds scale){p_end}
{synopt :{opt lin:init}}obtain initial values by first fitting a linear function of ln(time){p_end}
{synopt :{it:{help stpm2##maximize_options:maximize_options}}}control the maximization process; seldom used{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}
You must {cmd:stset} your data before using {cmd:stpm2}; see {manhelp stset ST}.{p_end}
{p 4 6 2}
{cmd:fweight}s, {cmd:iweight}s, and {cmd:pweight}s may be specified using
{cmd:stset}; see {manhelp stset ST}.{p_end}


{title:Description}

{pstd}
{cmd:stpm2} fits flexible parametric survival models (Royston-Parmar models).
{cmd:stpm2} can be used with single- or multiple-record or single- or
multiple-failure {cmd:st} data.  Survival models can be fit on the log
cumulative-hazard scale, the log cumulative-odds scale, the standard normal
deviate (probit) scale, or on a scale defined by the value of {cmd:theta()} by
using the Aranda-Ordaz family of link functions.

{pstd}
{cmd:stpm2} can fit the same models as {cmd:stpm} but is more
flexible in that it does not force the knots for time-dependent effects
to be the same as those used for the baseline distribution function.  In
addition, {cmd:stpm2} can fit relative survival models with the
{cmd:bhazard()} option.  Postestimation commands have been extended over
what is available in {cmd:stpm}.  {cmd:stpm2} is noticeably faster than
{cmd:stpm}.

{pstd}
See {manhelp streg ST} for other (standard) parametric survival models.


{title:Options}

{phang}
{opth bhazard(varname)} is used when fitting relative survival models.
{it:varname} gives the expected mortality rate at the time of death or
censoring.  {cmd:stpm2} gives an error message when there are missing values
of {it:varname}, because this usually indicates that an error has occurred
when merging the expected mortality rates.

{phang}
{opt bknots(knotslist)} is a two-element {it:knotslist} giving the boundary
knots.  By default, these are located at the minimum and maximum of the
uncensored survival times.  They are specified on the scale defined by
{cmd:knscale()}.

{phang}
{opt bknotstvc(knotslist)} gives the boundary knots for any time-dependent
effects.  By default, these are the same as for the {cmd:bknots()} option.
They are specified on the scale defined by {cmd:knscale()}.

{pmore}
For example, {cmd:bknotstvc(x1 0.01 10 x2 0.01 8)}.

{phang}
{opt cure} is used when fitting cure models.  It forces the cumulative hazard
to be constant after the last knot.  When the {cmd:df()} option is used
together with the {cmd:cure} option, the internal knots are placed evenly
according to centiles of the distribution of the uncensored log survival-times
except one, which is placed at the 95th centile.  Alternative knot locations
can be selected using the {cmd:knots()} option.  Cure models can only be used
when modeling on the log cumulative-hazard scale ({cmd:scale(hazard)}).

{phang}
{opt df(#)} specifies the degrees of freedom for the restricted cubic spline
function used for the baseline function.  {it:#} must be between 1 and 10, but
usually a value between 1 and 4 is sufficient, with 3 being the default.  The
{cmd:knots()} option is not applicable if the {cmd:df()} option is specified.
The knots are placed at the following centiles of the distribution of the
uncensored log survival times:

        {hline 60}
        df  knots  Centile positions
        {hline 60}
         1    0    (no knots)
         2    1    50
         3    2    33 67
         4    3    25 50 75
         5    4    20 40 60 80
         6    5    17 33 50 67 83
         7    6    14 29 43 57 71 86
         8    7    12.5 25 37.5 50 62.5 75 87.5
         9    8    11.1 22.2 33.3 44.4 55.6 66.7 77.8 88.9
        10    9    10 20 30 40 50 60 70 80 90
        {hline 60}
	
{pmore}
These are interior knots; there are also boundary knots placed at the minimum
and maximum of the distribution of uncensored survival times.

{pmore}
When the {cmd:cure} option is used, {cmd:df()} must be between 3 and 11, and
the default location of the knots are as follows:

        {hline 60}
        df  knots  Centile positions
        {hline 60}
         3    2    50 95
         4    3    33 67 95
         5    4    25 50 75 95
         6    5    20 40 60 80 95
         7    6    17 33 50 67 83 95
         8    7    14 29 43 57 71 86 95
         9    8    12.5 25 37.5 50 62.5 75 87.5 95
        10    9    11.1 22.2 33.3 44.4 55.6 66.7 77.8 88.9 95
        11   10    10 20 30 40 50 60 70 80 90 95		
        {hline 60}

{phang}
{opt dftvc(df_list)} gives the degrees of freedom for time-dependent effects
in {it:df_list}.  The potential degrees of freedom are listed under the
{opt df()} option.  With 1 degree of freedom, a linear effect of log time is
fit.  If there is more than one time-dependent effect and if different degrees
of freedom are requested for each time-dependent effect, then the following
syntax applies:

{pmore}
{cmd:dftvc(x1:3 x2:2 1)}

{pmore}
This will use 3 degrees of freedom for {cmd:x1}, 2 degrees of freedom for
{cmd:x2}, and 1 degree of freedom for all remaining time-dependent effects.

{phang}
{opt failconvlininit} automatically tries the {cmd:lininit} option if the
model fails to converge.

{phang}
{opth knots(numlist)} specifies knot locations for the baseline distribution
function, as opposed to the default locations set by {cmd:df()}.  The
locations of the knots are placed on the scale defined by {cmd:knscale()}.
However, the scale used by the restricted cubic spline function is always log
time.  Default knot positions are determined by the {opt df()} option.

{phang}
{opt knotstvc(knotslist)} defines {it:knotslist} as the location of the
interior knots for time-dependent effects.  If different knots are required
for different time-dependent effects, then the option is specified, for
example, as follows:

{pmore}
{cmd:knotstvc(x1 1 2 3 x2 1.5 3.5)}

{phang}
{opt knscale(scale)} sets the scale on which user-defined knots are specified.
{cmd:knscale(time)} denotes the original time scale, and {cmd:knscale(log)}
denotes the log time scale.  {cmd:knscale(centile)} specifies that the knots
be taken as centile positions in the distribution of the uncensored log
survival times.  The default is {cmd:knscale(time)}.

{phang}
{opt noconstant} suppresses the constant term (intercept) in the model.

{phang}
{cmd:rcsbaseoff} drops baseline spline variables from the model.  With this
option, you will generally want to specify your baseline separately in two or
more strata.  For example, the following code will fit a separate baseline
hazard for males and for females:

{pmore}
{cmd:stpm2 males females, scale(hazard) tvc(males females) dftvc(3) nocons rcsbaseoff}

{pmore}
Note that identical fitted values would be obtained if we instead used the
following:

{pmore}
{cmd:stpm2 females, df(3) scale(hazard) tvc(females) dftvc(3)}

{phang}
{cmd:noorthog} suppresses the orthogonal transformation of spline variables.

{phang}
{opt scale(scalename)} specifies the scale on which the survival model is to
be fit.

{pmore}
{cmd:scale({ul:h}azard)} fits a model on the log cumulative-hazard scale, that
is, the scale of ln[-ln{S(t)}].  If no time-dependent effects are specified,
the resulting model has proportional hazards.

{pmore}
{cmd:scale({ul:o}dds)} fits a model on the log cumulative-odds scale, that is,
ln[{1 - S(t)}/S(t)].  If no time-dependent effects are specified, then this
gives a proportional odds model.

{pmore}
{cmd:scale({ul:n}ormal)} fits a model on the normal equivalent deviate scale,
that is, a probit link for the survival function, invnorm{1 - S(t)}.

{pmore}
{cmd:scale({ul:t}heta)} fits a model on a scale defined by the value of theta
for the Aranda-Ordaz family of link functions, that is, ln[{S(t)^(-theta) -
1}/theta].  theta = 1 corresponds to a proportional-odds model, and theta = 0
corresponds to a proportional cumulative-hazards model.

{phang}
{opth stratify(varlist)} is provided for compatibility with {helpb stpm}.
Members of {it:varlist} are modeled with time-dependent effects.  See the
{opt tvc()} and {opt dftvc()} options for {cmd:stpm2}'s way of specifying
time-dependent effects.

{phang}
{cmd:theta(}{cmd:est}|{it:#}{cmd:)} is provided for compatibility with
{helpb stpm}.  {cmd:est} requests that theta be estimated, whereas {it:#}
fixes theta to {it:#}.  See {opt constheta()} and {cmd:inittheta()} for
{cmd:stpm2}'s way of specifying theta.

{phang}
{opth tvc(varlist)} gives the name of the variables that are time dependent.
Time-dependent effects are fit using restricted cubic splines.  The degrees of
freedom are specified using the {cmd:dftvc()} option.

{phang}
{opt alleq} reports all equations used by {cmd:ml}.  The models are fit by
using various constraints for parameters associated with the derivatives of
the spline functions.  These parameters are generally not of interest and thus
are not shown by default.  In addition, an extra equation is used when fitting
delayed-entry models, and again this is not shown by default.

{phang}
{opt eform} reports the exponentiated coefficients.  For models on the log
cumulative-hazard scale {cmd:scale(hazard)}, this gives hazard ratios if the
covariate is not time dependent.  Similarly, for models on the log
cumulative-odds scale {cmd:scale(odds)}, this option will give odds ratios for
non-time-dependent effects.

{phang}
{opt keepcons} prevents the constraints imposed by {cmd:stpm2} on the
derivatives of the spline function when fitting delayed-entry models from
being dropped.  By default, the constraints are dropped.

{phang}
{opt level(#)} specifies the confidence level, as a percentage, for confidence
intervals.  The default is {cmd:level(95)} or as set by {helpb set level}.

{phang}
{opt showcons} lists the output and the constraints used by {cmd:stpm2} for
the derivatives of the spline function and when fitting delayed-entry models;
the default is to not list them.

{phang}
{opt constheta(#)} constrains the value of theta; that is, it is treated as a
known constant.

{phang}
{opt inittheta(#)} gives an initial value for theta in the Aranda-Ordaz family
of link functions.

{phang}
{opt lininit} obtains initial values by fitting only the first spline basis
function (that is, a linear function of log survival time).  This option is
seldom needed.

{phang}{it:maximize_options}:  {opt dif:ficult}, 
{opt tech:nique(algorithm_spec)}, {opt iter:ate(#)}, 
[{cmdab:no:}]{opt lo:g}, {opt tr:ace}, {opt grad:ient}, {opt showstep},
{opt hess:ian}, {opt shownr:tolerance}, {opt tol:erance(#)}, 
{opt ltol:erance(#)}, {opt gtol:erance(#)}, {opt nrtol:erance(#)}, 
{opt nonrtol:erance}, and {opt from(init_specs)}; see 
{manhelp maximize R}.  These options are seldom used, but the 
{opt difficult} option may be useful if there are convergence problems
when fitting models that use the Aranda-Ordaz family of link functions.


{title:Remarks}

{pstd}
Let t denote time.  {cmd:stpm2} works by first calculating the survival
function after fitting a Cox proportional hazards model.  The procedure is
illustrated for proportional hazards models, specified by the
{cmd:scale(hazard)} option.  S(t) is converted to an estimate of the log
cumulative-hazard function, Z(t), by the formula

{pin}
Z(t) = ln[-ln{S(t)}]

{pstd}
This estimate of Z(t) is then smoothed on ln(t) by using regression splines
with knots placed at certain quantiles of the distribution of t.  The knot
positions are chosen automatically if the spline complexity is specified by
the {cmd:df()} option or manually by way of the {cmd:knots()} option.  (The
knots are placed on values of ln(t), not t.)  Denote the predicted values of
the log cumulative-hazard function by Z_hat(t).  The density function, f(t),
is

{pin}
f(t) = -dS(t)/dt = dS/dZ_hat dZ_hat/dt = S(t) exp(Z_hat) dZ_hat(t)/dt

{pstd}
dZ_hat(t)/dt is computed from the regression coefficients of the fitted spline
function.  The estimated survival function is calculated as

{pin}
S_hat(t) = exp{-exp Z_hat(t)}

{pstd}
The hazard function is calculated as f(t)/S_hat(t).

{pstd}
If {it:varlist} is specified, the baseline survival function (that is, at zero
values of the covariates) is used instead of the survival function of the raw
observations.  With {cmd:df(1)}, a Weibull model is fit.

{pstd}
With {cmd:scale(normal)}, smoothing is of the normal quantile function,
invnorm{1 - S(t)}, instead of the log cumulative-hazard function.  With
{cmd:df(1)}, a lognormal model is fit.

{pstd}
With {cmd:scale(odds)}, smoothing is of the log odds-of-failure function,
ln[{1 - S(t)}/S(t)], instead of the log cumulative-hazard function.  With
{cmd:df(1)}, a loglogistic model is fit.

{pstd}
Estimation is performed by maximum likelihood.  Optimization uses the default
technique, {cmd:nr} (Stata's version of Newton-Raphson iteration).


{title:Examples}

{pstd}Setup{p_end}
{phang2}{cmd:.} {bf:{stata "webuse brcancer"}}{p_end}
{phang2}{cmd:.} {bf:{stata "stset rectime, failure(censrec = 1)"}}{p_end}

{pstd}Proportional hazards model{p_end}
{phang2}{cmd:.} {bf:{stata "stpm2 hormon, scale(hazard) df(4) eform"}}{p_end}

{pstd}Proportional odds model{p_end}
{phang2}{cmd:.} {bf:{stata "stpm2 hormon, scale(odds) df(4) eform"}}{p_end}

{pstd}Time-dependent effects on a cumulative hazard scale{p_end}
{phang2}{cmd:.} {bf:{stata "stpm2 hormon, scale(hazard) df(4) tvc(hormon) dftvc(3)"}}{p_end}

{pstd}User-defined knots at centiles of uncensored event times{p_end}
{phang2}{cmd:.} {bf:{stata "stpm2 hormon, scale(hazard)  knots(20 50 80) knscale(centile)"}}{p_end}


{title:Authors}

{pstd}
Paul Lambert{p_end}
{pstd}University of Leicester{p_end}
{pstd}Leicester, UK{p_end}
{pstd}{browse "mailto:paul.lambert@leicester.ac.uk":paul.lambert@leicester.ac.uk}{p_end}

{pstd}The option-to-fit cure models were implemented by:{p_end}
{pstd}Therese Andersson{p_end}
{pstd}Karolinska Institutet{p_end}
{pstd}Stockholm, Sweden{p_end}
{pstd}{browse "mailto:therese.m-l.andersson@ki.se":therese.m-l.andersson@ki.se}{p_end}

{pstd}Various other additions and suggestions by:{p_end}
{pstd}Patrick Royston{p_end}
{pstd}MRC Clinical Trials Unit{p_end}
{pstd}London, UK{p_end}
{pstd}{browse "mailto:pr@ctu.mrc.ac.uk":pr@ctu.mrc.ac.uk}{p_end}


{title:Also see}

{p 4 14 2}Article:  {it:Stata Journal}, volume 12, number 4: {browse "http://www.stata-journal.com/article.html?article=st0165_1":st0165_1},{break}
                    {it:Stata Journal}, volume 9, number 2: {browse "http://www.stata-journal.com/article.html?article=st0165":st0165}

{p 7 14 2}Help:  {helpb stpm2 postestimation}, {manhelp stset ST}, {helpb stpm} (if installed)
{p_end}
