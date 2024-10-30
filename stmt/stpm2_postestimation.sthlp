{smcl}
{* *! version 1.4.4 10Mar2012}{...}
{cmd:help stpm2 postestimation}{right: ({browse "http://www.stata-journal.com/article.html?article=st0165_1":SJ12-4: st0165_1})}
{hline}

{title:Title}

{p2colset 5 29 31 2}{...}
{p2col :{hi:stpm2 postestimation} {hline 2}}Postestimation tools for stpm2{p_end}
{p2colreset}{...}


{title:Description}

{pstd}
The following standard postestimation commands are available after
{cmd:stpm2}:

{synoptset 13}{...}
{p2coldent :Command}Description{p_end}
{synoptline}
INCLUDE help post_adjust2
INCLUDE help post_estat
INCLUDE help post_estimates
INCLUDE help post_lincom
INCLUDE help post_lrtest
INCLUDE help post_nlcom
{p2col :{helpb stpm2 postestimation##predict:predict}}predictions, residuals, influence statistics, and other diagnostic measures{p_end}
INCLUDE help post_predictnl
INCLUDE help post_test
INCLUDE help post_testnl
{synoptline}
{p2colreset}{...}


{marker predict}{...}
{title:Syntax for predict}

{p 8 16 2}
{cmd:predict} {newvar} {ifin} [{cmd:,} {it:statistic}]

{synoptset 30}{...}
{synopthdr :statistic}
{synoptline}
{synopt :{opt ab:c}}area between log hazard-ratio curves{p_end}
{synopt :{cmd:at(}{varname} {it:#} [...]{cmd:)}}predict at values of specified covariates{p_end}
{synopt :{opt cen:tile(#)}}request {it:#}th centile of survival distribution{p_end}
{synopt :{opt ci}}confidence interval{p_end}
{synopt :{opt cumh:azard}}cumulative hazard{p_end}
{synopt :{opt cumo:dds}}cumulative odds{p_end}
{synopt :{opt c:ure}}cure proportion{p_end}
{synopt :{opt dens:ity}}density function{p_end}
{synopt :{opt fail:ure}}failure function{p_end}
{synopt :{opt h:azard}}hazard function{p_end}
{synopt :{cmdab:hrn:umerator(}{varname} {it:#} [...]{cmd:)}}numerator for (time-dependent) hazard ratio{p_end}
{synopt :{cmdab:hrd:enominator(}{varname} {it:#} [...]{cmd:)}}denominator for
(time-dependent) hazard ratio{p_end}
{synopt :{cmd:hdiff1(}{varname} {it:#} [...]{cmd:)}}first hazard function for
difference in hazard functions{p_end}
{synopt :{cmd:hdiff2(}{varname} {it:#} [...]{cmd:)}}second hazard function for
difference in hazard functions{p_end}
{synopt :{opt mart:ingale}}martingale residuals{p_end}
{synopt :{opt means:urv}}population-averaged survival function{p_end}
{synopt :{opt n(#)}}number of evaluation points{p_end}
{synopt :{opt nor:mal}}standard normal deviate of survival function{p_end}
{synopt :{opt per(#)}}express hazard rates (and differences) per # person-years{p_end}
{synopt :{opt rm:st}}restricted mean survival time{p_end}
{synopt :{opt rsd:st}}standard deviation of restricted survival time{p_end}
{synopt :{cmd:sdiff1(}{varname} {it:#} [...]{cmd:)}}first survival curve for difference in survival functions{p_end}
{synopt :{cmd:sdiff2(}{varname} {it:#} [...]{cmd:)}}second survival curve for difference in survival functions{p_end}
{synopt :{opt stdp}}standard error of predicted function{p_end}
{synopt :{opt s:urvival}}survival function{p_end}
{synopt :{opth time:var(varname)}}time variable used for predictions; default is {cmd:timevar(_t)}{p_end}
{synopt :{opt tma:x(#)}}set upper bound of time for {opt rmst} and {opt abc} options{p_end}
{synopt :{opt tmi:n(#)}}set lower bound of time for {opt rmst} and {opt abc} options{p_end}
{synopt :{opth tvc(varname)}}time-varying coefficient for {it:varname}{p_end}
{synopt :{opt unc:ured}}survival and hazard functions for the "uncured"{p_end}
{synopt :{opt xb}}linear predictor{p_end}
{synopt :{opt xbnob:aseline}}linear predictor, excluding the spline function{p_end}
{synopt :{opt zero:s}}set all covariates to 0 (baseline prediction){p_end}
{synopt :{opt centol(#)}}define tolerance level when estimating centile{p_end}
{synopt :{opt dev:iance}}deviance residuals{p_end}
{synopt :{opt dxb}}derivative of linear predictor{p_end}
{synopt :{opt lev:el(#)}}set confidence level; default is {cmd:level(95)}{p_end}
{synopt :{opt startunc(#)}}set starting value for Newton-Raphson algorithm for
estimating a centile of the survival distribution of "uncured"{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2} 
Statistics are available both in and out of sample; type
{cmd:predict} {it:...} {cmd:if e(sample)} {it:...} if wanted only for the
estimation sample.{p_end}
{p 4 6 2} 


{title:Options for predict}

{pstd}
If a relative survival model has been fit by use of the {cmd:bhazard()}
option, then survival refers to relative survival and hazard refers to excess
hazard.

{phang}
{opt abc} evaluates the area between a constant log hazard-ratio and a
time-dependent log hazard-ratio curve.  It integrates the difference between a
log hazard-ratio curve and a constant log hazard-ratio over the time range
between {opt tmin()} and {opt tmax()}.  The constant hazard-ratio is supplied
by the {opt hr0()} option.  The time-dependent log hazard-ratio curve is
determined according to {cmd:hrnumerator()}, which must therefore be
specified.  You may also specify {opt hrdenominator()}.  The {opt n()},
{opt at()}, and {opt zeros} options are valid with {opt abc}.

{phang}
{cmd:at(}{varname} {it:#} [{it:varname #} ...]{cmd:)} requests that the
covariates specified by {it:varname} be set to {it:#}.  For example,
{cmd:at(x1 1 x3 50)} would evaluate predictions at {cmd:x1} = 1 and {cmd:x3} =
50.  This is a useful way to obtain out-of-sample predictions.  Note that if
{opt at()} is used together with {opt zeros}, all covariates not listed in
{opt at()} are set to 0.  If {opt at()} is used without {opt zeros}, then all
covariates not listed in {opt at()} are set to their sample values.

{phang}
{opt centile(#)} requests the {it:#}th centile of the survival
time distribution, calculated using a Newton-Raphson algorithm.

{phang}
{opt ci} calculates a confidence interval for the requested statistic and
stores the confidence limits in {it:newvar}{cmd:_lci} and
{it:newvar}{cmd:_uci}.

{phang}
{opt cumhazard} predicts the cumulative hazard function.

{phang}
{opt cumodds} predicts the cumulative odds-of-failure function.

{phang}
{opt cure} predicts the cure proportion after fitting a cure model.

{phang}
{opt density} predicts the density function.

{phang}
{opt failure} predicts the failure function, that is, F(t) = 1 - S(t).

{phang}
{opt hazard} predicts the hazard function.

{phang}
{cmd:hrnumerator(}{varname} {it:#} [{it:varname #} ...]{cmd:)}
specifies the numerator of the (time-dependent) hazard ratio.  By
default, all covariates not specified using this option are set to 0.
Setting the remaining values of the covariates to 0 may not always be
sensible, particularly with models other than those on the cumulative
hazard scale or when more than one variable has a time-dependent effect.
If {it:#} is set to missing ({cmd:.}), then {it:varname} has the values
defined in the dataset.

{phang}
{cmd:hrdenominator(}{varname} {it:#} [{it:varname #} ...]{cmd:)}
specifies the denominator of the hazard ratio.  By default, all
covariates not specified using this option are set to 0.  See the
cautionary note in {opt hrnumerator()}.  If {it:#} is set to missing
({cmd:.}), then {it:varname} has the values defined in the dataset.

{phang}
{cmd:hdiff1(}{varname} {it:#} [{it:varname #} ...]{cmd:)} and
{cmd:hdiff2(}{it:varname #} [{it:varname #} ...]{cmd:)} predict the
difference in hazard functions, with the first hazard function defined
by the covariate values listed for {opt hdiff1()} and the second by
those listed for {opt hdiff2()}.  By default, covariates not specified
using either option are set to 0.  Setting the remaining values of the
covariates to 0 may not always be sensible.  If {it:#} is set to missing
({cmd:.}), then {it:varname} has the values defined in the dataset.

{pmore}Example: {cmd:hdiff1(hormon 1)} (without specifying
{cmd:hdiff2()}) computes the difference in predicted hazard functions at
{cmd:hormon} = 1 compared with {cmd:hormon} = 0.

{pmore}Example: {cmd:hdiff1(hormon 2) hdiff2(hormon 1)} computes the
difference in predicted hazard functions at {cmd:hormon} = 2 compared
with {cmd:hormon} = 1.

{pmore}Example: {cmd:hdiff1(hormon 2 age 50) hdiff2(hormon 1 age 30)}
computes the difference in predicted hazard functions at {cmd:hormon} =
2 and {cmd:age} = 50 compared with {cmd:hormon} = 1 and {cmd:age} = 30.

{phang}
{opt martingale} calculates martingale residuals.

{phang}
{opt meansurv} calculates the population-averaged survival curve.  This
differs from the predicted survival curve at the mean of all the covariates in
the model.  A predicted survival curve is obtained for each subject, and all
the survival curves in a population are averaged.  The process can be
computationally intensive.  It is recommended that the {opt timevar()} option
be used to reduce the number of survival times at which the survival curves
are averaged.  Combining {cmd:meansurv} with the {cmd:at()} option enables
adjusted survival curves to be estimated.

{phang}
{opt n(#)} defines the number of evaluation
points for integrating the estimated survival function(s) with respect
to time.  The larger {it:#} is, the more accurate is the estimated
restricted mean survival time, but the longer the calculation takes.
There is no gain by setting {it:#} above 5000.  The default is
{cmd:n(1000)}.  {cmd:n()} is available only with {opt rmst}.

{phang}
{opt normal} predicts the standard normal deviate of the survival function.

{phang}
{opt per(#)} expresses hazard rates and difference in hazard rates per {it:#}
person-years.

{phang}
{opt rmst} evaluates the mean or restricted mean survival time.
This is done by integrating the predicted survival curve from 0 to 
{opt tmax(#)}; see also the {opt n()} and {opt tmax()} options.  Note
that the {opt at()} and {opt zeros} options are valid with {opt rmst}.

{phang}
{opt rsdst} evaluates the standard deviation of the (restricted)
survival time.  For a single sample, the standard error of the
restricted mean survival time may be estimated by dividing the standard
deviation by the square root of the number of observations.  See also
the {opt rmst}, {opt n()}, and {cmd:tmax()} options.  Note that the 
{cmd:at()} and {opt zeros} options are valid with {opt rsdst}.

{phang}
{cmd:sdiff1(}{varname} {it:#} [{it:varname #} ...]{cmd:)} and
{cmd:sdiff2(}{it:varname #} [{it:varname #} ...]{cmd:)} predict the
difference in survival curves with the first survival curve defined by
the covariate values listed for {opt sdiff1()} and the second by those
listed for {opt sdiff2()}.  By default, covariates not specified using
either option are set to 0.  Setting the remaining values of the
covariates to 0 may not always be sensible.  If {it:#} is set to missing
({cmd:.}), then {it:varname} has the values defined in the dataset.

{pmore}Example: {cmd:sdiff1(hormon 1)} (without specifying
{cmd:sdiff2()}) computes the difference in predicted survival curves at
{cmd:hormon} = 1 compared with {cmd:hormon} = 0.

{pmore}Example: {cmd:sdiff1(hormon 2) sdiff2(hormon 1)} computes the
difference in predicted survival curves at {cmd:hormon} = 2 compared
with {cmd:hormon} = 1.

{pmore}Example: {cmd:sdiff1(hormon 2 age 50) sdiff2(hormon 1 age 30)}
computes the difference in predicted survival curves at {cmd:hormon} = 2
and {cmd:age} = 50 compared with {cmd:hormon} = 1 and {cmd:age} = 30.

{phang}
{opt stdp} calculates the standard error of prediction and stores
it in {newvar}{cmd:_se}.  {cmd:stdp} is available only with the {cmd:xb}
and {cmd:dxb} options.

{phang}
{opt survival} predicts the survival function.

{phang}
{opth timevar(varname)} defines the variable used as time in the
predictions.  The default is {cmd:timevar(_t)}.  This is useful for
large datasets where for plotting purposes, predictions are needed for
only 200 observations, for example.  Some caution should be taken when
using this option because predictions may be made at whatever covariate
values are in the first 200 rows of data.  This can be avoided by using
the {cmd:at()} option or the {cmd:zeros} option to define the covariate
patterns for which you require the predictions.

{phang}
{opt tmax(#)} defines the upper limit of time over which the integration of
the estimated survival function is to be conducted.  The default is
{cmd:tmax(0)}, meaning an upper limit as close to t = infinity as is
reasonable (in fact, using the estimated 99.999999th centile of the survival
distribution).  {cmd:tmax()} is available only with {opt rmst} or {opt abc}.

{phang}
{opt tmin(#)} defines the lower bound of time over which the integration of
the estimated survival function is to be conducted.  The default is
{cmd:tmin(-1)}, taken as 0 and meaning a lower bound of 0.  {cmd:tmin()} is
available only with {opt rmst} or {opt abc}.

{phang}
{opth tvc(varname)} stands for time-varying coefficient and
computes the estimated coefficient for {it:varname}, a covariate in
{cmd:stpm2}'s {it:varlist}.  If {it:varname} is time-fixed, then
{it:newvar} will be a constant.  If {it:varname} is included in the
{cmd:tvc()} option, then {it:newvar} will depend on {cmd:_t} and may be
interpreted as the time-varying effect of {it:varname} on the chosen
scale of the model (proportional hazards, proportional odds, or probit).
For example, in a hazard-scale model ({cmd:scale(hazard)}), {it:newvar}
multiplied by {it:varname} will be an estimate of the time-varying log
cumulative-hazard ratio for {it:varname} (compared with {it:varname} =
0) at every observed value of {it:varname}.  {it:newvar} alone will give
the log cumulative-hazard ratio for a one-unit change in {it:varname}.
Note that the time-varying log cumulative-hazard ratio for {it:varname}
will not be identical to the time-varying log hazard-ratio for
{it:varname}.

{phang}
{opt uncured} can be used after fitting a cure model.  It can be
used with the {cmd:survival}, {cmd:hazard}, and {opt centile()} options
to base predictions for the uncured group.

{phang}
{opt xb} predicts the linear predictor, including the spline function.

{phang}
{opt xbnobaseline} predicts the linear predictor, excluding the spline
function, that is, only the time-fixed part of the model.

{phang}
{opt zeros} sets all covariates to 0 (baseline prediction).  For example,
{cmd:predict s0, survival zeros} calculates the baseline survival function.

{phang}
{opt centol(#)} defines the tolerance when searching for the predicted
survival time at a given centile of the survival distribution.  The default is
{cmd:centol(0.0001)}.

{phang}
{opt deviance} calculates deviance residuals.

{phang}
{opt dxb} calculates the derivative of the linear predictor.

{phang}
{opt level(#)} specifies the confidence level, as a percentage, for confidence
intervals.  The default is {cmd:level(95)} or as set by {helpb set level}.

{phang}
{opt startunc(#)} sets the starting value for the Newton-Raphson
algorithm for estimating a centile of the survival time distribution of
the uncured; the default is the 12.5th centile of the observed follow-up
times.


{title:Examples}

{pstd}Setup{p_end}
{phang2}{cmd:.} {bf:{stata "webuse brcancer"}}{p_end}
{phang2}{cmd:.} {bf:{stata "stset rectime, failure(censrec = 1)"}}{p_end}

{pstd}Proportional hazards model{p_end}
{phang2}{cmd:.} {bf:{stata "stpm2 hormon, scale(hazard) df(4) eform"}}{p_end}
{phang2}{cmd:.} {bf:{stata "predict h, hazard ci"}}{p_end}
{phang2}{cmd:.} {bf:{stata "predict s, survival ci"}}{p_end}

{pstd}Time-dependent effects on cumulative hazard scale{p_end}
{phang2}{cmd:.} {bf:{stata "stpm2 hormon, scale(hazard) df(4) tvc(hormon) dftvc(3)"}}{p_end}
{phang2}{cmd:.} {bf:{stata "predict hr, hrnumerator(hormon 1) ci"}}{p_end}
{phang2}{cmd:.} {bf:{stata "predict survdiff, sdiff1(hormon 1) ci"}}{p_end}
{phang2}{cmd:.} {bf:{stata "predict hazarddiff, hdiff1(hormon 1) ci"}}{p_end}

{pstd}Use of the {cmd:at()} option{p_end}
{phang2}{cmd:.} {bf:{stata "stpm2 hormon x1, scale(hazard) df(4) tvc(hormon) dftvc(3)"}}{p_end}
{phang2}{cmd:.} {bf:{stata "predict s60h1, survival at(hormon 1 x1 60) ci"}}{p_end}


{title:Authors}

{pstd}Paul Lambert{p_end}
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

{p 7 14 2}Help:  {helpb stpm2} (if installed)
{p_end}
