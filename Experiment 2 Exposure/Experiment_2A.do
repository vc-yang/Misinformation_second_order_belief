clear * 
insheet using "Experiment_2A_data.csv"

// pre-treatment filtering
sort id
gen idmatch=id==id[_n-1]
 drop if idmatch==1
drop if screener1_nicks~=5

// very little dropout, not sig diff across treatments
tabulate finished condition, chi2
 
 reshape long _belief1_ _belief2_ _interesting_, i(id) j(item_num)

gen plausibility_bin=10-floor((item_num-1)/4)

sum _belief1_ _belief2_

gen repeated = _interesting_~=.
gen type=condition==2

gen repeatedXtype=repeated*type

rename _belief1_ belief
replace belief=_belief2_ if _belief2_~=.  

label define bla 0 1st 1 2nd
label values type bla

// define what is true and false
gen real=item_num<=20
gen repeatedXreal=repeated*real
gen repeatedXrealXtype=repeated*real*type


label define blaY 0 "Novel" 1 Repeated
label values repeated blaY

label define bla1 0 False 1 True
label values real bla1



preserve
collapse (mean) belief if repeated==0 & type==0, by(item_num)
rename belief plausibility
save plaus, replace
restore

merge m:1 item_num using plaus
gen plausibilityR=round(plausibility*5)/5
gen plausibility2=plausibility^2
gen plausibility_bin2=plausibility_bin^2


preserve
collapse (mean) belief if type==0, by(item_num repeated)
reshape wide belief, i(item_num) j(repeated)
gen avg_plausibility=(belief0+belief1)/2
gen ITE=belief1-belief0
twoway (scatter  ITE avg_plausibility) (qfit ITE avg_plausibility)


save plaus2, replace
restore

drop _merge
merge m:1 item_num using plaus2


encode id, generate(idN)
xtset idN
//////////////////
// observational part

// 2nd order beliefs higher than first order beliefs
xi: reg belief i.item_num type if repeated==0, cluster(id )
// diff betwen 1st and 2nd order beliefs varies quadratically
xi: reg belief type##c.plausibility type##c.plausibility2 if repeated==0, cluster(id )
//graph bar (mean) belief if repeated==0, over(type) over(plausibilityR) b1title(Plausibility) name(obs_bar)
preserve
xi: statsby _b _se, clear by(item_num plausibility): reg belief type if repeated==0, cluster(id )
//twoway (scatter plausibility _b_type) (qfit plausibility _b_type) 
twoway (qfitci _b_type plausibility ) (scatter _b_type plausibility, xtitle(1st Order Belief) ytitle(2nd order - 1st order) name(obs))  
restore


*/ get the regression coefficients for centered dummy 
gen realC=real-0.5
gen typeC=type-0.5


gen beliefC=belief-0.5
gen repeatedC = repeated-0.5

ivreg2 beliefC c.repeatedC##c.typeC##c.realC, cluster( id item_num)


* --- Effect plot: (Repeated − Novel) distinguished by 1st/2nd order ---
* Same model
reg belief i.repeated##i.type##i.real, vce(cluster id)

quietly tab id if e(sample)
display "Number of participants: " r(r)

* Average marginal effect of repeated within each type (collapses over real)
margins, dydx(repeated) over(type)

marginsplot, ///
    recast(scatter) recastci(rcap) ///
    horizontal ///
    xline(0, lpattern(dash)) ///
    plotopts(msymbol(O) msize(medium)) ///
    ciopts(lwidth(medthin)) ///
    xtitle("Exposure effect") ytitle("") ///
	aspectratio(0.1) ///
	title("Trivia")

	
*--- Compute Bayes Factor for the null result using the BIC approximation for large N --- 
* FULL model
ivreg2 beliefC c.repeatedC##c.typeC##c.realC, cluster(id item_num)
est store FULL
gen byte __S = e(sample)             // lock the sample

/* count unique participants  */
egen byte __tag_id = tag(id) if __S  // marks one observation per participant within the estimation sample
quietly summarize __tag_id if __tag_id==1
scalar N_part = r(N)                 // number of participants (clusters)

scalar k_F = e(df_m) + 1             // add 1 for constant
scalar RSS_F = e(rss)
scalar BIC_F = N_part*ln(RSS_F/N_part) + k_F*ln(N_part)

* Reduced model a: effect of type on belief is the same, aggregaing over real, combining true and false. (drop repeatedC#typeC and 3-way repeatedC#typeC#realC terms)
ivreg2 beliefC ///
    c.repeatedC c.typeC c.realC ///
    c.repeatedC#c.realC ///
    c.typeC#c.realC, cluster(id item_num)
	
est store REDUCEDa
scalar k_R   = e(df_m) + 1
scalar RSS_R = e(rss)
scalar BIC_Ra = N_part*ln(RSS_R/N_part) + k_R*ln(N_part)

* Bayes factor (BIC approximation)
scalar dBICa = BIC_F - BIC_Ra
scalar BF01a = exp(dBICa/2)

* Reduced model b: effect of type on belief is the same, for either true or false in real.  (drop the 3 way interaction)
ivreg2 beliefC ///
    c.repeatedC c.typeC c.realC ///
    c.repeatedC#c.realC  c.typeC#c.realC ///
    c.repeatedC#c.typeC, cluster(id item_num)


est store REDUCEDb

scalar k_R   = e(df_m) + 1
scalar RSS_R = e(rss)
scalar BIC_Rb = N_part*ln(RSS_R/N_part) + k_R*ln(N_part)


* Bayes factor (BIC approximation)
scalar dBICb = BIC_F - BIC_Rb
scalar BF01b = exp(dBICb/2)

* Print results
di as txt "N participants = " %9.0f N_part
di as txt "k_full=" %6.0f k_F "  k_reduced=" %6.0f k_R
di as txt "BIC_full   = " %10.3f BIC_F
di as txt "BIC_reduced_a= " %10.3f BIC_Ra
di as txt "BIC_reduced_b= " %10.3f BIC_Rb

di as txt "ΔBIC_a = BIC_full - BIC_reduced_a = " %9.3f dBICa
di as txt "ΔBIC_b = BIC_full - BIC_reduced_b = " %9.3f dBICb

di as res "BF01_a aggregaing over true and false for real ≈ " %9.3f BF01a
di as res "BF01_b for either true or false in real ≈ " %9.3f BF01b

