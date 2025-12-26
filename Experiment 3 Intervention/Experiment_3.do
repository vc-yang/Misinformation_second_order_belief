clear * 
/*
// Below are for pre-processing data, which is saved into data.dta. Following analysis calls cleaned data, "data.dta" for analysis. 


insheet using "Experiment_3_data.csv"

drop if status==1

gen id=_n

reshape long _belief1_ _belief2_, i(id) j(item_num)

rename _belief1_ belief1
rename _belief2_ belief2

gen false=item_num<51
gen misleading=item_num>50 & item_num<61

gen item_type="True"
replace item_type="False" if false==1
replace item_type="Misleading" if misleading==1

label define bla 0 Control 3 Norms 6 Tips
label values condition bla




gen id2=_n
reshape long belief, i(id2) j(beliefType)



gen norms=condition==3
gen tips=condition==6

gen treated=condition>0

gen beliefTypeX=beliefType

replace beliefType=beliefType-1.5
zscore false misleading 

label define blaX 1 FirstOrder 2 2ndOrder 
label values beliefTypeX blaX

//ivreg2 belief norms##beliefType##false##norms false##tips  misleading##norms misleading##tips, cluster( id item_num)
replace belief=(belief-1)/5


save data, replace
*/
use data

*/ bar plot of means

*/ graph bar (mean) belief, over(condition) over(item_type) over(beliefTypeX)
 
*=== Bars with 95% CI error bars ===*
* Encode string variables first
encode item_type, gen(item_type_num)

cibar belief, over1(condition) over2(item_type_num) over3(beliefTypeX) ///
    bargap(10) level(95)




gen beliefType_1st=0.5-beliefType
gen beliefType_2nd=0.5+beliefType


// Regression analysis
ivreg2 belief norms##tips##c.beliefType##false  norms##tips##c.beliefType##misleading, cluster( id item_num)


*--- Compute Bayes Factor for the null result using the BIC approximation for large N --- 
* -- false headlines -- 
* FULL model for false headlines: 
ivreg2 belief norms##tips##c.beliefType##false, cluster( id item_num)
est store FULL
gen byte __S = e(sample)             // lock the sample

/* count unique participants  */
egen byte __tag_id = tag(id) if __S  // marks one observation per participant within the estimation sample
quietly summarize __tag_id if __tag_id==1
scalar N_part = r(N)                 // number of participants (clusters)

scalar N   = e(N)
scalar k_F = e(df_m) + 1             // add 1 for constant
scalar RSS_F = e(rss)
scalar BIC_1_false = N_part*ln(RSS_F/N_part) + k_F*ln(N_part)

* Reduced model: drop terms with the 3 way interaction tips × false x beliefType
ivreg2 belief ///
    norms tips c.beliefType false ///
    norms#c.beliefType norms#false ///
    tips#c.beliefType tips#false false#c.beliefType ///
    norms#c.beliefType#false ///
    , cluster(id item_num)
	
est store REDUCED_false_tip
scalar k_R   = e(df_m) + 1
scalar RSS_R = e(rss)
scalar BIC_0_false_tip = N_part*ln(RSS_R/N_part) + k_R*ln(N_part)

* Reduced model: drop terms with the 3 way interaction norms × false x beliefType
ivreg2 belief ///
    norms tips c.beliefType false ///
    norms#c.beliefType norms#false ///
    tips#c.beliefType tips#false false#c.beliefType ///
    tips#c.beliefType#false ///
    , cluster(id item_num)
	
est store REDUCED_false_norm
scalar k_R   = e(df_m) + 1
scalar RSS_R = e(rss)
scalar BIC_0_false_norm = N_part*ln(RSS_R/N_part) + k_R*ln(N_part)

* Bayes factor (BIC approximation)
scalar BF01_false_tip = exp((BIC_1_false - BIC_0_false_tip)/2)
scalar BF01_false_norm = exp((BIC_1_false - BIC_0_false_norm)/2)

* Print results
di as res "BF01_false headlines_tips × false x beliefType ≈ " %9.3f BF01_false_tip
di as res "BF01_false headlines_norms × false x beliefType ≈ " %9.3f BF01_false_norm


* --- Misleading headlines--- 
* FULL model for misleading headlines: 
ivreg2 belief norms##tips##c.beliefType##misleading, cluster( id item_num)
est store FULL
gen byte __S2 = e(sample)             // lock the sample

/* count unique participants  */
egen byte __tag_id2 = tag(id) if __S2  // marks one observation per participant within the estimation sample
quietly summarize __tag_id2 if __tag_id2 ==1
scalar N_part = r(N)                 // number of participants (clusters)

scalar k_F = e(df_m) + 1             // add 1 for constant
scalar RSS_F = e(rss)
scalar BIC_1_mislead = N_part*ln(RSS_F/N_part) + k_F*ln(N_part)

* Reduced model: drop terms with the 3 way interaction tips × misleading x beliefType
ivreg2 belief ///
    norms tips c.beliefType misleading ///
    norms#c.beliefType norms#misleading ///
    tips#c.beliefType tips#misleading misleading#c.beliefType ///
    norms#c.beliefType#misleading ///
    , cluster(id item_num)

est store REDUCED_false
scalar k_R   = e(df_m) + 1
scalar RSS_R = e(rss)
scalar BIC_0_mislead_tips = N_part*ln(RSS_R/N_part) + k_R*ln(N_part)

* Bayes factor (BIC approximation)
scalar BF01_mislead_tips = exp((BIC_1_mislead - BIC_0_mislead_tips)/2)

* Reduced model: drop terms with the 3 way interaction norms × misleading x beliefType
ivreg2 belief ///
    norms tips c.beliefType misleading ///
    norms#c.beliefType norms#misleading ///
    tips#c.beliefType tips#misleading misleading#c.beliefType ///
    tips#c.beliefType#misleading ///
    , cluster(id item_num)

est store REDUCED_false
scalar k_R   = e(df_m) + 1
scalar RSS_R = e(rss)
scalar BIC_0_mislead_norms = N_part*ln(RSS_R/N_part) + k_R*ln(N_part)

* Bayes factor (BIC approximation)
scalar BF01_mislead_norms = exp((BIC_1_mislead - BIC_0_mislead_norms)/2)


* Print results
di as res "BF01_misleading headlines, tips × misleading x beliefType ≈ " %9.3f BF01_mislead_tips
di as res "BF01_misleading headlines, norms × misleading x beliefType ≈ " %9.3f BF01_mislead_norms



