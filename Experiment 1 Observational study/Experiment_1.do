clear *
insheet using "Experiment_1_data.csv"

drop if status==1
gen attentive=(screener_3_4==1 & screener_3_5==1 & screener_3_1==. & screener_3_2==. & screener_3_3==.)

gen demrep=demrep_c >3
replace demrep=. if  	demrep_c ==4 | demrep_c==.
gen demrep3=demrep
replace demrep3=3 if demrep_c==4
replace age=. if age<18
replace age=. if age>2222
replace age=2023-age if age>1000
gen lie_yes=i_lie_to_benefit<4
replace lie_yes=. if i_lie_to_benefit==.
gen college =education>=5
replace college=. if education==.

gen party3=party
replace party3=3 if party==4

xtile incomeM=income

label define bla 1 Democrat 2 Republican 3 Independent
label values party bla
label values party3 bla



gen id=_n
*/ count number of participants
count
display "Number of participants: " r(N)

reshape long _firstorder_ _secondorder_, i(id) j(item_num)

replace _firstorder_=(_firstorder_-1)/5
replace _secondorder_=(_secondorder_-1)/5

gen real=item_num>40


* regression to see the interaction of belief type and 
* Collapse first- and second-order beliefs into one variable
gen belief = .
replace belief = _firstorder_ if !missing(_firstorder_)
replace belief = _secondorder_ if !missing(_secondorder_)

* Indicator for belief type: 0 = first order, 1 = second order
gen belief_type = .
replace belief_type = 0 if !missing(_firstorder_)
replace belief_type = 1 if !missing(_secondorder_)

* Regression with interaction,
reg belief i.belief_type##i.real




* bar plot
* ===== Grouped CI bar chart with cibar: True/False on x, First vs Second within =====


    * Plot: x-axis groups = real (False/True), within-group bars = belief_type
*encode item_type, gen(item_type_num)

cibar belief, over1(belief_type) over2(real)  bargap(10) level(95) 
* ===== End grouped CI bar chart =====


* Observational t-test comparing mean of 1st order and 2nd order

collapse (mean) _firstorder_ _secondorder_, by(item_num )
gen real=item_num>40

gen two_m_one=_secondorder_ - _firstorder_
ttest _firstorder_ = _secondorder_
bysort real: ttest _firstorder_ = _secondorder_

pwcorr two_m_one _firstorder_ , sig




