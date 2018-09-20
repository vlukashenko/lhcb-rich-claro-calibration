##########################################################################################
################ Instructions to use S_curve.cpp and calibration.cpp #####################
########################## by Valeriia Lukashenko 2018 ###################################
##########################################################################################

I. S_curve.cpp
Used for: initial preselection of channels (empty/ON), plotting all s-curves, s-curve fit, transition point estimation
Options: th_value defines transition level, by default 0.5
         two fits: Erf and logistic sigmoid, by default Erf
Note: if not 0.5 is used transition point error should be calculated with errorXlogistic, errorXerf; by default if 0.5 calculated as parameter error.

Input: DAC scan root files; uses TGraphAsymmErrors 
Output: "threshold.txt".
"threshold.txt" constituents: 
    "The threshold level:
     Attenuation: 
     Offset: 
     Threshold: 
     Channel ID      Middle value, DAC st      Error DOWN, DAC st      Error UP, DAC st      Middle value, ke-      logistic, DAC st      Error DOWN, DAC st      Error UP, DAC st      Chi2      Fit status      logistic, ke      erf, DAC st      Error DOWN, DAC st      Error UP, DAC st      chi2      Fit status      erf, ke     
     EC0PMT1CH1      ..."

Attention! Don't forget to change data:    int offset[n]={}; //values of offsets
    					   int attenuation[n]={}; //values of attenuation
    					   int thresholds[n]={}; //values of threholds
    					   string names[n]={}; //names of data files
                               !Everything has to be in agreement!
Also, change dir value to switch from board to board or to point to the data directory.

II. calibration.cpp
Used for: building a calibration curve
Options: using different transition points calculations. if opt==1 - middle point calculation; if opt==2 - logistic sigmoid fit; if opt==3 - Erf fit.

Attention! there are two hacks for data sample: line 583 and line 566;

Input: "threshold.txt" from S_curve.cpp
Ouput: "calibration_new.root", "avaraged_results.txt", "file_par.txt".
"calibration_new.root" constituents:
	2 trees: calibration06, calibration13
	calibration06:
		EC06
		PMT06
		CH06

		thresholdstep_06_att0 //off:1 att:0
		thresholdstep_06_att1 //off:1 att:1 
		thresholdstep_06_att2 //off:1 att:2
		thresholdstep_06_att3 //off:1 att:3
		thresholdstep_06_off0 //off:0 att:1
	
		offset_06_att0 //off:1 att:0
		offset_06_att1 //off:1 att:0
		offset_06_att2 //off:1 att:0
		offset_06_att3 //off:1 att:0
		offset_06_off0 //off:0 att:1

	calibration13:
		EC13
		PMT13
		CH13

		thresholdstep_13_att0 //off:1 att:0
		thresholdstep_13_att1 //off:1 att:1 
		thresholdstep_13_att2 //off:1 att:2
		thresholdstep_13_att3 //off:1 att:3
		thresholdstep_13_off0 //off:0 att:1
	
		offset_13_att0 //off:1 att:0
		offset_13_att1 //off:1 att:0
		offset_13_att2 //off:1 att:0
		offset_13_att3 //off:1 att:0
		offset_13_off0 //off:0 att:1

"avaraged_results.txt" constituents:
		Name                   Attenuation                        Offset                      all, ke-
		slop                             3                             1                       ...
		offs                             3                             1                       ...

"file_par.txt" constituents:
   		"Name      Attenuation      Offset      slope, ke-      |at zero|, ke-
		EC0_PMT0_Ch1 ...

 





