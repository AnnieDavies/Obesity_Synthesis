# Data
## Outcome data
The files 5_11_OutcomeData.xlsx and 12-18_OutcomeData.xlsx contain the outcome data for the analysis for the two age groups.

The columns in these files are defined as follows:
- study: identifying name of study
- CRCT: is the study a cluster randomized trial?
- CA: Does the study require SEs to be adjusted for clustering?
- measure: outcome scale on which data is reported (e.g. BMI)
- time: time point of measurment (short, mid, long)
- Aarm: Name of intervention arm
- Aint [not used]: Category of intervention arm (e.g. diet and activity)
- nA0: number of participants at baseline in intervention arm
- nA1: number of participants at follow-up in intervention arm
- nAc: number of clusters in intervention arm
- Barm, Bint, nB0, nB1, nBc: same as above for the reference arm
- calcICC: Did the trial report an ICC?
- ICC: value of ICC (if calcICC=Y)
- result: type of outcome data reported, either (i) baseline and follow-up (FU) means per arm, (ii) change from baseline per arm (CS), (iii) follow-up means per arm adjusted for baseline (Adjusted FU), (iv) follow-up means per arm (FU), (v) mean difference between arms (MD).
- AFmeanA, AFsdA: Adjusted FU mean & SD in intervention arm
- BFmeanA, BFsdA: Adjusted FU mean & SD in reference arm
- CSmeanA, CSsdA: Change score mean & SD in intervention arm
- CSmeanB, CSsdB: Change score mean & SD in intervention arm
- BFmean0A, BFsd0A: Baseline mean & SD in intervention arm
- BFmean0B, BFsd0B: Baseline mean & SD in reference arm
- BFmean1A, BFsd1A: Follow-up mean & SD in intervention arm
- BFmean1B, BFsd1B: Follow-up mean & SD in reference arm
- MDmean, MDse: Mean difference & SE between arms
- Map: Does this outcome require mapping onto zBMI?
- Map_ReqUnAdj: Does the data to be mapped need to be 'unadjusted' for clustering?
- Map_ReAdj: Once mapped, does the data need to be adjusted for clustering?
- Map_outcome: Outcome scale from which we are mapping
- Map_result: type of result data that is being mapped (as before e.g. baseline and FU etc...)
- MAP_AFmeanA, MAP_AFsdA, MAP_AFmeanB, MAP_AFsdB, MAP_CSmeanA, MAP_CSsdA, MAP_CSmeanB, MAP_CSsdB, MAP_BFmean0A, MAP_BFsd0A, MAP_BFmean0B, MAP_BFsd0B, MAP_BFmean1A, MAP_BFsd1A, MAP_BFmean1B, MAP_BFsd1B, MAP_MDmean, MAP_MDse: Equivalent data as above but this is the data to be mapped onto zBMI (e.g. with arm level baseline & follow-up extracted where available)
- Mean_Age, SD_Age: Mean and SD of age of participants in the trial
- Prop_Male: percentage of males in the trial
- fu_months: length of follow-up in months
- country: country in which the trial is conducted
- chart_ref: does the trial report an LMS reference chart?
- chart: which LMS reference chart to use

