# addVariables
Package installation:
(install.packages("devtools"))
library("devtools")
install_github("dkurtoic/addVariables")
library("addVariables")

This package is introduced to help synchronize data filtering on PERMIT project. In this package you can find these functions:
addBMI - This function takes your dataset (crea.rep or similar) and aligns BMI event.dates with creatinine measurement event.dates.
It will add a BMI column with NAs at the beginning and then fill it with corresponding BMI values while running.
It loops through each patient ID, finds all event.dates and CodeValues for BMI in this patient, finds the closest event.date
to creatinine event.date (not necessarily after the creatinine event), and saves this BMI value. It is saving a new dataframe row by row, in every loop.
addDiabetes - This function takes your dataset (crea.rep or similar) and finds diabetes information on each patient in crea dataset. For each patient
it first finds the earliest event.date of diabetes and then searches for the closest creatinine event.date (making sure that crea event.date
was chronologically *after* the diabetes date)
addBP- This function takes your dataset (crea.rep or similar) and finds blood pressure information on each patient in crea dataset.
It is up to you which ReadCode you want to analyse. You can only analyse full ReadCode (not the ones that only begin with some string).
This function also adds Flag column to the final data set. If the BP measure is more than 30 days old, the Flag column for this patient's event.date will have TRUE.

For each of these functions you therefore need a dataset with creatinine CodeValues and "event.date" column for date when the creatinine values were measured. Patient ID 
column name is "PatientID".
You also need corresponding datasets for other variables, e.g. BMI, BP or diabetes. Datasets for this analysis were gotten from SIR data set, using specific
codes to filter the variables. For BMI  all code begining with "22K" were extracted. The data was then cleaned removing duplicated rows and NAs.
BP dataset was also filtered out from SIR data, codes extracted were "246A." for diastolic and "2469." for systolic. 
For diabetes, codes begining with "C10" or "66A" were extracted.

Other potentially useful functions in the package are:
variableDataEditing - the files gotten from add()* function are without column names. variableDataEditing will append column names and add 
column "formerge" that will help merging all the files together.
mergeMe - function for quite pretty merging
