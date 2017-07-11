#' Adding BMI data to your final file
#'
#' This function takes your dataset (crea.rep or similar) and aligns BMI event.dates with creatinine measurement event.dates.
#' It will add a BMI column with NAs at the beginning and then fill it with corresponding BMI values while running.
#' It loops through each patient ID, finds all event.dates and CodeValues for BMI in this patient, finds the closest event.date
#' to creatinine event.date (not necessarily after the creatinine event), and saves this BMI value. It is saving a new dataframe row by row, in every loop.
#'
#'@details it uses parallel lapply - mclapply from "parallel" package. It's recommended to either run it in the background of your PC or on the cluster.
#'
#'@param crea.dataset dataset with creatinine ReadCodes, at least three creatinine measurements per patient.
#'
#'@param bmi.dataset file with extracted BMI codes from SIR data (codes starting with" 22K"); duplicates removed, NA free. All columns from SIR data can be present.
#'
#'@param NbOfCores how much cores do you want to use for your parallel process? Will be added to getOption() from mclapply. Default is 4 (should work on any PC), but can be changed to less or more.
#'
#'@param filename name of the final .txt file your output will be saved to. Can be a path also. Should be a txt file because of write.table function inside.
#'
#'@return this function returns nothig. It saves the output to your filename automatically
#'
#'@import parallel
#'
#'@importFrom ("install.packages", "installed.packages", "write.table")
#'
#' @examples
#'\dontrun{
#' mclapply_out <- addBMI(crea.dataset = crea.rep, bmi.dataset = bmiData, NbOfCores =  4L, filename = "~/mydoc/permit/creaWithBMI.txt")
#'}
#'
#' @export
addBMI <- function(crea.dataset, bmi.dataset, NbOfCores=4L, filename)
{
  pckgList<- "parallel"
  toInstall <- pckgList[!(pckgList %in% installed.packages()[,"Package"])]
  if(length(toInstall))
  {
    install.packages(toInstall)
    print("Package parallel is being installed")
  }
  NbOfCores <- as.integer(NbOfCores)
  crea.dataset$BMI <- NA
  crea.test <- crea.dataset
  bmiNoNA <- bmi.dataset

  mclapply(unique(crea.dataset$PatientID), function(x)
  {
    for (i in which(crea.test$PatientID == x)) #for each row of pat.ID X
    {
      if(sum(unique(bmiNoNA$PatientID) %in% x)==1)
      {
        BMIdate <- bmiNoNA[bmiNoNA$PatientID == x,"event.date"][which.min(abs(crea.test[i,"event.date"] - bmiNoNA[bmiNoNA$PatientID == x,"event.date"]))]
        BMIdate <- BMIdate[1]
        BMI <- min(bmiNoNA[bmiNoNA$event.date==BMIdate & bmiNoNA$PatientID==x,"CodeValue"])

        crea.test[i,"BMI"] <- BMI
        write.table(crea.test[i,], filename, row.names = F, col.names = F, append = T)
      }
      # if there is no BMI data for this Patient, NA is left in the BMI column
      else {write.table(crea.test[i,], filename, row.names = F, col.names = F, append = T)}
    }
  }, mc.cores = getOption("mc.cores",as.integer(NbOfCores)))
}

#' Adding diabetes data to your final file
#'
#' This function takes your dataset (crea.rep or similar) and finds diabetes information on each patient in crea dataset. For each patient
#' it first finds the earliest event.date of diabetes and then searches for the closest creatinine event.date (making sure that crea event.date
#' was chronologically *after* the diabetes date)
#'
#'@details it uses parallel lapply - mclapply from "parallel" package. It's recommended to either run it in the background of your PC or on the cluster.
#'
#'@param crea.dataset dataset with creatinine ReadCodes, at least three creatinine measurements per patient.
#'
#'@param diabetes.dataset file with extracted BMI codes from SIR data (codes starting with" 22K"); duplicates removed, NA free
#'
#'@param NbOfCores how much cores do you want to use for your parallel process? Will be added to getOption() from mclapply. Default is 4.
#'
#'@param filename name of the final .txt file your output will be saved to. Can be a path also. Should be a txt file because of write.table function inside.
#'
#'@return this function returns out file of mclapply function. Can be useful to save the output to a file in case an error occurrs. You can find error details in that output file. The creatinine data with diabetes is automatically saved, row by row.
#'
#'@import parallel
#'
#' @importFrom ("install.packages", "installed.packages", "write.table")
#'
#'@details if the parallel package is not installed already, it will be.
#'
#' @examples
#'\dontrun{
#' mclapply_out <- addBMI(crea.dataset = crea.rep, diabetes.dataset = diabetes, NbOfCores=4L, filename = "~/mydoc/permit/creaWithDiabetes.txt")
#'}
#'
#' @export
#'
addDiabetes <- function(crea.dataset, diabetes.dataset, NbOfCores=4L, filename)
{
  # install parallel if not installed already
  pckgList<- "parallel"
  toInstall <- pckgList[!(pckgList %in% installed.packages()[,"Package"])]
  if(length(toInstall))
  {
    install.packages(toInstall)
    print("Package parallel is being installed")
  }
  # prepare the data
  NbOfCores <- as.integer(NbOfCores)
  crea.test <- crea.dataset
  crea.test$Diabetes <- 0
  diabetes <- diabetes.dataset

  # run parallel lapply
  mclapply(unique(crea.test$PatientID), function(x)
  {
    if(sum(unique(diabetes$PatientID) %in% x)==1)
    {
      diabetesStartDate <- min(diabetes[diabetes$PatientID == x,"event.date"]) # find minimal diabetes event.date (diabetes Start Date)
      #subtract all the crea event.dates for this patient from diabetes StartDate to find the "switch" date, i.e. the crea event.date from which on the patient has diabetes
      alldiabdates <- crea.test[which(crea.test$PatientID==x), "event.date"]-diabetesStartDate
      names(alldiabdates) <- which(crea.test$PatientID == x) # to find the correct row number
      diabPos <- (alldiabdates[alldiabdates >= 0]) #take only positive differences, i.e. crea dates after the diabetesStartDate

      creaPos <- as.numeric(names(diabPos[which.min(diabPos)])) # take the crea event.date closest to the diabetesStartDate (min positive difference, event.date after diabetesStartDate)
      creaDiabStart <- crea.test[creaPos,"event.date"] # take the first crea date from which on the patient has diabetes
      crea.test[crea.test$event.date >= creaDiabStart & crea.test$PatientID==x,"Diabetes"] <- as.factor(1)
      #add 1 to those rows where the date is bigger or equal than creaDiabStart


      #calculate the time how long (in years) has patient been having diabetes (time.since.diagnosis)
      crea.test[crea.test$PatientID==x & crea.test$Diabetes==1, "time.since.diagnosis"] <- format(as.numeric((crea.test[crea.test$PatientID==x & crea.test$Diabetes==1,"event.date"] - diabetesStartDate)/365), digits=5)

      write.table(crea.test[crea.test$PatientID==x,], filename, row.names = F, col.names = F, append = T)
    }
    else #if there is no diabetes data on this Patient, just save without editing. They will have 0 in the Diabetes column.
    {
      write.table(crea.test[crea.test$PatientID==x,], filename, row.names = F, col.names = F, append = T)
    }
  }, mc.cores = getOption("mc.cores",NbOfCores))
}

#' Adding blood pressure  data to your final file
#'
#' This function takes your dataset (crea.rep or similar) and finds blood pressure information on each patient in crea dataset.
#' It is up to you which ReadCode you want to analyse. You can only analyse full ReadCode (not the ones that only begin with some string).
#' This function also adds Flag column to the final data set. If the BP measure is more than 30 days old, the Flag column for this patient's event.date will have TRUE.
#'
#'@details it uses lapply, not mclapply.
#'
#'@param crea.dataset dataset with creatinine ReadCodes, at least three creatinine measurements per patient.
#'
#'@param bpdata blood pressure data filtered out from SIR data or similar. It has to be clean data, no NA, no duplicates. Has to contain ReadCodes for blood pressure type, event.date of blood pressure measurement and patient ID.
#'
#'@param BPReadCode blood pressure ReadCode, character vector of length one (see examples).
#'
#'@param filename name of the final .txt file your output will be saved to. Can be a path also. Should be a txt file because of write.table function inside.
#'
#'@return this function returns out file of lapply function. Can be useful to save the output to a file in case an error occurrs. You can find error details in that output file. The creatinine data with diabetes is automatically saved, row by row.
#'
#'@import dplyr
#'
#'@importFrom ("install.packages", "installed.packages", "write.table")
#'
#'@details if the dplyr package is not installed already, it will be.
#'
#' @examples
#'\dontrun{
#' lapply_out <- addBMI(crea.dataset=crea.rep, bpdata, BPReadCode=c("246A."), filename="~/mydoc/permit/creaWithDiabetes.txt")
#'}
#'
#' @export
#'
addBP <- function(crea.dataset=crea.rep, bpdata=bpdata, BPReadCode, filename)
{

  #install dplyr if not already installed
  pckgList<- "dplyr"
  toInstall <- pckgList[!(pckgList %in% installed.packages()[,"Package"])]
  if(length(toInstall))
  {
    install.packages(toInstall)
    print("Package dplyr is being installed")
  }

  #prepare the data
  crea.test <- crea.rep
  crea.test$BP <- NA
  crea.test$Flag <- F

  #filter bp data. Analysing diastolic or systolic?
  bpdata_filtered <- bpdata %>% filter_(~ReadCode == ReadCode)
  bpdata_filtered$CodeValue <- as.numeric(as.character(bpdata_filtered$CodeValue)) # convert factor to numeric

  # run lapply
  lapply(unique(crea.test$PatientID), function(x)
  {

    for (i in which(crea.test$PatientID == x)) #for each row of pat.ID X
    {
      if(sum(unique(bpdata_filtered$PatientID) %in% x)==1) # whether pat. ID X has data on BP
      {
        # find the closest BP date to crea event.date
        BPdate <- bpdata_filtered[bpdata_filtered$PatientID == x,"event.date"][which.min(bpdata_filtered(dia[dia$PatientID == x,"event.date"]-crea.test[i,"event.date"]))]

        # if the clostest date is more than 30 days away from creatinine event.date add TRUE to Flag column of the final dataset
        # if not, the condition is skipped
        if (abs(crea.test[i, "event.date"]-BPdate) > 30)
        {
          crea.test[i,"Flag"] <- T
        }

        BPdate <- BPdate[1]
        #this is in case there are more than one BP measures for the same BPdate
        BP <- min(as.numeric(as.character((bpdata_filtered[bpdata_filtered$event.date==BPdate & bpdata_filtered$PatientID==x,"CodeValue"]))))
        #add in the BP
        crea.test[i,"BP"] <- BP
        write.table(crea.test[i,], filename, row.names = F, col.names = F, append = T)
      }
      else #if there is no BP data on pat.ID X, just save the row as it is (it already has NAs)
      {
        write.table(crea.test[i,], filename, row.names = F, col.names = F, append = T)
      }
    }
  })
}

#' Editing variable files from add*() functions
#'
#' This function will edit column names of the file gotten from add*() function. It will append column names from crea.rep dataset first, and then customly add the names
#' of the new column(s), depending on which file are you editing (BP, BMI or diabetes).
#' It will also add the "formerge" column where it creates a identifier for merging. This column is created by simply pasting PatientID, event.date and CodeValue of the variable file.
#'
#'@param crea.dataset dataset with creatinine ReadCodes, at least three creatinine measurements per patient.
#'
#'@param variable.data bmi, diabetes or blood pressure data gotten from one of the add*() functions in the package.
#'
#'@param BPtype a character string. E.g. "diastolic" or "systolic". How ever you type it, it will be converter to lower case. Default is NULL
#'
#'@param toedit A vector. Can be "BP", "BMI" or "diabetes". This will tell the function which file editing you want to do.
#'
#'@return it returns a edited BMI, diaebete, diastolic BP or systolic BP data frame.
#'
#'@details If toedit="BP", you have to specify BPtype ("diastolic" or "systolic"). This will create a edited BP dataframe - it renames last two columns to BP code value and "Flag". Regardless
#'of what you chose in the "toedit", function will also add the "formerge" column created for the purposes of merging all the variables to the final file.
#'I suggest saving the final files as RDS (saveRDS(editeddiastolicfile, "path/to/file.rds")), as it uses less space and loads faster.
#'
#' @examples
#'\dontrun{
#' editeddiastolicfile <- variableDataEditing(crea.rep, crea_dia.txt, BPtype="diastolic", toedit="BP")
#' editedBMIfile <- variableDataEditing(crea.rep, crea_dia.txt, toedit="BMI")
#'}
#'
#' @export
variableDataEditing <- function(crea.dataset = crea.rep, variable.data, BPtype=NULL, toedit=c("BP", "BMI", "diabetes"))

{
  if(toedit=="BP")
  {
    BPtype <- tolower(BPtype)
    BPtype_colnames <- paste0(BPtype, c("CodeValue", "Flag"))
    colnames(variable.data) <- colnames(crea.dataset)
    colnames(variable.data)[c(19,20)] <- BPtype_colnames
    variable.data$formerge <- paste(variable.data$PatientID, variable.data$event.date, variable.data$CodeValue, sep="_")
    return(variable.data)
  }

  if(toedit=="BMI")
  {
    colnames(variable.data) <- colnames(crea.dataset)
    colnames(variable.data)[grep("NA",colnames(variable.data))] <- "BMI"
    variable.data$formerge <- paste(variable.data$PatientID, variable.data$event.date, variable.data$CodeValue, sep="_")
    return(variable.data)
  }
  if(toedit=="diabetes")
  {
    colnames(variable.data) <- colnames(crea.dataset)
    colnames(variable.data)[grep("NA",colnames(variable.data))] <- c("Diabetes", "time.since.Diabetes.diagnosis")
    variable.data$formerge <- paste(variable.data$PatientID, variable.data$event.date, variable.data$CodeValue, sep="_")
    return(variable.data)
  }
}




