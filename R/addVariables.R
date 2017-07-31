#' Adding BMI data to your final file
#'
#' This function takes your dataset (crea.rep or similar) and aligns BMI event.dates with creatinine measurement event.dates.
#' It will add a BMI column with NAs at the beginning and then fill it with corresponding BMI values while running.
#' It loops through each patient ID, finds all event.dates and CodeValues for BMI in this patient, finds the closest event.date
#' to creatinine event.date (not necessarily after the creatinine event), and saves this BMI value. It is saving a new dataframe row by row, in every loop.
#'
#'@details This function uses parallel lapply - mclapply from "parallel" package. It's recommended to either run it in the background of your PC or on the cluster. The column used for time analysis is named "event.date".
#'
#'@param crea.dataset dataset with creatinine ReadCodes, at least three creatinine measurements per patient.
#'
#'@param bmi.dataset file with extracted BMI codes from SIR data (codes starting with" 22K"); duplicates removed, NA free. All columns from SIR data can be present.
#'
#'@param NbOfCores how much cores do you want to use for your parallel process? Will be added to getOption() from mclapply. Default is 4 (should work on any PC), but can be changed to less or more.
#'
#'@param filename name of the final .txt file your output will be saved to. Can be a path also. Should be a txt file because of write.table function inside.
#'
#'@return This function returns output file of mclapply function. Can be useful to save the output to a file in case an error occurrs because you can find error details in the output file. The creatinine data with BMI is automatically saved to a .txt file, row by row, and can later be accessed by using "read.table" function.
#'
#'@import parallel
#'
#'@importFrom utils write.table
#'
#' @examples
#'\dontrun{
#'addBMI(crea.dataset, bmi.dataset, NbOfCores =  4L, filename = "creaWithBMI.txt")
#'running the script in your background:
#'nohup Rscript myscript.R &
#'}
#'
#' @export
addBMI <- function (crea.dataset, bmi.dataset, NbOfCores = 4L, filename)
{
  crea.dataset$BMI <- NA
  crea.dataset$event.date <- as.Date.character(crea.dataset$event.date)
  bmi.dataset$event.date <- as.Date.character(bmi.dataset$event.date)
  bmi.dataset$CodeValue <- as.numeric(as.character(bmi.dataset$CodeValue))

  mclapply(unique(crea.dataset$PatientID), function(x)
  {
    bmi.x <-  bmi.dataset[bmi.dataset$PatientID == x, c("event.date", "CodeValue")]
    crea.rows <- which(crea.dataset$PatientID == x)
    crea.x <- crea.dataset[crea.rows,"event.date"]
    for (i in 1:length(crea.x))
    {
      if (x %in% (bmi.dataset$PatientID))
      {
        BMIdate <- bmi.x[which.min(abs(crea.x[i] - bmi.x$event.date))[1],"event.date"]
        BMI <- min(bmi.x[bmi.x$event.date ==
                           BMIdate, "CodeValue"])
        crea.dataset[crea.rows[i], "BMI"] <<- BMI
        write.table(crea.dataset[i, ], filename, row.names = F, col.names = F, append = T)
      }
      else
      {
        write.table(crea.dataset[i, ], filename, row.names = F, col.names = F, append = T)
      }
    }
  }, mc.cores = getOption("mc.cores", as.integer(NbOfCores)))
}

#' Adding diabetes data to your final file
#'
#' This function takes your dataset (crea.rep or similar) and finds diabetes information on each patient in crea dataset. For each patient
#' it first finds the earliest event.date of diabetes and then searches for the closest creatinine event.date (making sure that crea event.date
#' was chronologically *after* the diabetes date)
#'
#'@details This function uses parallel lapply - mclapply from "parallel" package. It's recommended to either run it in the background of your PC or on the cluster. The column used for time analysis is named "event.date".
#'
#'@param crea.dataset dataset with creatinine ReadCodes, at least three creatinine measurements per patient.
#'
#'@param diabetes.dataset file with extracted BMI codes from SIR data (codes starting with" 22K"); duplicates removed, NA free
#'
#'@param NbOfCores how much cores do you want to use for your parallel process? Will be added to getOption() from mclapply. Default is 4.
#'
#'@param filename name of the final .txt file your output will be saved to. Can be a path also. Should be a txt file because of write.table function inside.
#'
#'@return This function returns out file of mclapply function. Can be useful to save the output to a file in case an error occurrs. You can find error details in that output file. The creatinine data with diabetes is automatically saved to a .txt file, row by row, and can later be accessed by using "read.table" function.
#'
#'@import parallel
#'
#'@importFrom utils write.table
#'
#' @examples
#'\dontrun{
#'addBMI(crea.dataset, diabetes.dataset = diabetes, NbOfCores=4L, filename = "creaWithDiabetes.txt")
#'}
#'
#' @export
addDiabetes <- function(crea.dataset, diabetes.dataset, NbOfCores=4L, filename)
{
  # prepare the data
  NbOfCores <- as.integer(NbOfCores)
  crea.dataset$Diabetes <- 0
  crea.dataset$diabetes.exposure <- 0

  diabetes$event.date <- as.Date(diabetes$event.date)

  # run parallel lapply
  mclapply(unique(crea.dataset$PatientID), function(x)
  {
    diab.x <- diabetes[diabetes$PatientID==x, "event.date"]
    crea.x <- crea.dataset[crea.dataset$PatientID==x,]
    crea.x$rows <- which(crea.dataset$PatientID==x)

    if(length(diab.x) != 0)
    {
      diabetesStartDate <- min(diab.x) # find minimal diabetes event.date (diabetes Start Date)
      #subtract all the crea event.dates for this patient from diabetes StartDate to find the "switch" date, i.e. the crea event.date from which on the patient has diabetes
      alldiffs <- crea.x[, "event.date"]-diabetesStartDate
      names(alldiffs) <- crea.x$rows # to find the correct row number

      diabPos <- (alldiffs[alldiffs >= 0]) #take only positive differences, i.e. crea dates after the diabetesStartDate

      creaPos <- as.numeric(names(diabPos[which.min(diabPos)])) # take the crea event.date closest to the diabetesStartDate (min positive difference, event.date after diabetesStartDate)
      creaDiabStart <- crea.dataset[creaPos,"event.date"] # take the first crea date from which on the patient has diabetes
      crea.x[crea.x$event.date >= creaDiabStart, "Diabetes"] <- as.factor(1)
      #add 1 to those rows where the date is bigger or equal than creaDiabStart


      #calculate the time how long (in years) has patient been having diabetes (diabetes.exposure)
      crea.x[crea.x$Diabetes==1, "diabetes.exposure"] <- format(as.numeric((crea.x[crea.x$Diabetes==1,"event.date"] - diabetesStartDate)/365), digits=5)

      write.table(crea.x, filename, row.names = F, col.names = F, append = T)
    }
    else #if there is no diabetes data on this Patient, just save without editing. They will have 0 in the Diabetes column.
    {
      write.table(crea.x, filename, row.names = F, col.names = F, append = T)
    }
  }, mc.cores = getOption("mc.cores", as.integer(NbOfCores)))
}

#' Adding blood pressure  data to your final file
#'
#' This function takes your dataset (crea.rep or similar) and finds blood pressure information on each patient in crea dataset.
#' It is up to you which ReadCode you want to analyse. You can only analyse full ReadCode (not the ones that only begin with some string).
#' This function also adds Flag column to the final data set. If the BP measure is more than 30 days old, the Flag column for this patient's event.date will have TRUE.
#'
#'@details This function uses lapply, not mclapply. Even more recommended to run it in the background or on cluster. The column used for time analysis is named "event.date".
#'
#'@param crea.datasett dataset with creatinine ReadCodes, at least three creatinine measurements per patient.
#'
#'@param bpdata blood pressure data filtered out from SIR data or similar. It has to be clean data, no NA, no duplicates. Has to contain ReadCodes for blood pressure type, event.date of blood pressure measurement and patient ID.
#'
#'@param BPReadCode blood pressure ReadCode, character vector of length one (see examples).
#'
#'@param NbOfCores how much cores do you want to use for your parallel process? Will be added to getOption() from mclapply. Default is 4.
#'
#'@param filename name of the final .txt file your output will be saved to. Can be a path also. Should be a txt file because of write.table function inside.
#'
#'@return This function returns out file of lapply function. Can be useful to save the output to a file in case an error occurrs. You can find error details in that output file. The creatinine data with diabetes is automatically saved to a .txt file, row by row, and can later be accessed by using "read.table" function.
#'
#'@import dplyr
#'
#'@import parallel
#'
#'@importFrom utils write.table
#'
#' @examples
#'\dontrun{
#'addBP(crea.dataset, bpdata, BPReadCode=c("246A."), filename="creaWithDiabetes.txt")
#'}
#'
#' @export
addBP <- function(crea.datasett, bpdata=bpdata, BPReadCode, NbOfCores=4L, filename)
{

  #prepare the data
  crea.dataset$BP <- NA
  crea.dataset$Flag <- F

  #filter bp data. Analysing diastolic or systolic?
  bpdata_filtered <- bpdata %>% filter_(~ReadCode == ReadCode)
  bpdata_filtered$CodeValue <- as.numeric(as.character(bpdata_filtered$CodeValue)) # convert factor to numeric

  # run lapply
  mclapply(unique(crea.dataset$PatientID), function(x)
  {
    crea.x <- crea.dataset[crea.dataset$PatientID==x,c("event.date", "BP", "Flag")]
    bp.x <- bpdata[bpdata$PatientID==x,c("event.date", "CodeValue")]

    for (i in 1:length(crea.x)) #for each row of pat.ID X
    {
      if(dim(bp.x)[1] != 0) # whether pat. ID X has data on BP
      {
        # find the closest BP date to crea event.date
        BPdate <- bp.x[,"event.date"][which.min(bp.x[,"event.date"]-crea.x[i, "event.date"])]

        # if the clostest date is more than 30 days away from creatinine event.date add TRUE to Flag column of the final dataset
        # if not, the condition is skipped
        if (abs(crea.x[i, "event.date"]-BPdate) > 30)
        {
          crea.x[i,"Flag"] <- T
        }

        BPdate <- BPdate[1]
        #this is in case there are more than one BP measures for the same BPdate
        BP <- min(as.numeric(as.character((bp.x[bp.x$event.date==BPdate,"CodeValue"]))))
        #add in the BP
        crea.x[i,"BP"] <- BP

        write.table(crea.x[i,], filename, row.names = F, col.names = F, append = T)
      }
      else #if there is no BP data on pat.ID X, just save the row as it is
      {
        write.table(crea.x[i,], filename, row.names = F, col.names = F, append = T)
      }
    }
  }, mc.cores=getOption("mc.cores", as.integer(NbOfCores)))
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
#'@param BPtype a character string. E.g. "diastolic" or "systolic", or simply "D" and "S". How ever you type it, it will be converter to upper case. Default is NULL
#'
#'@param toedit A vector. Can be "BP", "BMI" or "diabetes". This will tell the function which file editing you want to do.
#'
#'@return it returns a edited crea-BMI, crea-diabetes, crea-diastolic BP or crea-systolic BP data frame.
#'
#'@details If toedit="BP", you have to specify BPtype ("diastolic" or "systolic"). This will create a edited BP dataframe - it renames last two columns to BP code value and "Flag". Regardless
#'of what you chose in the "toedit", function will also add the "formerge" column created for the purposes of merging all the variables to the final file.
#'I suggest saving the final files as RDS (saveRDS(editeddiastolicfile, "path/to/file.rds")), as it uses less space and loads faster.
#'The column used for time analysis is named "event.time".
#'
#' @examples
#'\dontrun{
#' editeddiastolicfile <- variableDataEditing(crea.rep, crea_dia.txt, BPtype="diastolic", toedit="BP")
#' editedBMIfile <- variableDataEditing(crea.rep, crea_dia.txt, toedit="BMI")
#'}
#'
#' @export
variableDataEditing <- function(crea.dataset, variable.data = variable.data, BPtype=NULL, toedit=c("BP", "BMI", "diabetes"))

{
  if(toedit=="BP")
  {
    if(class(BPtype) != "character") {print("BPtype has to be a character vector.")}
    BPtype <- toupper(substr(BPtype, 1,1))
    BPtype_colnames <- paste0(BPtype, c("BP", "BP30d"))
    colnames(variable.data) <- colnames(crea.dataset)
    colnames(variable.data)[is.na(colnames(variable.data))] <- BPtype_colnames
    variable.data$formerge <- paste(variable.data$PatientID, variable.data$event.date, variable.data$CodeValue, sep="_")
    return(variable.data)
  }

  if(toedit=="BMI")
  {
    colnames(variable.data) <- colnames(crea.dataset)
    colnames(variable.data)[is.na(colnames(variable.data))] <- "BMI"
    variable.data$formerge <- paste(variable.data$PatientID, variable.data$event.date, variable.data$CodeValue, sep="_")
    return(variable.data)
  }
  if(toedit=="diabetes")
  {
    colnames(variable.data) <- colnames(crea.dataset)
    colnames(variable.data)[is.na(colnames(variable.data))] <- c("Diabetes", "diabetes.exposure")
    variable.data$formerge <- paste(variable.data$PatientID, variable.data$event.date, variable.data$CodeValue, sep="_")
    return(variable.data)
  }
}

#' Merging it all together
#'
#' Helps to merge four datasets (BMI, Systolic BP, Diastolic BP and Diabetes). First one you submit (df1) will be the "main", meaning that all the columns
#' of this data set will be preserved. In the second dataset only some columns will be preserved (according to your demands). You customize that
#' in the colsToSelect argument. "formerge" column is always preserved, and you specify the others.
#'
#' @param df1 one of the datasets from variableDataEditing() function.
#' @param df2 one of the datasets from variableDataEditing() function.
#' @param colsToSelect Default is NULL. A character vector of columns which you want to save from df2. No quotation marks. If you leave it as NULL, only formerge column will be added to the final data frame.
#'
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' mergeMe(bmidata, sysdata, colsToSelect=c(systolicCodeValue, systolicFlag))
#' mergeMe(bmisys, diastdata, colsToSelect=c(diastolicCodeValue, diastolicFlag))
#' mergeMe(bmisysdiast, diabetes.data, colsToSelect=c("Diabetes","diabetes.exposure"))
#' }
#'
#' @export

mergeMe <- function(df1, df2, colsToSelect=NULL)
{
  df2 <- df2 %>% select_(~formerge, ~colsToSelect)
  df1df2 <- merge(df1, df2, by="formerge")
  return(df1df2)
}




