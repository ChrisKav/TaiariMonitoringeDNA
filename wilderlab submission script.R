# author: Christopher Kavazos
# date: 7/07/2022 
# aim: to transform data from collector app (http://libs.mobi/t/6018069850750976) to format required to upload to Wilderlab 

#Step 1: install relevant packages
library(readxl)
library(openxlsx)
library(splitstackshape)
library(data.table)
library(lubridate)
library(tidyverse)

#Step 2: Read relevant data sheet. Check getwd() or setwd() to correct filepath
data <- read_excel(".xlsx") #Inser directory path

#Step 3: Split UID column so that each UID has its own row
df <- cSplit(data, "UID", sep = ",", direction = "long")

#Step 4: Calculate deployment time
df$`Retrieval date` <- as.Date(df$`Retrieval date`, format = "%y/%m/%d")
df$`Deployment date` <- as.Date(df$`Deployment date`, format = "%y/%m/%d")

Retrieval = data.frame(
  date=df$`Retrieval date`,
  time=format(df$`Retrieval time`, "%H:%M")
)

Deployment = data.frame(
  date=df$`Deployment date`,
  time=format(df$`Deployment time`, "%H:%M")
)

Deployment <- as.POSIXct(paste(Deployment$date, Deployment$time), format="%Y-%m-%d %H:%M")
Retrieval <- as.POSIXct(paste(Retrieval$date, Retrieval$time), format="%Y-%m-%d %H:%M")

df$period <- round(difftime(Retrieval, Deployment, units = "hours"),1)

#Step 5: Format correctly for submission to wilderlab.co.nz/submit-samples
df2 <- cbind(df$UID, df$Site, df$Observer, df$`Deployment date`, df$Latitude, df$Longitude, df$`Volume filtered`,
             df$period, df$`Environment type`, df$notes)

colnames(df2) <- c("UID", "Reference", "Collector", "Date collected", "Latitude", "Longitude", "Volume filtered",
                   "Hours deployed", "Environment type", "Extra notes")

#Step 6: Save as an excel spreadsheet
write.xlsx(df2, "wilderlab_submission.xlsx")