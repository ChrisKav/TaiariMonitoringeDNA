library(readxl)
library(openxlsx)
library(splitstackshape)
library(data.table)
library(lubridate)
library(tidyverse)

data <- read_excel(".xlsx") #Inser directory path

df <- cSplit(data, "UID", sep = ",", direction = "long")

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

df <- df[which(df$UID %in% coa$UID),]

df2 <- cbind(df$UID, df$Site, df$Observer, df$`Deployment date`, df$Latitude, df$Longitude, df$`Volume filtered`,
             df$period, df$`Environment type`, df$notes)

colnames(df2) <- c("UID", "Reference", "Collector", "Date collected", "Latitude", "Longitude", "Volume filtered",
                   "Hours deployed", "Environment type", "Extra notes")

write.xlsx(df2, "wilderlab_submission.xlsx")