options(stringsAsFactors = FALSE)

library("grofit")
library("doBy")
library("RColorBrewer")
library("ggplot2")
library("gdata")
library("plyr")
library("reshape")
library("plotrix")
library("corrplot")


setwd("C:/Users/Dana/OneDrive - UW-Madison/Y1000PhenotypingRawData/")

PlateInfo = read.csv("AllTreatmentInfo.csv", header = TRUE, check.names = FALSE)
StoragePlates = read.csv("AllSpeciesLayout.csv", header = TRUE, check.names = FALSE)
SpeciesLayout = read.csv("PlateInfo.csv", header = TRUE, check.names = FALSE)
FolderList = read.csv("ExperimentList.csv", header = TRUE, check.names = FALSE)

Directories = FolderList[,1]
Directories=Directories[-which(Directories %in% c("S2R2"))]
Directories2 = c("S2R2")

for(i in 1:length(Directories)){
  DIR = paste("C:/Users/Dana/OneDrive - UW-Madison/Y1000PhenotypingRawData/", Directories[[i]], sep = "")
  setwd(DIR)
  DataPlates = list.files()
  DataPlates=strsplit(DataPlates, "[.]") # Remove end file file name
  DataPlate_list = list()
 
   for(j in 1:length(DataPlates)){
    DataPlate_list[[j]] = DataPlates[[j]][1]
  }
  
  for(j in 1:length(DataPlate_list)){
    filename = paste(DataPlate_list[[j]]) # Creates an object name for each plate
    assign(filename, read.csv(paste(DataPlate_list[j], ".csv", sep = ""), header = TRUE, check.names = T, row.names = NULL, skip = 11)) # Reads in files from the directory based on the file names created above
  }
  
  EXPERIMENT = strsplit(Directories[[i]], split =  "")
  Set = rep(EXPERIMENT[[1]][2], times = 96)
  Rep = rep(EXPERIMENT[[1]][4], times = 96)
  
  for(j in 1:length(DataPlate_list)){
    filename = DataPlate_list[[j]]

    TEMP = strsplit(DataPlate_list[[j]], "_")[[1]][2]
    TEMP2 = strsplit(TEMP, split = "T")
    TEMP3 = strsplit(TEMP2[[1]][1], split = "P")[[1]][2]
    
    Treatment = rep(TEMP2[[1]][2], times = 96)
    Plate = rep(TEMP3, times = 96)
    
    DROPCOL = which(colnames(get(filename)) == "Content")
    assign(filename, get(filename)[-1,])
    assign(filename, get(filename)[,-DROPCOL])
    assign(filename, cbind(Set, get(filename)))
    assign(filename, cbind(Rep, get(filename)))
    assign(filename, cbind(Treatment, get(filename)))
    assign(filename, cbind(Plate, get(filename)))
  }
}

for(i in 1:length(Directories2)){
  DIR = paste("C:/Users/Dana/OneDrive - UW-Madison/Y1000PhenotypingRawData/", Directories2[[i]], sep = "")
  setwd(DIR)
  DataPlates = list.files()
  DataPlates=strsplit(DataPlates, "[.]") # Remove end file file name
  DataPlate_list = list()
  
  for(j in 1:length(DataPlates)){
    DataPlate_list[[j]] = DataPlates[[j]][1]
  }
  
  for(j in 1:length(DataPlate_list)){
    filename = paste(DataPlate_list[[j]]) # Creates an object name for each plate
    assign(filename, read.csv(paste(DataPlate_list[j], ".csv", sep = ""), header = TRUE, check.names = T, row.names = NULL, skip = 5)) # Reads in files from the directory based on the file names created above
  }
  
  EXPERIMENT = strsplit(Directories2[[i]], split =  "")
  Set = rep(EXPERIMENT[[1]][2], times = 96)
  Rep = rep(EXPERIMENT[[1]][4], times = 96)
  
  for(j in 1:length(DataPlate_list)){
    filename = DataPlate_list[[j]]
    
    TEMP = strsplit(DataPlate_list[[j]], "_")[[1]][2]
    TEMP2 = strsplit(TEMP, split = "T")
    TEMP3 = strsplit(TEMP2[[1]][1], split = "P")[[1]][2]
    
    Treatment = rep(TEMP2[[1]][2], times = 96)
    Plate = rep(TEMP3, times = 96)
    
    DROPCOL = which(colnames(get(filename)) == "Content")
    
    assign(filename, get(filename)[-1,])
    assign(filename, get(filename)[,-DROPCOL])
    
    if(nrow(get(filename)) > 96){
      assign(filename, get(filename)[-1,])
    }
    assign(filename, cbind(Set, get(filename)))
    assign(filename, cbind(Rep, get(filename)))
    assign(filename, cbind(Treatment, get(filename)))
    assign(filename, cbind(Plate, get(filename)))
  }
  
  
}

StoragePlates = StoragePlates[,-which(colnames(StoragePlates) %in% c("Row", "Col"))]
SpeciesInfo = merge(SpeciesLayout, StoragePlates, by = "Sp#")

SpeciesCount = data.frame(Species = character(), 
                          Count = numeric())
SPECIES = unique(SpeciesInfo$Species)

for(i in 1:length(SPECIES)){
  SpeciesCount[i,1] = SPECIES[[i]]
  SpeciesCount[i,2] = length(which(SpeciesInfo$Species == SPECIES[[i]]))
}

AllInfo = merge(SpeciesInfo, PlateInfo, by = c("Row", "Column", "Set", "Rep"))
setwd("C:/Users/Dana/OneDrive - UW-Madison/Y1000PhenotypingRawData/")
write.csv(AllInfo, "MergedPlateInfo.csv")

PlateNames = read.csv("AllPlates.csv", header = TRUE)
AllPlates = PlateNames[,1]

FixPlates = read.csv("FixPlates.csv", header = TRUE)
WellCol = read.csv("WellCol.csv", header = TRUE)

FixPlates = FixPlates[,1]

for(i in 1:length(FixPlates)){
  filename = FixPlates[[i]]
  assign(filename, get(filename)[,-which(colnames(get(filename)) == "Well")])
  assign(filename, cbind(WellCol, get(filename)))
}

AllPlates_list = list()
for(i in 1:length(AllPlates)){
  AllPlates_list[[i]] = list()
  filename = AllPlates[[i]]
  SET = as.numeric(unique(get(filename)$Set))
  REP = as.numeric(unique(get(filename)$Rep))
  TREAT = as.numeric(unique(get(filename)$Treatment))
  PLATE = as.numeric(unique(get(filename)$Plate))
  if(length(which(colnames(get(filename)) == "Well.Row") )< 1){
    TEMP = AllInfo[which(AllInfo$Set == SET & AllInfo$Rep == REP & AllInfo$Treatment == TREAT & AllInfo$Plate == PLATE),]
    assign(filename, merge(TEMP, get(filename), by.x = c("Row", "Column"), by.y = c("X", "X.1")))
    AllPlates_list[[i]] = get(filename)
    
  }else{
    TEMP = AllInfo[which(AllInfo$Set == SET & AllInfo$Rep == REP & AllInfo$Treatment == TREAT & AllInfo$Plate == PLATE),]
    assign(filename, merge(TEMP, get(filename), by.x = c("Row", "Column"), by.y = c("Well.Row", "Well.Col")))
    AllPlates_list[[i]] = get(filename)
  }
}

ncol_length = list()
for(i in 1:length(AllPlates)){
  filename = AllPlates[[i]]
  ncol_length[[i]] = ncol(get(filename))
}
toCOL=min(unlist(ncol_length))

ExperimentalData = list()
for(i in 1:length(AllPlates)){
  filename = AllPlates[[i]]
  ExperimentalData[[i]] = get(filename)
}

Experiment_edited = list()
for(i in 1:length(AllPlates)){
  filename = AllPlates[[i]]
  Experiment_edited[[i]] = get(filename)[,1:toCOL]
}

Experiment_df = data.frame()
for(i in 1:length(Experiment_edited)){
  Experiment_df = rbind(Experiment_df, Experiment_edited[[i]])
}

write.csv(Experiment_df, file = "Experiment_df.csv")

GroFitResults = data.frame(Row = character(),
                           Column = character(),
                           Set = numeric(),
                           Replicate = character(), 
                           SpNum = numeric(),
                           PUNum = numeric(),
                           Ship = numeric(),
                           Species = character(),
                           IDNum = numeric(),
                           Plate = character(), 
                           Treatment = character(),
                           Media = character(), 
                           Lag = numeric(), 
                           Growth = numeric(), 
                           Saturation = numeric())
for(i in 1:nrow(Experiment_df)){
  GroFitResults[i,1] = Experiment_df[i,1]
  GroFitResults[i,2] = Experiment_df[i,2] # MODIFY TO RELEVANT PLATE INFORMATION
  GroFitResults[i,3] = Experiment_df[i,3] # MODIFY TO RELEVANT PLATE INFORMATION
  GroFitResults[i,4] = Experiment_df[i,4]
  GroFitResults[i,5] = Experiment_df[i,5] # MODIFY TO RELEVANT PLATE INFORMATION
  GroFitResults[i,6] = Experiment_df[i,6]
  GroFitResults[i,7] = Experiment_df[i,7]
  GroFitResults[i,8] = Experiment_df[i,8]
  GroFitResults[i,9] = Experiment_df[i,9]
  GroFitResults[i,10] = Experiment_df[i,10]
  GroFitResults[i,11] = Experiment_df[i,11]
  GroFitResults[i,12] = Experiment_df[i,12]
  TimeData <- t(as.matrix(cbind(seq(0, 2*(ncol(Experiment_df)-17), 2))))
  GrowthData <- Experiment_df[i,17:ncol(Experiment_df)]
  info = data.frame(matrix(0,ncol = 3))
  GrowthData <- cbind(info, GrowthData)
  GrowthData <- t(data.frame(as.numeric(GrowthData)))
  TEMP = grofit(TimeData, GrowthData, control=grofit.control(suppress.messages = TRUE, fit.opt = "m", interactive = FALSE, model.type=c("logistic"), nboot.gc= 0, smooth.gc  = 5))
  GroFitResults[i,13] = TEMP[["gcFit"]][["gcTable"]][["lambda.model"]]
  GroFitResults[i,14] = TEMP[["gcFit"]][["gcTable"]][["mu.model"]]
  GroFitResults[i,15] = TEMP[["gcFit"]][["gcTable"]][["A.model"]]
}

write.csv(GroFitResults, file =  paste("GroFitResults", "_", Sys.Date(), ".csv", sep = ""), row.names = FALSE) # # unique(GroFit_df$Replicate), "_", Sys.Date(), ".csv", sep = ""), 


TimeData <- t(as.matrix(cbind(seq(0, 2*(ncol(Experiment_df)-17), 2))))
SPECIES <- unique(Experiment_df$`PU#`)
CARBON <- unique(Experiment_df$Media)

pdf(file = paste("GrowthCurves2", "_", Sys.Date(), ".pdf", sep = ""), height = 30, width = 30) 
for(i in 1:length(SPECIES)){
  rows <-  which(Experiment_df$`PU#`==SPECIES[[i]])
  tempData <- Experiment_df[rows,]
  colnames(tempData)[17:ncol(tempData)] = TimeData
  tempData = tempData[,-c(1,2,3,5,6,7,9,10,11,13,14,15,16)]
  ROWS =  which(GroFitResults$PUNum==SPECIES[[i]])
  groTemp <- GroFitResults[rows,]
  TEMP = merge(groTemp, tempData, by.y = c("Media", "Rep.x"), by.x = c("Media", "Replicate"))
  TEMP = TEMP[,-c(3, 4, 5, 6, 8, 10, 11, 12, 16)]
  TEMPMELT = melt(TEMP, id.vars = c("Media", "Replicate", "PUNum", "Species.x", "Lag", "Growth", "Saturation"))
  colnames(TEMPMELT) = c("Media","Rep", "PU", "Species", "Lag", "Growth", "Saturation", "Time", "OD") 
  TEMPMELT$OD = as.numeric(TEMPMELT$OD)
  
  NOCARBON = TEMPMELT[which(TEMPMELT$Media == "No Carbon"),]
  NOCARBON = NOCARBON[,-which(colnames(NOCARBON) == "Media")]
  
  curve = ggplot(TEMPMELT, aes(x = Time, y = OD, colour=factor(Rep)))
  curve = curve + geom_point(pch = 19)
  curve = curve + facet_wrap(Media~., scales = "free_y")
  curve = curve + geom_hline(data = TEMPMELT, aes(yintercept= Saturation, colour = factor(Rep))) 
  curve = curve + geom_abline(data = TEMPMELT, aes(slope = Growth, intercept=-Growth*Lag, colour = factor(Rep)))
  curve = curve + geom_point(data = NOCARBON, aes(x = Time, y = OD), pch = 19, color = "black")
  curve = curve + xlab("Time")+ylab("Absorbance")
  curve = curve + ggtitle(paste(unique(TEMPMELT$Species), unique(TEMPMELT$PU), sep = " "))
  curve = curve + theme(axis.title.x = element_text(size = 10, face = "bold"),
                               axis.title.y = element_text(size = 8, face = "bold"),
                               panel.background = element_rect(fill = "white", color = "black"),
                               panel.grid = element_blank(), 
                               axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5), 
                               axis.text.y = element_text(size = 8),
                               strip.background =element_blank(), 
                               strip.text.x = element_text(size = 12)
                        )
  print(curve)
}
dev.off()

NOGrowthPU = c(34742, 34650, 40952, 35006, 34683, 34761, 34663, 34589, 34607, 41811, 37898, 37889, 34727, 37897, 41819, 35642, 37835, 40986, 34755, 41845, 35032, 26124, 26126, 34959, 35261, 38082, 35297, 37871, 38385, 35665, 35022, 40969, 40964, 38369, 40993, 35710, 35302, 40990, 37848, 35279, 34605, 36579, 40960, 37862, 34731, 40957, 34655, 41847, 34591, 34860, 26279, 37836, 34751, 35638, 26197, 41786, 37906, 26049, 38052, 40961, 34741, 37920, 38399, 26117, 34707, 38317, 34872, 34760, 34935, 40962, 34717, 34682, 34688, 37869, 41785, 38394, 41813, 34867, 34894, 38050, 34704, 34698, 35711, 34748, 34891, 37903, 34687, 40971, 35038, 34759, 35323, 41783, 34597, 38342, 35280, 38401, 40945, 34645, 34602, 34686, 41797, 35281, 37922, 37879, 41799, 37901, 34667, 41719, 40944, 37921, 40915, 34757, 35658, 35720, 38315, 40943, 37886, 38047, 35678, 41796, 35322, 35693, 40998, 35266, 37849, 38051, 40951, 34647, 34746, 38379, 40936, 34666, 34750, 34752, 41784, 34753, 34673, 41787, 40920, 34571, 34984, 35626, 34975, 34749, 34582, 37929, 37850, 34706, 34649, 40966, 40938, 40955, 34744, 35278, 34874, 34635, 34652, 26240, 35629, 40996, 37923, 41802, 38361, 41733, 35295, 37876, 40974, 34580, 38376, 40991, 35270, 37918, 41770, 35690, 41855, 41794, 40992, 35654, 41862, 35294, 34596, 41822, 34976, 37891, 34747, 35308, 34670, 34575)

pdf(file = paste("NoGrowPU", "_", Sys.Date(), ".pdf", sep = ""), height = 30, width = 30) 
for(i in 1:length(NOGrowthPU)){
  rows <-  which(Experiment_df$`PU#`==NOGrowthPU[[i]])
  tempData <- Experiment_df[rows,]
  colnames(tempData)[17:ncol(tempData)] = TimeData
  tempData = tempData[,-c(1,2,3,5,6,7,9,10,11,13,14,15,16)]
  ROWS =  which(GroFitResults$PUNum==NOGrowthPU[[i]])
  groTemp <- GroFitResults[rows,]
  TEMP = merge(groTemp, tempData, by.y = c("Media", "Rep.x"), by.x = c("Media", "Replicate"))
  TEMP = TEMP[,-c(3, 4, 5, 6, 8, 10, 11, 12, 16)]
  TEMPMELT = melt(TEMP, id.vars = c("Media", "Replicate", "PUNum", "Species.x", "Lag", "Growth", "Saturation"))
  colnames(TEMPMELT) = c("Media","Rep", "PU", "Species", "Lag", "Growth", "Saturation", "Time", "OD") 
  TEMPMELT$OD = as.numeric(TEMPMELT$OD)
  curve = ggplot(TEMPMELT, aes(x = Time, y = OD, colour=factor(Rep)))
  curve = curve + geom_point(pch = 19)
  curve = curve + facet_wrap(Media~., scales = "free_y")
  curve = curve + geom_hline(data = TEMPMELT, aes(yintercept= Saturation, colour = factor(Rep))) 
#  curve = curve + geom_abline(data = TEMPMELT, aes(slope = Growth, intercept=-Growth*Lag, colour = factor(Rep)))
  curve = curve + xlab("Time")+ylab("Absorbance")
  curve = curve + ggtitle(paste(unique(TEMPMELT$Species), unique(TEMPMELT$PU), sep = " "))
  curve = curve + theme(axis.title.x = element_text(size = 10, face = "bold"),
                        axis.title.y = element_text(size = 8, face = "bold"),
                        panel.background = element_rect(fill = "white", color = "black"),
                        panel.grid = element_blank(), 
                        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5), 
                        axis.text.y = element_text(size = 8),
                        strip.background =element_blank(), 
                        strip.text.x = element_text(size = 12)
  )
  print(curve)
}
dev.off()

NoCarbGrowth = c(41834, 34610, 34700, 35335, 34614, 38105, 37873, 
                 35689, 41701, 41850, 38076, 35681, 35703, 35336, 
                 35659, 38057, 38056, 34576, 35251, 38375, 34991, 
                 26145, 34933, 38083, 35323, 41801, 35268, 38063, 
                 34932, 34968, 38038, 26190, 34695, 37857, 34615, 
                 26167, 38409, 40976, 37926, 35318, 38334, 34924, 
                 35337, 38373, 35262, 40963, 38058, 35293, 38344, 
                 38059, 35286, 26143, 35244, 34710, 37872, 26240, 
                 38041, 34996, 35708, 38350, 26073, 35679, 37855, 
                 41817, 35715, 34938, 37852, 40940, 26137, 38364, 
                 38078, 34897)

pdf(file = paste("NoCarbGrowth", "_", Sys.Date(), ".pdf", sep = ""), height = 30, width = 30) 
for(i in 1:length(NoCarbGrowth)){
  rows <-  which(Experiment_df$`PU#`==NoCarbGrowth[[i]])
  tempData <- Experiment_df[rows,]
  colnames(tempData)[17:ncol(tempData)] = TimeData
  tempData = tempData[,-c(1,2,3,5,6,7,9,10,11,13,14,15,16)]
  ROWS =  which(GroFitResults$PUNum==NoCarbGrowth[[i]])
  groTemp <- GroFitResults[rows,]
  TEMP = merge(groTemp, tempData, by.y = c("Media", "Rep.x"), by.x = c("Media", "Replicate"))
  TEMP = TEMP[,-c(3, 4, 5, 6, 8, 10, 11, 12, 16)]
  TEMPMELT = melt(TEMP, id.vars = c("Media", "Replicate", "PUNum", "Species.x", "Lag", "Growth", "Saturation"))
  colnames(TEMPMELT) = c("Media","Rep", "PU", "Species", "Lag", "Growth", "Saturation", "Time", "OD") 
  TEMPMELT$OD = as.numeric(TEMPMELT$OD)
  NOCARBON = TEMPMELT[which(TEMPMELT$Media == "No Carbon"),]
  NOCARBON = NOCARBON[,-which(colnames(NOCARBON) == "Media")]
  curve = ggplot(TEMPMELT, aes(x = Time, y = OD, colour=factor(Rep)))
  curve = curve + geom_point(pch = 19)
  curve = curve + facet_wrap(Media~., scales = "free_y")
  curve = curve + geom_hline(data = TEMPMELT, aes(yintercept= Saturation, colour = factor(Rep))) 
  #  curve = curve + geom_abline(data = TEMPMELT, aes(slope = Growth, intercept=-Growth*Lag, colour = factor(Rep)))
  curve = curve + geom_point(data = NOCARBON, aes(x = Time, y = OD), pch = 19, color = "black")
  curve = curve + xlab("Time")+ylab("Absorbance")
  curve = curve + ggtitle(paste(unique(TEMPMELT$Species), unique(TEMPMELT$PU), sep = " "))
  curve = curve + theme(axis.title.x = element_text(size = 10, face = "bold"),
                        axis.title.y = element_text(size = 8, face = "bold"),
                        panel.background = element_rect(fill = "white", color = "black"),
                        panel.grid = element_blank(), 
                        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5), 
                        axis.text.y = element_text(size = 8),
                        strip.background =element_blank(), 
                        strip.text.x = element_text(size = 12)
  )
  print(curve)
}
dev.off()

Rechecks = c(34718, 26150, 34675, 34957, 34962, 34589, 35048, 35642, 34966, 38365, 37835, 40986, 41845, 37845, 26124, 34904, 38389, 38082, 37871, 35022, 38101, 35008, 35012, 40993, 34934, 41846, 35258, 34626, 38314, 36579, 41841, 38316, 40960, 34956, 40983, 34603, 34860, 26279, 35686, 34864, 37916, 41850, 37890, 38381, 41859, 26197, 38028, 26049, 35269, 37854, 35703, 35648, 38399, 26050, 26117, 34965, 38026, 26068, 34760, 34935, 40962, 34717, 34682, 38354, 34688, 37869, 38056, 41785, 38113, 38394, 34576, 34867, 34961, 34894, 34631, 40933, 38050, 40985, 34880, 34704, 34698, 34611, 35711, 35331, 37903, 34991, 34687, 34954, 40971, 35038, 35694, 34759, 35323, 41783, 34597, 41835, 38342, 35280, 38401, 34884, 34604, 37846, 41778, 34645, 38371, 41842, 34985, 38348, 34602, 38345, 34686, 41797, 37922, 38038, 41799, 34689, 37901, 26190, 34667, 37921, 34757, 38372, 38047, 35678, 35322, 34701, 40998, 34999, 37842, 35243, 40951, 38379, 40936, 34666, 34750, 34752, 41784, 34868, 34753, 34673, 41787, 40920, 35626, 34975, 34749, 41691, 38059, 37850, 34706, 40966, 40938, 40906, 34651, 37899, 34874, 34635, 34652, 35629, 37874, 40996, 37923, 41802, 41733, 35295, 37876, 40974, 40991, 35270, 37918, 41770, 35690, 41794, 40992, 35314, 35654, 34705, 41862, 34596, 41822, 37891)
pdf(file = paste("Rechecks", "_", Sys.Date(), ".pdf", sep = ""), height = 30, width = 30) 
for(i in 1:length(Rechecks)){
  rows <-  which(Experiment_df$`PU#`==Rechecks[[i]])
  tempData <- Experiment_df[rows,]
  colnames(tempData)[17:ncol(tempData)] = TimeData
  tempData = tempData[,-c(1,2,3,5,6,7,9,10,11,13,14,15,16)]
  ROWS =  which(GroFitResults$PUNum==Rechecks[[i]])
  groTemp <- GroFitResults[rows,]
  TEMP = merge(groTemp, tempData, by.y = c("Media", "Rep.x"), by.x = c("Media", "Replicate"))
  TEMP = TEMP[,-c(3, 4, 5, 6, 8, 10, 11, 12, 16)]
  TEMPMELT = melt(TEMP, id.vars = c("Media", "Replicate", "PUNum", "Species.x", "Lag", "Growth", "Saturation"))
  colnames(TEMPMELT) = c("Media","Rep", "PU", "Species", "Lag", "Growth", "Saturation", "Time", "OD") 
  TEMPMELT$OD = as.numeric(TEMPMELT$OD)
  NOCARBON = TEMPMELT[which(TEMPMELT$Media == "No Carbon"),]
  NOCARBON = NOCARBON[,-which(colnames(NOCARBON) == "Media")]
  curve = ggplot(TEMPMELT, aes(x = Time, y = OD, colour=factor(Rep)))
  curve = curve + geom_point(pch = 19)
  curve = curve + facet_wrap(Media~., scales = "free_y")
  curve = curve + geom_hline(data = TEMPMELT, aes(yintercept= Saturation, colour = factor(Rep))) 
  #  curve = curve + geom_abline(data = TEMPMELT, aes(slope = Growth, intercept=-Growth*Lag, colour = factor(Rep)))
  curve = curve + geom_point(data = NOCARBON, aes(x = Time, y = OD), pch = 19, color = "black")
  curve = curve + xlab("Time")+ylab("Absorbance")
  curve = curve + ggtitle(paste(unique(TEMPMELT$Species), unique(TEMPMELT$PU), sep = " "))
  curve = curve + theme(axis.title.x = element_text(size = 10, face = "bold"),
                        axis.title.y = element_text(size = 8, face = "bold"),
                        panel.background = element_rect(fill = "white", color = "black"),
                        panel.grid = element_blank(), 
                        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5), 
                        axis.text.y = element_text(size = 8),
                        strip.background =element_blank(), 
                        strip.text.x = element_text(size = 12)
  )
  print(curve)
}
dev.off()

SpeciesChecks = read.csv("SpeciesChecks.csv", header = TRUE, check.names = F)

pdf(file = paste("SpecificChecks", "_", Sys.Date(), ".pdf", sep = "")) 
for(i in 1:nrow(SpeciesChecks)){
  rows <-  which(Experiment_df$`PU#`== SpeciesChecks[i,"PU#"] & Experiment_df$Media == SpeciesChecks[i,"Media"])
  tempData <- Experiment_df[rows,]
  colnames(tempData)[17:ncol(tempData)] = TimeData
  tempData = tempData[,-c(1,2,3,5,6,7,9,10,11,13,14,15,16)]

  ROWS =  which(GroFitResults$PUNum == SpeciesChecks[i,"PU#"] & GroFitResults$Media == SpeciesChecks[i,"Media"])
  groTemp <- GroFitResults[ROWS,]
  TEMP = merge(groTemp, tempData, by.y = c("Media", "Rep.x"), by.x = c("Media", "Replicate"))
  TEMP = TEMP[,-c(3, 4, 5, 6, 8, 10, 11, 12, 16)]
  if(nrow(TEMP) > 0){
    TEMPMELT = melt(TEMP, id.vars = c("Media", "Replicate", "PUNum", "Species.x", "Lag", "Growth", "Saturation"))
    colnames(TEMPMELT) = c("Media","Rep", "PU", "Species", "Lag", "Growth", "Saturation", "Time", "OD") 
    TEMPMELT$OD = as.numeric(TEMPMELT$OD)
    
    NOCARBON = Experiment_df[which(Experiment_df$Media == "No Carbon" & Experiment_df$`PU#` == SpeciesChecks[i,"PU#"]),]
    colnames(NOCARBON)[17:ncol(NOCARBON)] = TimeData
    
    NOCARBON = NOCARBON[,-c(1,2,3,5,7,9,10,11,13,14,15,16)]
    NOCARBMELT = melt(NOCARBON, id.vars = c("Media", "Rep.x", "PU#", "Species"))
    colnames(NOCARBMELT)[5:6] = c("Time", "OD")
    NOCARBMELT$OD = as.numeric(NOCARBMELT$OD)
    
    curve = ggplot(TEMPMELT, aes(x = Time, y = OD, colour=factor(Rep)))
    curve = curve + geom_point(pch = 19)
    curve = curve + geom_hline(data = TEMPMELT, aes(yintercept= Saturation, colour = factor(Rep))) 
    curve = curve + geom_point(data = NOCARBMELT, aes(x = Time, y = OD), group = "Rep.x" ,pch = 19, color = "black")
    curve = curve + xlab("Time")+ylab("Absorbance")
    curve = curve + ggtitle(paste(unique(TEMPMELT$Species), unique(TEMPMELT$PU), unique(TEMPMELT$Media), sep = " "))
    curve = curve + theme(axis.title.x = element_text(size = 10, face = "bold"),
                          axis.title.y = element_text(size = 8, face = "bold"),
                          panel.background = element_rect(fill = "white", color = "black"),
                          panel.grid = element_blank(), 
                          axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5), 
                          axis.text.y = element_text(size = 8),
                          strip.background =element_blank(), 
                          strip.text.x = element_text(size = 12)
    )
    print(curve)
  }  
}
  
dev.off()
