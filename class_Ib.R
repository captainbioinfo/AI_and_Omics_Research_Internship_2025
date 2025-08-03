# --------------------------------------------------------------------------
#### Tasks ####

# 1. Set Working Directory
# Create a new folder on your computer "AI_Omics_Internship_2025".

# 2. Create Project Folder

# In RStudio, create a new project named "Module_I" in your "AI_Omics_Internship_2025" folder.

# Inside the project directory, create the following subfolders using R code:
# raw_data, clean_data, scripts, results or Tasks, plots etc

# ---------------------------------------------------------------------------
# 3. Download "patient_info.csv" dataset from GitHub repository

# load the dataset into your R environment

# Inspect the structure of the dataset using appropriate R functions

# Identify variables with incorrect or inconsistent data types.

# Convert variables to appropriate data types where needed

# Create a new variable for smoking status as a binary factor:

# 1 for "Yes", 0 for "No"

# Save the cleaned dataset in your clean_data folder with the name patient_info_clean.csv
# Save your R script in your script folder with name "class_Ib"
# Upload "class_Ib" R script into your GitHub repository
# Save your R workspace with name "YourName_Class_Ib_Assignment.RData" and upload the saved .RData file in the assignment submission form.

########################## STARTING ################################
getwd()

# Create directory
dir.create("AI_Omics_Internship_2025")

# Seting &  checking working directory 
setwd("AI_Omics_Internship_2025")
getwd()

# Create project folder
dir.create("Module_I")
setwd("Module_I")

# Creating directories inside the Module I folder
sapply(c("raw_data", "clean_data", "scripts", "results","plots"), dir.create)

setwd("raw_data")

## Downloading raw data file from Github
url <- "https://raw.githubusercontent.com/AI-Biotechnology-Bioinformatics/AI_and_Omics_Research_Internship_2025/main/patient_info.csv"
data <- read.csv(url)
#
download.file(url, destfile = "patient_info.csv", mode = "wb")

# Inspect the structure of the dataset using appropriate R functions
head(data)
tail(data)
colnames(data)
rownames(data)
data$patient_id
data$age
data$smoker
#### Check structure of each column
str(data)  
selected_data <- data[2, 1:4]

####
# Convert gender and diagnosis to factors
data$gender <- as.factor(data$gender)
data$diagnosis <- as.factor(data$diagnosis)
data$smoker <- as.character(data$smoker)

str(data)
data$smoker_binary <- data$smoker
data$smoker_binary <- ifelse(tolower(data$smoker_binary) == "yes", 0,
                             ifelse(tolower(data$smoker_binary) == "no", 1, NA))

## Final result checking 
table(data$smoker_binary)
table(data$smoker)
data_clean <- data[-6]

# Saving the final result
setwd("..")
setwd("clean_data")
write.csv(data_clean, "patient_info_clean.csv", row.names = FALSE)

