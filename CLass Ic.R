#### Try It Yourself ####
# Practice Exercises 

# ----------------------------------------------------------------------------------------------------------------

# 1. Check Cholesterol level (using if) 
# Write an If statement to check cholesterol level is greater than 240, 
# if true, it will prints “High Cholesterol”

cholesterol <- 230

# ----------------------------------------------------------------------------------------------------------------

# 2. Blood Pressure Status (using if...else)
# Write an if…else statement to check if blood pressure is normal.
# If it’s less than 120, print: “Blood Pressure is normal”
# If false then print: “Blood Pressure is high”

Systolic_bp <- 130

# ----------------------------------------------------------------------------------------------------------------

# 3. Automating Data Type Conversion with for loop

# Use patient_info.csv data and metadata.csv
# Perform the following steps separately on each dataset (patient_info.csv data and metadata.csv)
# Create a copy of the dataset to work on.
# Identify all columns that should be converted to factor type.
# Store their names in a variable (factor_cols).

# Example: factor_cols <- c("gender", "smoking_status")

# Use a for loop to convert all the columns in factor_cols to factor type.
# Pass factor_cols to the loop as a vector.

# Hint:
# for (col in factor_cols) {
#   data[[col]] <- as.factor(data[[col]])  # Replace 'data' with the name of your dataset
# }

# ----------------------------------------------------------------------------------------------------------------

# 4. Converting Factors to Numeric Codes

# Choose one or more factor columns (e.g., smoking_status).
# Convert "Yes" to 1 and "No" to 0 using a for loop.

# Hint:
# binary_cols <- c("smoking_status")   # store column names in a vector
# use ifelse() condition inside the loop to replace Yes with 1 and No with 0
# for (col in binary_cols) {
#   data[[col]] <- # insert your ifelse() code here
# }

# ----------------------------------------------------------------------------------------------------------------

#  Verification:
#    Compare the original and modified datasets to confirm changes.
#    str(original_data)
#    str(data)

# ----------------------------------------------------------------------------------------------------------------

# Final Note:
# All instructions are written as comments.
# For actual code execution, remove the # symbol from each line you want to run.
############################# Practice Session #############################################################
# AI & Omics Internship - Module I
# Practice Exercises: If Statements, Data Conversion, and Loops
# ---------------------------------------------------------------------------

# 1. Check Cholesterol level (using if)
# Task:
#   - Write an If statement to check if cholesterol is greater than 240.
#   - If true, print “High Cholesterol”.
# ---------------------------------------------------------------------------
## 1.
cholesterol <- 230

if (cholesterol > 240) {
  print("High Cholesterol")
} 

## To check
cholesterol <- 260

if (cholesterol > 240) {
  print("High Cholesterol")
} 

## 2. 
# 2. Blood Pressure Status (using if...else)
# Task:
#   - Write an if…else statement to check if blood pressure is normal.
#   - If systolic blood pressure is less than 120, print: “Blood Pressure is normal”.
#   - Otherwise, print: “Blood Pressure is high”.
# ---------------------------------------------------------------------------
Systolic_bp <- 130

if (Systolic_bp > 120) {
  print("Blood Pressure is high")
} else {
  print("Blood Pressure is normal")
}

### To check 

Systolic_bp <- 100

if (Systolic_bp > 120) {
  print("Blood Pressure is high")
} else {
  print("Blood Pressure is normal")
}

## 3.
# 3. Automating Data Type Conversion with for loop
# Task:
#   - Work with two datasets: patient_info.csv and Metadata.csv.
#   - Create a copy of each dataset.
#   - Identify categorical columns that should be converted to factors.
#   - Store their names in a variable (factor_cols).
#   - Use a for loop to convert these columns to factor type.
# ---------------------------------------------------------------------------
# Set working directory

setwd("~/AI_Omics_Internship_2025/Module_I")
# Load patient_info dataset
data <- read.csv("raw_data/patient_info.csv")

head(data)

## Downloading Metadata data file from Github
url <- "https://raw.githubusercontent.com/AI-Biotechnology-Bioinformatics/AI_and_Omics_Research_Internship_2025/main/Metadata.csv"
download.file(url, destfile = "raw_data/Metadata.csv", mode = "wb")

Metadata <- read.csv("raw_data/Metadata.csv")

head(Metadata)
##
# -----------------------------
# Patient Info Dataset
# -----------------------------

# Copy of dataset
patient_copy <- data

# Identify categorical columns
factor_cols_patient <- c("patient_id", "gender", "diagnosis", "smoker")

# Convert to factors using a for loop
for (col in factor_cols_patient) {
  patient_copy[[col]] <- as.factor(patient_copy[[col]])
}

# Check structure
str(patient_copy)


# -----------------------------
# Metadata Dataset
# -----------------------------

# Copy of dataset
metadata_copy <- Metadata

# Identify categorical columns
factor_cols_metadata <- c("name", "height", "gender")

# Convert to factors using a for loop
for (col in factor_cols_metadata) {
  metadata_copy[[col]] <- as.factor(metadata_copy[[col]])
}

# Check structure
str(metadata_copy)

## 4. 
# 4. Converting Factors to Numeric Codes
# Task:
#   - Choose one or more factor columns (e.g., smoker).
#   - Convert "Yes" to 1 and "No" to 0 using a for loop.
#   - Use ifelse() inside the loop for replacement.
# ---------------------------------------------------------------------------

# Make a copy to preserve original
patient_copy2 <- data

# Columns to convert
binary_cols <- c("smoker")

# Loop through and convert Yes/No → 1/0
for (col in binary_cols) {
  patient_copy2[[col]] <- ifelse(patient_copy2[[col]] == "Yes", 1, 0)
}

# Verification

# Compare before and after
str(data)          # original dataset
str(patient_copy2) # modified dataset
head(patient_copy2)

#
# Save modified patient_info dataset
#write.csv(patient_copy2, "clean_data/patient_info_modified.csv", row.names = FALSE)

# Save modified metadata dataset
#write.csv(metadata_copy, "clean_data/metadata_modified.csv", row.names = FALSE)

########################## Session ended #####################################
