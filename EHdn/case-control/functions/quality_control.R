# Function that creates the age distribution across case and control datasets, given a R environment that has defined under `case-control`
plotting_age_distribution <- function(working_directory, environment_file){
  # libraries
  library(dplyr)
  library(ggplot2)
  
  # defining working directory
  setwd(working_directory)
  
  # Load ANALYSIS case-control environment Rdata
  load(environment_file)
  
  # Pilot and Main cases - taking platekey, sex and age
  selected_pilot_cases = pilot_cases %>% 
    filter(plateKey %in% l_pilot_cases) %>%
    select(plateKey, sex, yearOfBirth)
  dim(selected_pilot_cases)
  
  selected_main_cases = main_cases %>%
    filter(platekey %in% l_main_cases) %>%
    select(platekey, participant_phenotypic_sex, year_of_birth)
  dim(selected_main_cases)
  
  selected_main_cases = unique(selected_main_cases)
  dim(selected_main_cases)
  
  colnames(selected_pilot_cases) = c("platekey", "sex", "YOB")
  colnames(selected_main_cases) = c("platekey", "sex", "YOB")
  
  selected_cases = rbind(selected_pilot_cases,
                         selected_main_cases)
  
  dim(selected_cases)
  
  # Pilot and Main controls
  selected_pilot_controls = pilot_controls %>%
    filter(plateKey %in% l_pilot_controls) %>%
    select(plateKey, sex, yearOfBirth)
  dim(selected_pilot_controls)
  
  selected_pilot_controls = unique(selected_pilot_controls)
  dim(selected_pilot_controls)
  
  selected_main_controls = main_controls %>%
    filter(platekey %in% l_main_controls) %>%
    select(platekey, participant_phenotypic_sex, year_of_birth)
  dim(selected_main_controls)
  
  selected_main_controls = unique(selected_main_controls)
  dim(selected_main_controls)
  
  colnames(selected_pilot_controls) = c("platekey", "sex", "YOB")
  colnames(selected_main_controls) = c("platekey", "sex", "YOB")
  
  selected_controls = rbind(selected_pilot_controls,
                            selected_main_controls)
  dim(selected_controls)
  
  # Define each `selected_cases` and `selected_controls` with the group name
  selected_cases$group = rep("case", length(selected_cases$platekey))
  selected_controls$group = rep("control", length(selected_controls$platekey))
  
  merged_pat = rbind(selected_cases,
                     selected_controls)
  
  dim(merged_pat)
  
  merged_pat = merged_pat %>%
    group_by(platekey) %>%
    mutate(age = 2020 - YOB) %>%
    ungroup() %>%
    as.data.frame()
  
  # Plot age distribution
  age_distribution = ggplot(merged_pat,
         aes(x = age, fill = group)) + geom_density(alpha=0.25)
  
  plot_output_name = paste(working_directory, "age_distribution.png", sep = "output/")
  png(plot_output_name)
  print(age_distribution)
  dev.off()


}
