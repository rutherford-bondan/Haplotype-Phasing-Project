##### Authors: Audrey Magsig and Austin Bond
##### Date: 03/19/2020
##### CS M124 Class Project

##### PACKAGES #####

packages <- c("readr", "tibble", "plyr", "dplyr", "tidyr", "stringr")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}

library(readr)
library(tibble)
library(plyr)
library(dplyr)
library(tidyr)
#library(ggfortify)
#library(mice)
library(stringr)

##### FUNCTIONS #####

try_to_phase <- function(sub, heterozyg, homozyg){
  
  phased <- sub
  
  haplo <- c()
  for (i in homozyg$index) {
    phased[sub[,i] == 0, i] = "0|0"
    phased[sub[,i] == 2, i] = "1|1"
    
    name <- subjects[i]
    first <- paste0(name, "L")
    second <- paste0(name, "R")
    
    partial <- phased %>% separate(!!name, into = c(first, second)) %>% select(first, second)
    
    #new <- str_c(partial[,1], collapse = "")
    # if (!(new %in% haplo)) {
    #   haplo <- c(haplo, str_c(partial[,1], collapse = ""))
    # }
    new <- list(as.numeric(partial[,1]))
    
    haplo <- c(haplo, new)
  }
  
  #long <- haplo
  
  haplo <- unique(haplo)
  
  all_potential <- c()
  
  temp <- phased
  
  phased <- temp
  unphased <- c()
  
  set.seed(1234)
  
  for(i in heterozyg$index){
    is_phased = FALSE
    
    potential <- c()
    
    for(j in 1:length(haplo)){
      
      other <- list(sub[,i] - unlist(haplo[j]))
      if (other %in% haplo) {
        phased[,i] = paste(unlist(haplo[j]), unlist(other), sep = "|")
        is_phased = TRUE
        break
      }else {
        
        if(all(unlist(other) %in% c(0, 1))){
          potential <- c(potential, other) 
        }
      }
      
    }
    
    # if (length(potential) > 1) {
    #   all_potential <- c(all_potential, potential)
    #   #unphased <- c(unphased, i)
    #   #next
    #   
    # }
    
    if(is_phased == FALSE){
      
      pick <- sample(1:length(potential), 1)
      
      orig <- list(sub[,i] - unlist(potential[pick]))
      if(orig %in% haplo){
        phased[,i] = paste(unlist(orig), unlist(potential[pick]), sep = "|")
        haplo <- c(haplo, potential[pick])
      } else{
        unphased <- c(unphased, i)
      }
      
      #unphased <- c(unphased, i)
    }
    
  }
  
  #first_block <- phased
  
  # possible <- c()
  # for (i in 1:length(all_potential)) {
  #   new <- str_c(unlist(all_potential[i]), collapse = "")
  #   possible <- c(possible, new)
  # }
  #possible <- sort(possible)
  
  redo <- unphased
  
  unphased <- c()
  
  if (length(redo) > 0) {
    
    for(i in redo){
      is_phased = FALSE
      potential <- c()
      
      for(j in 1:length(haplo)){
        
        other <- list(sub[,i] - unlist(haplo[j]))
        if (other %in% haplo) {
          phased[,i] = paste(unlist(haplo[j]), unlist(other), sep = "|")
          is_phased = TRUE
          break
        }else {
          
          if(all(unlist(other) %in% c(0, 1))){
            potential <- c(potential, other)
          }
        }
        
      }
      
      # if (length(potential) > 1) {
      #   unphased <- c(unphased, i)
      #   next
      # }
      
      if(is_phased == FALSE){
        
        # if (length(potential) == 0 & length(sub[sub[,i] == 1, i]) == 1) {
        #   ### FIX THIS LATER ###
        #   phased[sub[,i] == 1, i] = "1|0"
        #   phased[sub[,i] == 0, i] = "0|0"
        #   phased[sub[,i] == 2, i] = "1|1"
        #   break
        # }
        
        pick <- sample(1:length(potential), 1)
        
        orig <- list(sub[,i] - unlist(potential[pick]))
        if(orig %in% haplo){
          phased[,i] = paste(unlist(orig), unlist(potential[pick]), sep = "|")
          haplo <- c(haplo, potential[pick])
        } else{
          unphased <- c(unphased, i)
        }
        
        #unphased <- c(unphased, i)
      }
    }
    
  }
  
  if (length(unphased) > 0) {
    return(NULL)
  } else{
    return(phased)
  }
  
}

fix_for_unique_snp <- function(sub, heterozyg, homozyg){
  
  phased <- sub
  
  haplo <- c()
  for (i in homozyg$index) {
    phased[sub[,i] == 0, i] = "0|0"
    phased[sub[,i] == 2, i] = "1|1"
    
    name <- subjects[i]
    first <- paste0(name, "L")
    second <- paste0(name, "R")
    
    partial <- phased %>% separate(!!name, into = c(first, second)) %>% select(first, second)
    
    #new <- str_c(partial[,1], collapse = "")
    # if (!(new %in% haplo)) {
    #   haplo <- c(haplo, str_c(partial[,1], collapse = ""))
    # }
    new <- list(as.numeric(partial[,1]))
    
    haplo <- c(haplo, new)
  }
  
  #long <- haplo
  
  haplo <- unique(haplo)
  
  all_potential <- c()
  
  temp <- phased
  
  phased <- temp
  unphased <- c()
  
  set.seed(1234)
  
  for(i in heterozyg$index){
    is_phased = FALSE
    
    potential <- c()
    
    for(j in 1:length(haplo)){
      
      other <- list(sub[,i] - unlist(haplo[j]))
      if (other %in% haplo) {
        phased[,i] = paste(unlist(haplo[j]), unlist(other), sep = "|")
        is_phased = TRUE
        break
      }else {
        
        if(all(unlist(other) %in% c(0, 1))){
          potential <- c(potential, other) 
        }
      }
      
    }
    
    # if (length(potential) > 1) {
    #   all_potential <- c(all_potential, potential)
    #   #unphased <- c(unphased, i)
    #   #next
    #   
    # }
    
    if(is_phased == FALSE){
      
      pick <- sample(1:length(potential), 1)
      
      orig <- list(sub[,i] - unlist(potential[pick]))
      if(orig %in% haplo){
        phased[,i] = paste(unlist(orig), unlist(potential[pick]), sep = "|")
        haplo <- c(haplo, potential[pick])
      } else{
        unphased <- c(unphased, i)
      }
      
      #unphased <- c(unphased, i)
    }
    
  }
  
  #first_block <- phased
  
  # possible <- c()
  # for (i in 1:length(all_potential)) {
  #   new <- str_c(unlist(all_potential[i]), collapse = "")
  #   possible <- c(possible, new)
  # }
  #possible <- sort(possible)
  
  redo <- unphased
  
  unphased <- c()
  
  if (length(redo) > 0) {
    
    for(i in redo){
      is_phased = FALSE
      potential <- c()
      
      for(j in 1:length(haplo)){
        
        other <- list(sub[,i] - unlist(haplo[j]))
        if (other %in% haplo) {
          phased[,i] = paste(unlist(haplo[j]), unlist(other), sep = "|")
          is_phased = TRUE
          break
        }else {
          
          if(all(unlist(other) %in% c(0, 1))){
            potential <- c(potential, other)
          }
        }
        
      }
      
      # if (length(potential) > 1) {
      #   unphased <- c(unphased, i)
      #   next
      # }
      
      if(is_phased == FALSE){
        
        if (length(potential) == 0) {
          for (k in which(sub[,i] == 1)) {
            temp = phased[k, sub[k,] == 1]
            temp = temp[temp != 1]
            if (length(temp) > 0) {
              het = getmode(temp)
            } else{
              het = '1|0'
            }
            phased[k, i] = het
          }
          #phased[sub[,i] == 1, i] = "1|0"
          phased[sub[,i] == 0, i] = "0|0"
          phased[sub[,i] == 2, i] = "1|1"
          next
        }
        
        pick <- sample(1:length(potential), 1)
        
        orig <- list(sub[,i] - unlist(potential[pick]))
        if(orig %in% haplo){
          phased[,i] = paste(unlist(orig), unlist(potential[pick]), sep = "|")
          haplo <- c(haplo, potential[pick])
        } else{
          unphased <- c(unphased, i)
        }
        
        #unphased <- c(unphased, i)
      }
    }
    
  }
  
  if (length(unphased) > 0) {
    return(NULL)
  } else{
    return(phased)
  }
  
}

data <- read.table("test_data_masked.txt", header = F, stringsAsFactors = F)

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

test <- data

for (i in 1:nrow(test)) {
  guess = getmode(as.numeric(test[i, test[i,] %in% c(0,2)]))
  
  test[i, test[i, ] == '*'] = guess
  
}

data <- test

rm(test)

data <- data %>%
  mutate_all(as.numeric)

subjects <- c()
for (i in 1:ncol(data)) {
  subjects <- c(subjects, paste0("indiv", i))
}

names(data) <- subjects

totalsnps = nrow(data)
  
all_phased <- c()

i = 1

full_window = 10
window = full_window
use_unique = FALSE

while(i <= totalsnps){
  
  if (window == 0) {
    
    if(use_unique){
      break
    } else{
      window = full_window
      use_unique = TRUE
    }
    
  }
  
  front = i
  end = i + window
  
  if (end > totalsnps) {
    end = totalsnps
  }
  
  sub <- data[front:end,]
  
  df <- data.frame(numhet = colSums(sub==1), index = seq(1:ncol(data))) %>%
    arrange(numhet)
  
  homozyg <- df %>% filter(numhet == 0)
  
  heterozyg <- df %>% filter(numhet > 0)
  
  if (nrow(homozyg) > 0) {
    
    if (use_unique) {
      phased = fix_for_unique_snp(sub, heterozyg, homozyg)
    } else{
      phased = try_to_phase(sub, heterozyg, homozyg)
    }
    
    if (length(phased) == 0) {
      window = window - 1
    } else{
      all_phased <- bind_rows(all_phased, phased)
      i = end+1
      window = full_window
      use_unique = FALSE
    }
  } else{
    window = window - 1
  }
  
}

if (any(all_phased == 1)) {
  print("Error: Not completely phased")
}

done <- all_phased

cols <- names(done)

for(i in 1:length(cols)){
  
  name <- cols[i]
  first <- paste0(name, "A")
  second <- paste0(name, "B")
  
  done <- done %>% separate(!!name, into = c(first, second))
  
}

done <- done %>%
  mutate_all(as.numeric)

filename <- paste0("test_data_sol.txt")

write.table(done, file = filename, row.names = FALSE, col.names = FALSE)



