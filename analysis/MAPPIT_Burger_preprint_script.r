setwd()

library(dplyr)
library(stringr)
library(readr)
library(qgraph)
library(data.table)
library(yarrr)
library(R.utils)
library(ggplot2)
library(gridExtra)

######################################################################################
#################################### Settings ########################################
######################################################################################

# Identify which items have been present in X% of the assessments
incl <- 1/3 # Inclusion criterion for symptoms (how often do they need to be included?)
m <- 1000 # Number of bootstraps

# Filter rows of participants who have at least 20 observations
nobs <- 20
cut.off <- .7 # Cut-off for evaluating satisfactory stability

######################################################################################
################################## Data handling #####################################
######################################################################################

### Read data
data.raw <- read.csv2("fulldata.csv")
col.data <- read.csv2("data_for_cols2.csv")
phq_data <- read.csv2("phq9.csv")[,c(1,11)]

data <- data.raw
colnames(data) <- colnames(col.data)

# Vector of participants with PHQ9 of at least 10 at beginning
phq <- phq_data %>%
  filter(PHQ9 < 10) %>%
  select(ID)

# Vector of participants who have a stable core of at least 3 nodes
core <- numeric()
ID <- unique(data$ID)
for(i in 1:length(unique(data$ID))){
  data_current <- data[data$ID == ID[i],]
  select_data <- select(data_current, starts_with("select_"))
  select <- names(which(colSums(select_data[,-1], na.rm = TRUE) >= nrow(select_data)*incl))
  if(length(select) < 3){
    core <- c(core, ID[i])
  }
}

# Filter out people with PHQ9 < 10, stable core < 3, and less than 20 days of observations
data <- data %>%
  group_by(ID) %>%
  filter(n() >= nobs) %>%
  filter(!(ID %in% unique(c(unlist(phq), "Plina", "wert", "1992", "0078", core))))

ID <- unique(data$ID)
length(ID) # 21 participants have at least 20 observations and a PHQ of at least 10

# Labels of all symptoms
labels <- gsub("_C_1","", colnames(select(data, ends_with("_C_1")))[-1])


######################################################################################
############################### Create objects #######################################
######################################################################################

combined.matrices <- list()
combined.boot <- data.frame()
df.combined_full <- data.frame()
df.combined <- data.frame()

cause.matrices <- list()
cause.boot <- data.frame()
df.cause_full <- data.frame()
df.cause <- data.frame()

effect.matrices <- list()
effect.boot <- data.frame()
df.effect_full <- data.frame()
df.effect <- data.frame()

missIDs <- numeric()
comb.plots <-list()

# Code stability at day X for each participant
stability.at.X.cause <- numeric(length(ID))
stability.at.X.effect <- numeric(length(ID))
stability.at.X.combined <- numeric(length(ID))

# Code number of days it took to reach X stability
days.until.X.stability.cause <- numeric(length(ID))
days.until.X.stability.effect <- numeric(length(ID))
days.until.X.stability.combined <- numeric(length(ID))

# Code agreement between cause and effect data
agreement <- numeric()
agreement_2 <- numeric()
core.size <- numeric()

# Slopes for causal drift
slopes <- numeric()
  
######################################################################################
################################ Start bootstrap #####################################
######################################################################################

# Select Participant for Figure 1
# i <- 2
# Select Participant with good stability: Participant #15
# i <- 15
# Select Participant with bad stability: Participant #9
# i <- 9


# Loop over each participant
  for(i in 1:length(unique(data$ID))){
    
    data_current <- data[data$ID == ID[i],]
    
    select_data <- select(data_current, starts_with("select_"))
    
    # Vector of selected variables
    select <- names(which(colSums(select_data[,-1], na.rm = TRUE) >= nrow(select_data)*incl))
    
    # Indices of variables to use
    vec <- parse_number(select)
    core.size[i] <- length(vec)
      
      ## Cause matrix
      # Columns to select, encoding causes for the selected variables
      causes <- paste0(rep(labels[vec], each = length(vec)), "_C_", vec)
      effects <- paste0(rep(labels[vec], length(vec)), "_E_", rep(vec, each = length(vec)))
      
      # Cause + Effect data
      cause.data <- data_current[,causes]
      cause.data[is.na(cause.data)] <- 0
      names(cause.data) <- paste(rep(labels[vec], each = length(vec)), "CAUSED BY", rep(labels[vec], length(vec)))
      
      effect.data <- data_current[,effects]
      effect.data[is.na(effect.data)] <- 0
      names(effect.data) <- paste(rep(labels[vec], each = length(vec)), "CAUSED BY", rep(labels[vec], length(vec)))
      
      # Agreement between cause and effect data
      agreement[i] <- sum(cause.data == effect.data) / (nrow(cause.data) * ncol(cause.data))
      
      # Combine the datasets
      combined.data <- matrix(NA, ncol = ncol(cause.data), nrow = nrow(cause.data))
      
      for(j in 1:nrow(combined.data)){
      combined.data[j,] <- as.numeric(unlist(cause.data[j,])) + as.numeric(unlist(effect.data[j,]))
      }
      
      rownames(combined.data) <- NULL
      colnames(combined.data) <- colnames(cause.data)
      
      agreement_2[i] <- sum(c(combined.data) == 2) / (sum(c(combined.data) == 2) + sum(c(combined.data) == 1))
      
      combined.data <- ifelse(combined.data != 0, 1, 0)
      
      # This uses *absolute* frequency of cause.
      combined.matrix <- matrix(colSums(combined.data == 1), length(vec), length(vec)) 
      colnames(combined.matrix) <- rownames(combined.matrix) <- labels[vec]
      
      combined.matrices[[i]] <- combined.matrix
      
      # write.csv(combined.data, paste0("combined.data_p", i, ".csv"))
      # write.csv(combined.matrix, paste0("combined.matrix_p", i, ".csv"))
      
#### Bootstrap for combined data ####
      
      # Calculate causal drift
      drift.comb <- combn(1:nrow(combined.data), 2)
      cor.comb <- numeric()
      dist.comb <- numeric()
      
      for(k in 1:ncol(drift.comb)){
        # Matrices of the two currently selected days
        m.comb1 <- matrix(combined.data[drift.comb[1,k],], length(vec), length(vec))
        m.comb2 <- matrix(combined.data[drift.comb[2,k],], length(vec), length(vec))
        
        # Correlation between matrices
        cor.comb[k] <- cor(m.comb1, m.comb2)
        
        # Distance between time points
        dist.comb[k] <- abs(drift.comb[1,k] - drift.comb[2,k])
      }
      
      df.comb <- data.frame(correlation = unlist(cor.comb),
                            distance = dist.comb)
      
      # Store slope
      if(length(na.omit(unlist(cor.comb))) > 0){
        res <- lm(unlist(cor.comb) ~ dist.comb)
        slopes[i] <- res$coefficients[2]
      }
    
      comb.plots[[i]] <- ggplot(df.comb, aes(x=distance, y=correlation)) + 
        geom_point() +
        geom_smooth(method=lm) +
        theme(panel.background = element_blank()) +
        ggtitle(paste0("Participant ", i))
      
      # Actual bootstrap
      ind <- nrow(combined.data)
      cor.vec.combined <- numeric(m)
      count <- 0
      boot.mat.combined <- matrix(NA, 
                                  m,
                                  length(seq(2, nrow(combined.data), by = 2)))
      boot.vec.combined <- numeric()
      lower.combined <- numeric()
      upper.combined <- numeric()
      
      set.seed(123)
      for(a in seq(2, nrow(combined.data), by = 2)){
        
        count <- count + 1
        
        for(b in 1:m){
          
          # Create index-vector for two sub-samples
          v <- sample(1:ind, size = a, replace = FALSE)
          v1 <- sample(v, size = a/2, replace = FALSE)
          v2 <- v[!(v %in% v1)]
          
          # Construct matrix from both samples
          # Sample 1:
          if(length(v1) == 1){
            m1 <- matrix(as.matrix(combined.data)[v1,], length(vec), length(vec))
          } else {
            m1 <- matrix(colSums(as.matrix(combined.data)[v1,] == 1), length(vec), length(vec))
          }
          colnames(m1) <- rownames(m1) <- labels[vec]
          
          # Sample 2:
          if(length(v2) == 1){
            m2 <- matrix(as.matrix(combined.data)[v2,], length(vec), length(vec))
          } else {
            m2 <- matrix(colSums(as.matrix(combined.data)[v2,] == 1), length(vec), length(vec))
          }
          colnames(m2) <- rownames(m2) <- labels[vec]
          
          # Store correlation between the two matrices
          r <- cor.test(m1, m2, method = "spearman")
          cor.vec.combined[b] <- r$estimate
          
        }
        
        print(mean(cor.vec.combined, na.rm = TRUE)) # Print mean
        
        # Store bootstrap objects
        boot.mat.combined[,count] <- cor.vec.combined
        boot.vec.combined[count] <- mean(cor.vec.combined, na.rm = TRUE)
        lower.combined[count] <- boot.vec.combined[count] - 1.96 * (sd(cor.vec.combined, na.rm = TRUE) / sqrt(m))
        upper.combined[count] <- boot.vec.combined[count] + 1.96 * (sd(cor.vec.combined, na.rm = TRUE) / sqrt(m))
      }
      
      # Save stability characteristics
      stability.at.X.combined[i] <- boot.vec.combined[nobs/2] # stability at nobs (here: 20)
      days.until.X.stability.combined[i] <- which(boot.vec.combined >= cut.off)[1] # days until .7 stability is reached
      
#### Bootstrap for cause data ####
      ind <- nrow(cause.data)
      cor.vec.cause <- numeric(m)
      count <- 0
      boot.mat.cause <- matrix(NA, 
                               m,
                               length(seq(2, nrow(cause.data), by = 2)))
      boot.vec.cause <- numeric()
      lower.cause <- numeric()
      upper.cause <- numeric()
      
      set.seed(123)
      for(a in seq(2, nrow(cause.data), by = 2)){
        
        count <- count + 1
        
        for(b in 1:m){
          
          # Create index-vector for two sub-samples
          v <- sample(1:ind, size = a, replace = FALSE)
          v1 <- sample(v, size = a/2, replace = FALSE)
          v2 <- v[!(v %in% v1)]
          
          # Construct matrix from both samples
          # Sample 1:
          m1 <- matrix(colSums(cause.data[v1,] == 1), length(vec), length(vec))
          colnames(m1) <- rownames(m1) <- labels[vec]
          
          # Sample 2:
          m2 <- matrix(colSums(cause.data[v2,] == 1), length(vec), length(vec))
          colnames(m2) <- rownames(m2) <- labels[vec]
          
          # Store correlation between the two matrices
          r <- cor.test(m1, m2, method = "spearman")
          cor.vec.cause[b] <- r$estimate
          
        }
        
        print(mean(cor.vec.cause, na.rm = TRUE)) # Print mean
        
        # Store bootstrap objects
        boot.mat.cause[,count] <- cor.vec.cause
        boot.vec.cause[count] <- mean(cor.vec.cause, na.rm = TRUE)
        lower.cause[count] <- boot.vec.cause[count] - 1.96 * (sd(cor.vec.cause, na.rm = TRUE) / sqrt(m))
        upper.cause[count] <- boot.vec.cause[count] + 1.96 * (sd(cor.vec.cause, na.rm = TRUE) / sqrt(m))
        
      }
      
      # This uses *absolute* frequency of cause.
      cause.matrix <- matrix(colSums(cause.data == 1), length(vec), length(vec)) 
      colnames(cause.matrix) <- rownames(cause.matrix) <- labels[vec]
      
      # Save stability characteristics
      stability.at.X.cause[i] <- boot.vec.cause[nobs/2] # stability at nobs (here: 20)
      days.until.X.stability.cause[i] <- which(boot.vec.cause >= cut.off)[1] # days until .7 stability is reached
      
      
#### Bootstrap for effect data ####
      ind <- nrow(effect.data)
      cor.vec.effect <- numeric(m)
      count <- 0
      boot.mat.effect <- matrix(NA, 
                                m,
                                length(seq(2, nrow(effect.data), by = 2)))
      boot.vec.effect <- numeric()
      lower.effect <- numeric()
      upper.effect <- numeric()
      
      set.seed(123)
      for(a in seq(2, nrow(effect.data), by = 2)){
        
        count <- count + 1
        
        for(b in 1:m){
          
          # Create index-vector for two sub-samples
          v <- sample(1:ind, size = a, replace = FALSE)
          v1 <- sample(v, size = a/2, replace = FALSE)
          v2 <- v[!(v %in% v1)]
          
          # Construct matrix from both samples
          # Sample 1:
          m1 <- matrix(colSums(effect.data[v1,] == 1), length(vec), length(vec))
          colnames(m1) <- rownames(m1) <- labels[vec]
          
          # Sample 2:
          m2 <- matrix(colSums(effect.data[v2,] == 1), length(vec), length(vec))
          colnames(m2) <- rownames(m2) <- labels[vec]
          
          # Store correlation between the two matrices
          r <- cor.test(m1, m2, method = "spearman")
          cor.vec.effect[b] <- r$estimate
          
        }
        
        print(mean(cor.vec.effect, na.rm = TRUE)) # Print mean
        
        # Store bootstrap objects
        boot.mat.effect[,count] <- cor.vec.effect
        boot.vec.effect[count] <- mean(cor.vec.effect, na.rm = TRUE)
        lower.effect[count] <- boot.vec.effect[count] - 1.96 * (sd(cor.vec.effect, na.rm = TRUE) / sqrt(m))
        upper.effect[count] <- boot.vec.effect[count] + 1.96 * (sd(cor.vec.effect, na.rm = TRUE) / sqrt(m))
        
      }
      
      # This uses *absolute* frequency of cause.
      effect.matrix <- matrix(colSums(effect.data == 1), length(vec), length(vec)) 
      colnames(effect.matrix) <- rownames(effect.matrix) <- labels[vec]
      
      # Save stability characteristics
      stability.at.X.effect[i] <- boot.vec.effect[nobs/2] # stability at nobs (here: 20)
      days.until.X.stability.effect[i] <- which(boot.vec.effect >= cut.off)[1] # days until .7 stability is reached
      
      # Save all results
      # write.csv(boot.mat.cause, paste0("p_", i, "_cause.boot.csv"))
      # write.csv(boot.mat.effect, paste0("p_", i, "_effect.boot.csv"))
      # write.csv(boot.mat.combined, paste0("p_", i, "_combined.boot.csv"))

  }


##### Results section
phq.vec <- c(22,
             13, 
             20,
             15, 
             11, 
             18, 
             13, 
             18, 
             14, 
             18, 
             17, 
             21, 
             16, 
             15, 
             19, 
             12, 
             14, 
             12, 
             24, 
             16
)

mean(phq.vec)
mean(agreement_2, na.rm = TRUE)
min(agreement_2, na.rm = TRUE)
max(agreement_2, na.rm = TRUE)

# Save results
results.data.frame <- data.frame(#ID = ID,
                                 core.size = core.size,
                                 stability.at.X.cause = stability.at.X.cause,
                                 stability.at.X.effect = stability.at.X.effect,
                                 stability.at.X.combined = stability.at.X.combined,
                                 days.until.X.stability.cause = days.until.X.stability.cause,
                                 days.until.X.stability.effect = days.until.X.stability.effect,
                                 days.until.X.stability.combined = days.until.X.stability.combined,
                                 phq = phq.vec,
                                 single.multiple = c(rep(0,10), rep(1,10)))

# write.csv(results.data.frame, "results.full.csv")

# FIGURE 1 (example network)
setwd()
matrix_p2 <- read.csv2("combined.matrix_p2.csv", sep = ",")[,-1]
rownames(matrix_p2) <- colnames(matrix_p2)

qgraph(matrix_p2,
       labels = c("procr.", "tired", "unfoc."),
       theme = "colorblind",
       vsize = 20,
       mar = rep(5,4))

### Descriptives

## Core size
# Overall 
mean(core.size) # 6.2
sd(core.size) # 3.955
# Single condition 
mean(core.size[1:10]) # 6.4
sd(core.size[1:10]) # 5.103
# Multiple condition 
mean(core.size[11:20]) # 6
sd(core.size[11:20]) # 2.625
t.test(core.size[1:10], core.size[11:20]) # t(13.45) = .22, p = .829

## Stability at Day 20
# Overall combined
mean(stability.at.X.combined) # 0.595
sd(stability.at.X.combined) # 0.201
# Single condition combined
mean(stability.at.X.combined[1:10]) # 0.562
sd(stability.at.X.combined[1:10]) # 0.248
# Multiple condition combined
mean(stability.at.X.combined[11:20]) # 0.629
sd(stability.at.X.combined[11:20]) # 0.146

# Overall cause
mean(stability.at.X.cause) # 0.588
sd(stability.at.X.cause) # 0.221
# Single condition cause
mean(stability.at.X.cause[1:10]) # 0.597
sd(stability.at.X.cause[1:10]) # 0.263
# Multiple condition cause
mean(stability.at.X.cause[11:20]) # 0.580
sd(stability.at.X.cause[11:20]) # 0.183

# Overall effect
mean(stability.at.X.effect) # 0.548
sd(stability.at.X.effect) # 0.212
# Single condition effect
mean(stability.at.X.effect[1:10]) # 0.472
sd(stability.at.X.effect[1:10]) # 0.236
# Multiple condition effect
mean(stability.at.X.effect[11:20]) # 0.623
sd(stability.at.X.effect[11:20]) # 0.163

## Days until .70 stability
# Overall combined
mean(days.until.X.stability.combined*2, na.rm = TRUE) # 15.25
sd(days.until.X.stability.combined*2, na.rm = TRUE) # 3.196
# Single condition combined
mean(days.until.X.stability.combined[1:10]*2, na.rm = TRUE) # 13.333
sd(days.until.X.stability.combined[1:10]*2, na.rm = TRUE) # 2.309
# Multiple condition combined
mean(days.until.X.stability.combined[11:20]*2, na.rm = TRUE) # 16.4
sd(days.until.X.stability.combined[11:20]*2, na.rm = TRUE) # 3.286

# Overall cause
mean(days.until.X.stability.cause*2, na.rm = TRUE) # 14
sd(days.until.X.stability.cause*2, na.rm = TRUE) # 7.659
# Single condition cause
mean(days.until.X.stability.cause[1:10]*2, na.rm = TRUE) # 12.5
sd(days.until.X.stability.cause[1:10]*2, na.rm = TRUE) # 9.983
# Multiple condition cause
mean(days.until.X.stability.cause[11:20]*2, na.rm = TRUE) # 16
sd(days.until.X.stability.cause[11:20]*2, na.rm = TRUE) # 4

# Causal drift plots
# All plots for the appendix
# grid.arrange(comb.plots[[1]], comb.plots[[2]], comb.plots[[3]], comb.plots[[4]],
#              comb.plots[[5]], comb.plots[[6]], comb.plots[[7]], comb.plots[[8]],
#              comb.plots[[9]], comb.plots[[10]], comb.plots[[11]], comb.plots[[12]],
#              comb.plots[[13]], comb.plots[[14]], comb.plots[[15]], comb.plots[[16]],
#              comb.plots[[17]], comb.plots[[18]], comb.plots[[19]], comb.plots[[20]],
#              comb.plots[[21]],
#              nrow = 5)

# Main plots for paper to illustrate 3 cases
grid.arrange(comb.plots[[12]], comb.plots[[7]], comb.plots[[20]],
             nrow = 1)

# Slope descriptives
mean(slopes, na.rm = TRUE) # mean = -.015
min(slopes, na.rm = TRUE) # min = -.123
max(slopes, na.rm = TRUE) # max = .076
sd(slopes, na.rm = TRUE) # SD = .041

### Correlates
cor.df <- data.frame(Core.size = core.size,
                     Stability = stability.at.X.combined,
                     PHQ = phq.vec)
cor.test(cor.df[,2], cor.df[,3])

mean(phq.vec[is.na(days.until.X.stability.combined)])
mean(phq.vec[!is.na(days.until.X.stability.combined)])
