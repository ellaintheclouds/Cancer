# Set up------------------------------------------------------------------------
#setwd("C:/Users/El Richardson/OneDrive - Lancaster University/Biomedicine/Year 
#      3/353 Cancer/Coursework")
if (!require("ggplot2", quietly = TRUE)) { install.packages("ggplot2") } 
if (!require("plotly", quietly = TRUE)) { install.packages("plotly") } 
if (!require("RColorBrewer", quietly = TRUE)) { install.packages("RColorBrewer") } 

library(ggplot2)
library(plotly)
library(RColorBrewer)


# Experiment 1 -----------------------------------------------------------------
ex1_data <- read.csv("Input/Experiment 1.csv")
View(ex1_data)

# Proportion of all mutations that the most frequently mutated amino acid 
# accounts for
ex1_data$proportion_freq_aa <- NA

for(i in 1:10){
  ex1_data[i, "proportion_freq_aa"] <- 
    ex1_data[i, "n_mut_at_freq"]/ex1_data[i, "n_mut"]
}

# The proportion of truncating mutations, for each gene
for(i in 1:10){
  ex1_data[i, "proportion_trunc_mut"] <- 
    ex1_data[i, "n_truncating_mut"]/ex1_data[i, "n_mut"]
}

# Plot----------
# Scatterplot of these two sets of values plotted against each other
title_label <- "Novel Cancer Mutational Profile: 
Proportion of Truncations vs Mutational 'Hot-Spots'"

ex1_plot <- 
  ggplot(ex1_data, aes(x = proportion_freq_aa, y = proportion_trunc_mut)) +
  geom_point(aes(colour  = gene_name), size = 2) +
  xlab("Proportion of Mutations that Most Frequently Mutated Amino Acid Accounts For") + 
  ylab("Proportion of Truncating Mutations") +
  labs(title = title_label,
       colour = "Gene") + 
  scale_colour_brewer(palette = "Spectral") + 
  theme_bw()

# Interactive plot
ggplotly(ex1_plot)

# Save
ggsave(ex1_plot, filename = "Output/Experiment 1 scatterplot.png", 
       width = 7, height = 6)


# Experiment 2 -----------------------------------------------------------------
# Loading data for A
ex2a_data <- read.csv("Input/Experiment 2 B353-A.csv")
rownames(ex2a_data) <- ex2a_data[,1]
ex2a_data <- ex2a_data[,-1]

# Loading data for V
ex2v_data <- read.csv("Input/Experiment 2 B353-V.csv")
rownames(ex2v_data) <- ex2v_data[,1]
ex2v_data <- ex2v_data[,-1]

# Putting A and V data in to a list
ex2_data <- list(ex2a_data, ex2v_data)

# Calculations required to find log survival----------
for(i in 1:2){
  # Calculating mean
  ex2_data[[i]]["mean",] <- apply(ex2_data[[i]][3:11,], 2, mean)
  
  # Efficiency (number of colonies/number of cells plated)
  ex2_data[[i]]["efficiency",] <- ex2_data[[i]]["mean",]/ex2_data[[i]]["cell_count",]
  
  # % Survival (efficiency at each dose of radiation/efficiency in untreated)
  ex2_data[[i]]["survival",] <- (ex2_data[[i]]["efficiency",]/ex2_data[[i]]["efficiency", 1])*100
  
  # Log % survival
  ex2_data[[i]]["log_survival",] <- log(ex2_data[[i]]["survival",])
}

View(ex2_data[[1]])
View(ex2_data[[2]])

# Plot----------
ex2_title <- "Survival of B353-V and B353-A Cell Lines Exposed to Ionising Radiation"

plot_df <- data.frame("dose" = c(0, 1, 5, 10, 20),
                      "a_cells" = as.numeric(ex2_data[[1]]["log_survival",]),
                      "v_cells" = as.numeric(ex2_data[[2]]["log_survival",]))

ex2_plot <- ggplot(plot_df) + 
  geom_line(aes(x = dose, y = a_cells, colour = "a_cells")) + 
  geom_point(aes(x = dose, y = a_cells, colour = "a_cells"), size = 2) +
  geom_line(aes(x = dose, y = v_cells, colour = "v_cells")) + 
  geom_point(aes(x = dose, y = v_cells, colour = "v_cells"), size = 2) + 
  scale_color_manual(name = "Cell Line", labels = c("B353-V", "B353-A"), 
                     values = c("a_cells" = "cornflowerblue", 
                                 "v_cells" = "deeppink4")) + 
  xlab("Dose of Ionising Mutation (Grays)") + ylab("log(% Survival of Cells)") + 
  labs(title = ex2_title) +
  theme_bw()

# Interactive Plot
ggplotly(ex2_plot)

# Save
ggsave(ex2_plot, filename = "Output/Experiment 2 survival rplot.png", 
       width = 7, height = 6)
