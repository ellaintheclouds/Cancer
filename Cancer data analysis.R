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
ex1_data$proportion_freq_aa <- ex1_data$n_mut_at_freq/ex1_data$n_mut

# The proportion of truncating mutations, for each gene
ex1_data$proportion_trunc_mut <- ex1_data$n_truncating_mut/ex1_data$n_mut

# Plot----------
# Scatterplot of these two sets of values plotted against each other
ex1_plot <- 
  ggplot(ex1_data, aes(x = proportion_freq_aa, y = proportion_trunc_mut)) +
  geom_point(aes(colour  = gene_name), size = 2) +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0))) + 
  scale_y_continuous(limits = c(0, 0.25), expand = expansion(mult = c(0, 0))) + 
  xlab("Proportion of Total Mutations Accounted for by 
       Most Frequently Mutated Amino Acid") + 
  ylab("Proportion of Truncating Mutations") +
  labs(title = "Rare Cancer Mutational Profile: Proportion of 
       Truncations vs Mutational 'Hot-Spots'", colour = "Gene") + 
  scale_colour_brewer(palette = "Paired") + 
  theme_bw() +
  theme(plot.title = element_text(size=12, face="bold", 
                                  margin = margin(10, 0, 10, 0)))

# Save
ggsave(ex1_plot, filename = "Output/Experiment 1 scatterplot.png", 
       width = 4.5, height = 4)

# Interactive plot
ggplotly(ex1_plot)


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
ex2_data <- list("A" = ex2a_data, "V" = ex2v_data)

# Calculations required to find log survival----------
for(i in 1:2){
  # Calculating mean
  ex2_data[[i]]["Exp 1 mean",] <- apply(ex2_data[[i]][3:5,], 2, mean)
  ex2_data[[i]]["Exp 2 mean",] <- apply(ex2_data[[i]][6:8,], 2, mean)
  ex2_data[[i]]["Exp 3 mean",] <- apply(ex2_data[[i]][9:11,], 2, mean)
  
  ex2_data[[i]]["mean",] <- apply(ex2_data[[i]][12:14,], 2, mean)
  
  # Efficiency (number of colonies/number of cells plated)
  ex2_data[[i]]["Exp 1 efficiency",] <- ex2_data[[i]]["Exp 1 mean",]/ex2_data[[i]]["cell_count",]
  ex2_data[[i]]["Exp 2 efficiency",] <- ex2_data[[i]]["Exp 2 mean",]/ex2_data[[i]]["cell_count",]
  ex2_data[[i]]["Exp 3 efficiency",] <- ex2_data[[i]]["Exp 3 mean",]/ex2_data[[i]]["cell_count",]
  
  ex2_data[[i]]["efficiency",] <- ex2_data[[i]]["mean",]/ex2_data[[i]]["cell_count",]
  
  # % Survival (efficiency at each dose of radiation/efficiency in untreated)
  ex2_data[[i]]["Exp 1 survival",] <- (ex2_data[[i]]["Exp 1 efficiency",]/ex2_data[[i]]["Exp 1 efficiency", 1])*100
  ex2_data[[i]]["Exp 2 survival",] <- (ex2_data[[i]]["Exp 2 efficiency",]/ex2_data[[i]]["Exp 2 efficiency", 1])*100
  ex2_data[[i]]["Exp 3 survival",] <- (ex2_data[[i]]["Exp 3 efficiency",]/ex2_data[[i]]["Exp 3 efficiency", 1])*100
  
  ex2_data[[i]]["survival",] <- (ex2_data[[i]]["efficiency",]/ex2_data[[i]]["efficiency", 1])*100
  
  
  # Standard error of % survival (for SE bars)
  ex2_data[[i]]["sd",] <- apply(ex2_data[[i]][20:22,], 2, sd)
}

View(ex2_data[["A"]])
View(ex2_data[["V"]])

# Plot----------
plot_df <- data.frame("dose" = c(0, 1, 5, 10, 20),
                      "a_cells" = as.numeric(ex2_data[["A"]]["survival",]),
                      "v_cells" = as.numeric(ex2_data[["V"]]["survival",]),
                      "a_sd" = as.numeric(ex2_data[["A"]]["sd",]),
                      "v_sd" = as.numeric(ex2_data[["V"]]["sd",]))

ex2_plot <- ggplot(plot_df) + 
  geom_line(aes(x = dose, y = a_cells, colour = "a_cells")) + 
  geom_point(aes(x = dose, y = a_cells, colour = "a_cells"), size = 2) +
  geom_errorbar(aes(x = dose, ymin = a_cells - a_sd, ymax = a_cells + a_sd, width = 0.4), colour = "darkblue") +
  
  geom_line(aes(x = dose, y = v_cells, colour = "v_cells")) + 
  geom_point(aes(x = dose, y = v_cells, colour = "v_cells"), size = 2) +
  geom_errorbar(aes(x = dose, ymin = v_cells - v_sd, ymax = v_cells + v_sd, width = 0.4), colour = "magenta4") +
  
  scale_color_manual(name = "Cell Line", labels = c("B353-A", "B353-V"), 
                     values = c("a_cells" = "cornflowerblue", 
                                 "v_cells" = "plum")) + 
  xlab("Dose of Ionising Mutation (Gy)") + ylab("Clonogenic Survival (%)") + 
  labs(title = "Survival of B353-V and B353-A Cell Lines
       Exposed to Ionising Radiation") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.005))) + 
  scale_y_continuous(trans="log10", limits = c(NA, 100), 
                     expand = expansion(mult = c(0, 0))) + 
  annotation_logticks(sides = "l") +
  theme_bw() + 
  theme(plot.title = element_text(size=12, face="bold", 
                                  margin = margin(10, 0, 10, 0)))

# Save
ggsave(ex2_plot, filename = "Output/Experiment 2 survival rplot.png", 
       width = 5, height = 4)

# Interactive Plot
ggplotly(ex2_plot)
