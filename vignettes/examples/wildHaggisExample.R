
library(momentuHMM)
library(ggplot2)
theme_set(theme_bw())

# Get haggis data from moveHMM
raw <- moveHMM::haggis_data
data <- prepData(data = raw, type = "UTM", covNames = "slope")

# Fit 2-state model with quadratic effect of slope
Par0 <- list(step = c(1, 5, 1, 5), angle = c(pi, 0, 1, 3))
dists <- list(step = "gamma", angle = "vm")
mod <- fitHMM(data = data, 
              nbStates = 2, 
              dist = dists, 
              Par0 = Par0, 
              estAngleMean = list(angle = TRUE),
              formula = ~ slope + I(slope ^ 2))

# Get plotting data for transition probabilities
plot_data <- plot(mod, 
                  plotCI = TRUE, 
                  plotTracks = FALSE, 
                  ask = FALSE, 
                  return = TRUE)
tpm_data_list <- plot_data$estimates$beta$slope

# Add transition probability name as column
tpm_data_list <- lapply(1:length(tpm_data_list), function(i) 
    cbind(tpm_data_list[[i]], names(tpm_data_list)[i]))

# Combine into a single data frame for plotting
tpm_data <- do.call(rbind, tpm_data_list)
colnames(tpm_data)[6] <- "name"

# Create plot of all transition probabilities
p <- ggplot(tpm_data, aes(slope, est)) +
    geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2) +
    geom_line() +
    facet_wrap("name", nrow = 2) +
    ylim(c(0, 1)) +
    labs(y = "transition probability")

ggsave("../plot_wildHaggis.pdf", plot = p, width = 5, height = 4)
