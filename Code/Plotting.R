###############################################################################
#                                                                     April '19
#          Discrete Choice Stan Model - RBT Prey Selection
#                             Turbidity
#
#  Notes:
#  *
#
###############################################################################
library(ggmcmc)
library(gridExtra)
library(ggthemes)
library(fishR)

setwd("U:/Desktop/Fish_Git/DCTurbidity/working_runs/")

name.key = data.frame(num = c(1:7),
                      name = c("Snails", "Gammarus",
                               "Black Fly Adult", "Black Fly Larva",
                               "Midge Adult", "Midge Larva","Worms"),
                      name2 = c("Snails", "Gammarus",
                                "Black Fly \nAdult / Pupae", "Black Fly \nLarva",
                                "Midge \nAdult / Pupae", "Midge \nLarva","Worms"))

# all these look very similar
load("Model_Turb_v1_Length_more_parms_500iter.RData")
load("Model_Turb_v1_Width_more_parms_500iter_tree_prob.RData")
load("Model_Turb_v1_Area_more_parms_500iter.RData")
load("Model_Turb_v1_Mass_more_parms_500iter.RData")
fit.1 = fit

#-----------------------------------------------------------------------------#
# the beta size and beta sp effects... run this twice to get both terms

all = bayes_summary(fit.1, par.name = "beta_sp_turb")
all = bayes_summary(fit.1, par.name = "beta_sz_turb")
all$id = as.numeric(gsub("[^0-9]", "", as.character(all$Parameter)))
all$taxa = name.key[match(all$id, name.key$num),3]

#-----------------------------------------------------------------------------#
# verson for the pub...

windows(xpos = 25, record = T, height = 24, width = 30)

p = ggplot(all2, aes(y = my.mean, x = taxa)) +
  geom_hline(yintercept = 0, color = "gray", linetype = 2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, size = 1) +
  # labs(y = "Slopes for Prey \nSelection - Turbidity Relationship \n(+ - 95% CRI)", x = "") +
  labs(y = "Slopes for Size \nEffect - Turbidity Relationship \n(+ - 95% CRI)", x = "") +
  # annotate("text", x = .6, y = 1, label = "A", size = 6)
  annotate("text", x = .6, y = .25, label = "B", size = 6)
p

my_theme =  theme(axis.title.x = element_text(size = 18, vjust = -.1),
                  axis.title.y = element_text(size = 18, vjust = 1),
                  # axis.text.x = element_blank(),
                  axis.text.x = element_text(size = 16, colour = "black"),
                  axis.text.y = element_text(size = 16, colour = "black"),
                  panel.background = element_rect(fill = "white"),
                  panel.grid.minor = element_line(colour = "white"),
                  panel.grid.major = element_line(colour = "white"),
                  panel.border = element_rect(colour = "black", fill = NA),
                  # panel.spacing = unit(.8, "lines"),
                  # plot.margin = unit(c(.8,.8,0,.5), "lines"),
                  plot.margin = unit(c(0,.8,.5,.5), "lines"),
                  title = element_text(size = 18),
                  legend.text = element_text(size = 16),
                  legend.key = element_rect(fill = "white"),
                  legend.title = element_text(size = 16))
g = p + my_theme
g

# g1 = g
# g2 = g
grid.arrange(g1, g2, nrow = 2, heights = c(1/2, 1/2))
#-----------------------------------------------------------------------------#







