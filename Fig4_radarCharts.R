# getwd()
gc()
setwd("~/Documents/immGen/GAM-clustering/10_32_04_191021_Todorov_32/Clustering/_current/Fig4/")
library(tidyverse)


# load('a.rda')
# anno <- a
anno <- read_csv("~/Documents/immGen/GAM-clustering/10_32_04_191021_Todorov_32/Clustering/_current/Fig4/anno.csv")

load("~/Documents/immGen/GAM-clustering/10_32_04_191021_Todorov_32/Clustering/_current/Fig4/eyegene.rda")
df <- eyegene %>% as_tibble()


# Glimpse on annotation
anno %>%
  count(value) %>%
  arrange(desc(n))

# Calculate per-module statistics
stats <- df %>%
  mutate(module = row_number()) %>%
  gather('sample', 'value', -module) %>%
  group_by(module) %>%
  summarize(max = max(value),
            min = min(value))
stats


###################### FUNCTIONS #######################

plot_radarchart <- function(long) {
  long %>%
    spread(module, value) %>%
    as.data.frame() %>%
    column_to_rownames('sample') %>% {
      # bind in following order: max, min, data, mean
      rbind(stats$max, # + 1, # max(stats$max)
            stats$min, # - 1, # min(stats$min)
            # 0,
            .,
            means)
    } %>% {
      .[, c(1, ncol(.):2)]
    } %>%
    fmsb::radarchart(
      pty = 32,  # point symbol (16 filled circle, 32 nothing)
      pcol = line_color,
      plwd = line_width,
      plty = line_type,
      cglcol = 'gray22', # grid color
      cglwd = 1,       # grid width
      cglty = 2,       # grid type (1 solid, 2 dashed, 3 dotted)
      axistype = 0,    # just 0
      axislabcol = "black", # labels color
      vlcex = 1        # font size magnification for labels
    )
}

plot_radarchart_pdf <- function(long, filename) {
  pdf(filename, width = 4, height = 4)
  plot_radarchart(long)
  dev.off()
}

plot_radarchart_png <- function(long, filename) {
  png(filename)
  plot_radarchart(long)
  dev.off()
}

####################################################



# Choose specific meta-sample
load("~/Documents/immGen/GAM-clustering/10_32_04_191021_Todorov_32/Clustering/annotation_colors.Rda")
load("~/Documents/immGen/GAM-clustering/10_32_04_191021_Todorov_32/Clustering/_current/Fig4/annotation_colors.Rda")
View(annotation_colors$metaSample)
names(annotation_colors$metaSample)


for(metasample in names(annotation_colors$metaSample)){
  print(metasample)

  # Pick specific meta-sample
  filtered <- df %>%
    select(select_vars(names(df),
                       anno %>%
                         filter(value == metasample) %>%
                         pull(samples)))
  # Convert to long format
  long <- filtered %>%
    # mutate(module = str_c('Mod.', row_number())) %>%
    mutate(module = row_number()) %>%
    gather('sample', 'value', -module)
  # Calculate means
  means <- long %>%
    group_by(module) %>%
    summarize(mean = mean(value)) %>%
    pull(mean)

  # Number of samples
  k <- ncol(filtered) %>% print()

  { # Setup plot line color
    color_all <- 'ivory4'  # all samples
    color_mean <- unname(annotation_colors$metaSample[metasample])  # mean of sample # navy
    line_color <- c(rep(color_all, k), color_mean)
  }
  { # Setup plot line width
    width_all <- 2
    width_mean <- 4
    line_width <- c(rep(width_all, k), width_mean)
  }
  { # Setup plot line type (1 solid, 2 dashed, 3 dotted)
    type_all <- 1
    type_mean <- 1
    line_type <- c(rep(type_all, k), type_mean)
  }

  # plot_radarchart(long)
  plot_radarchart_pdf(long, str_interp('spiderplot-${metasample}.pdf'))
  # plot_radarchart_png(long, str_interp('spiderplot-${metasample}.png'))
}
