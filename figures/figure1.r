
# ---- Load required libraries ----
library(ape)           # for reading and manipulating phylogenetic trees
library(ggplot2)       # core plotting package
library(ggtree)        # for visualizing phylogenetic trees with ggplot2
library(ggtreeExtra)   # for adding extra layers (like scatterplots) to ggtree plots
library(dplyr)         # for data wrangling and transformation
library(purrr)
library(stringr)







# ---- 2. Load haploid chromosome data from CSV files ----
file_list <- list.files(path = "../data/chrome/", pattern = "\\.csv$", full.names = TRUE)

# Initialize list to store data
all_data <- list()
## Reading the data
for (file in file_list) {
  df <- read.csv(file)
  if (!"haploid" %in% names(df)) {
    cat("Missing 'haploid' column in:", basename(file), "\n")
    next
  }
  file_name <- tools::file_path_sans_ext(basename(file))
  temp_df <- data.frame(clade = file_name, haploid = df$haploid)
  all_data[[length(all_data) + 1]] <- temp_df
}
dat <- do.call(rbind, all_data)
rm(list=ls()[-2])

# ---- 1. Load tree ----
tree <- read.tree(text="(((Chondrichthyes:462.4,((((Anabantaria:103.76124,(Cichlidae:91.7483,Nothobranchiidae:91.7483)14:12.01294)13:8.47282,Gobiidae:112.23406)25:111.7836,((Siluriformes:114.41554,Characiformes:114.41554)37:27.69607,Cyprinidae:142.11161)36:81.90605)35:204.98234,((Caudata:291.41676,Anura:291.41676)34:60.26978,((((Carnivora:76.0,Cetaceans:76.0)43:5.02605,Chiroptera:81.02605)33:12.97395,((Cricetidae:26.17238,Muridae:26.17238)51:61.02762,Primates:87.2)50:6.8)49:224.95,((Gekkonidae:189.28994,(Scincidae:173.54365,(Serpentes:161.03541,Iguania:161.03541)57:12.50824)60:15.74629)56:90.54156,((Accipitriformes:70.41822,Passeriformes:70.41822)48:190.95178,Testudines:261.37)63:18.4615)47:39.1185)68:32.73654)73:77.31346)72:33.4)84:245.2881164,((Araneae:397.0,Scorpiones:397.0)87:166.4032059,(((Blattodea:299.2948,Orthoptera:299.2948)83:83.97117038,(Hymenoptera:343.6917304,(((Dytiscidae:210.0767,Carabidae:210.0767)82:72.54179062,((Coccinellidae:219.2270108,(((Chrysomelidae:191.9352563,(Curculionidae:97.6957875,Scolytidae:97.6957875):94.23946875)80:10.34233208,Cerambycidae:202.2775883):13.76914042,(Elateridae:108.81373,Tenebrionidae:108.81373):107.2329988)79:3.180282083)94:37.53032417,Buprestidae:256.757335)92:25.86115562)78:50.24717588,(Drosophilidae:299.39011,Lepidoptera:299.39011)98:33.4755565)97:10.82606392)77:39.57423997)71:17.81903728,Odonata:401.0850077)67:162.3181983)66:144.2849105)46:890.7934047,(((((((Solanaceae:93.90774,Rubiaceae:93.90774)32:17.25109,Asteraceae:111.15883)106:13.96077,(Brassicaceae:121.36364,(Passiflora:117.06486,Fabaceae:117.06486)104:4.29878)103:3.75596)31:34.49316,(Magnoliaceae:142.05,(Orchidaceae:124.50176,Liliaceae:124.50176)30:17.5555):17.5555)111:122.63268,Gymnospermae:282.235):122.63268,Polypodiaceae:404.87812)29:46.50214,Embryophytes:451.38026)28:1147.101261);")
# make all tip labels lowercase for consistency with data
tree$tip.label <- tolower(tree$tip.label) 
plot(tree, cex=.8)

# Fix the haploid column by parsing multiple formats (commas, dashes, ranges)
dat_clean <- dat %>%
  mutate(
    haploid = str_trim(haploid),
    haploid = map_dbl(haploid, ~ {
      vals <- str_split(.x, ",|–|-")[[1]] |> str_trim()
      nums <- suppressWarnings(as.numeric(vals))
      nums <- nums[!is.na(nums) & nums > 0]
      if (length(nums) == 0) return(NA_real_)
      if (length(nums) == 1) return(nums)
      if (length(nums) == 2 && diff(nums) > 0) return(runif(1, nums[1], nums[2]))
      return(sample(nums, 1))  # fallback: sample any valid value
    })
  ) %>%
  filter(is.finite(haploid), haploid > 0)
dat <- dat_clean
rm(dat_clean)

# ---- 3. Match data to tree tips and log-transform haploid counts ----
final_df <- dat %>%
  filter(clade %in% tree$tip.label) %>%                           # keep only clades that exist in the tree
  mutate(
    clade = factor(clade, levels = tree$tip.label),              # preserve tree tip order for plotting
    log_haploid = log10(haploid)                                 # log-transform haploid counts for better visualization
  )

# ---- 4. Calculate data counts ----
counts <- final_df %>%
  count(clade, name = "counts")                                 # count number of species per clade (species richness)

# TODO Delete if at all possible
# Merge richness back into main dataframe
#final_df <- final_df %>% left_join(counts, by = "clade")

# ---- 5. Create baseline and label dataframes for visualization ----
x_min <- min(final_df$log_haploid)                                # minimum log haploid value (for baseline start)
x_max <- max(final_df$log_haploid)                                # maximum log haploid value (for baseline end)

# Baseline segments for each clade (horizontal reference lines)
baseline_df <- data.frame(
  label = tree$tip.label,                                         # match tree tip order
  y_pos = seq_along(tree$tip.label)                               # numeric y-position for each clade
)

# Create labels with species richness included (e.g., "rosales (n=25)")
#label_df <- richness_df %>%
#  mutate(label_text = paste0(label, " (n=", richness, ")"))

# ---- 6. Build tree and add layered visualization ----
p <- ggtree(tree) +                                               # basic phylogeny plot
  theme_tree() +                                                  # clean tree theme
  ggtitle("Haploid Chromosome Number Distributions per Clade") +  # main title
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

# Add beeswarm points, baselines, and clade labels
p_final <- p +
  # Horizontal baseline lines per clade for visual reference
  geom_facet(
    data = baseline_df,
    panel = "Haploid counts",
    mapping = aes(y = y_pos, yend = y_pos, x = x_min, xend = x_max),
    geom = geom_segment,
    linewidth = 0.35,
    color = "grey55",
    alpha = 0.6
  ) +
  # Beeswarm-style points showing haploid chromosome counts for species in each clade
  geom_facet(
    data = final_df,
    panel = "Haploid counts",
    mapping = aes(x = log_haploid),#, color = richness),
    geom = geom_point,
    position = position_jitter(width = 0.2, height = 0.2),          # jitter points horizontally to reduce overlap
    size = 1.2,
    alpha = 0.1                                                  # transparency so density patterns emerge
  ) +
  # Add clade labels with species richness
  geom_facet(
    data = label_df,
    panel = "Clade labels",
    mapping = aes(x = 0, label = label_text),
    geom = geom_text,
    hjust = 0,                                                   # left-align text
    size = 3,
    fontface = "italic"
  ) +
  # Color points by species richness using a viridis "plasma" gradient on a log scale
  scale_color_gradientn(
    colors = plasma(100),
    name = "Species richness",
    trans = "log10",                                             # log scale for wide range of richness values
    breaks = c(1, 10, 100, 1000),
    labels = c("1", "10", "100", "1000+")
  ) +
  # Use unlogged chromosome counts for x-axis labels
  scale_x_continuous(
    breaks = log10(c(2, 5, 10, 20, 50, 100, 200)),
    labels = c(2, 5, 10, 20, 50, 100, 200)
  ) +
  labs(x = "Haploid chromosome number (log10)") +                # x-axis label
  theme(
    legend.position = "right",                                   # place legend on the right
    panel.grid = element_blank(),                                # remove background gridlines
    axis.title.x = element_text(size = 12, margin = margin(t = 6)),
    axis.text.x  = element_text(size = 10)
  ) +
  expand_limits(x = x_max + 0.2)                                 # add a bit of padding on the x-axis

# ---- 7. Render the final plot ----
print(p_final)








