# ============================
# Phylogeny + Beeswarm Points + Clade Labels (colored by higher classification)
# ============================

# ---- 1. Load libraries ----
library(ape)           # read and manipulate phylogenetic trees (e.g., .tree or Newick format)
library(ggplot2)       # general-purpose plotting library
library(ggtree)        # visualize phylogenetic trees using ggplot2 syntax
library(ggtreeExtra)   # add additional panels ("facets") like scatterplots to tree plots
library(dplyr)         # data manipulation (filtering, counting, merging, etc.)
library(stringr)       # string handling (cleaning, trimming, extracting patterns)

# ---- 2. Load haploid chromosome data ----
# Identify all CSV files in the data directory that contain haploid chromosome counts
file_list <- list.files(path = "../data/chrome", pattern = "\\.csv$", full.names = TRUE)

# Read and combine all CSVs into a single data frame
dat <- lapply(file_list, function(file) {
  df <- read.csv(file, stringsAsFactors = FALSE)              # read CSV with character strings preserved
  if (!"haploid" %in% names(df)) return(NULL)                # skip file if it lacks a 'haploid' column
  data.frame(
    clade = tolower(tools::file_path_sans_ext(basename(file))), # use lowercase filename (without extension) as clade name
    haploid = as.character(df$haploid),                         # convert haploid values to character (prevents bind_rows errors)
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()                                            # combine all individual data frames into one

# ---- 3. Load tree ----
# Read the phylogenetic tree from a Newick string and standardize tip labels to lowercase
tree <- read.tree(text = "(((Chondrichthyes:462.4,((((Anabantaria:103.76124,(Cichlidae:91.7483,Nothobranchiidae:91.7483)14:12.01294)13:8.47282,Gobiidae:112.23406)25:111.7836,((Siluriformes:114.41554,Characiformes:114.41554)37:27.69607,Cyprinidae:142.11161)36:81.90605)35:204.98234,((Caudata:291.41676,Anura:291.41676)34:60.26978,((((Carnivora:76.0,Cetaceans:76.0)43:5.02605,Chiroptera:81.02605)33:12.97395,((Cricetidae:26.17238,Muridae:26.17238)51:61.02762,Primates:87.2)50:6.8)49:224.95,((Gekkonidae:189.28994,(Scincidae:173.54365,(Serpentes:161.03541,Iguania:161.03541)57:12.50824)60:15.74629)56:90.54156,((Accipitriformes:70.41822,Passeriformes:70.41822)48:190.95178,Testudines:261.37)63:18.4615)47:39.1185)68:32.73654)73:77.31346)72:33.4)84:245.2881164,((Araneae:397.0,Scorpiones:397.0)87:166.4032059,(((Blattodea:299.2948,Orthoptera:299.2948)83:83.97117038,(Hymenoptera:343.6917304,(((Dytiscidae:210.0767,Carabidae:210.0767)82:72.54179062,((Coccinellidae:219.2270108,(((Chrysomelidae:191.9352563,(Curculionidae:97.6957875,Scolytidae:97.6957875):94.23946875)80:10.34233208,Cerambycidae:202.2775883):13.76914042,(Elateridae:108.81373,Tenebrionidae:108.81373):107.2329988)79:3.180282083)94:37.53032417,Buprestidae:256.757335)92:25.86115562)78:50.24717588,(Drosophilidae:299.39011,Lepidoptera:299.39011)98:33.4755565)97:10.82606392)77:39.57423997)71:17.81903728,Odonata:401.0850077)67:162.3181983)66:144.2849105)46:890.7934047,(((((((Solanaceae:93.90774,Rubiaceae:93.90774)32:17.25109,Asteraceae:111.15883)106:13.96077,(Brassicaceae:121.36364,(Passiflora:117.06486,Fabaceae:117.06486)104:4.29878)103:3.75596)31:34.49316,(Magnoliaceae:142.05,(Orchidaceae:124.50176,Liliaceae:124.50176)30:17.5555):17.5555)111:122.63268,Gymnospermae:282.235):122.63268,Tracheophytes:404.87812)29:46.50214,Bryophytes:451.38026)28:1147.101261);")
tree$tip.label <- tolower(tree$tip.label)                     # ensure tip names match data frame clade names

# ---- 4. Clean haploid values ----
# Extract the first numeric value from each haploid entry and remove invalid rows
dat_clean <- dat %>%
  mutate(haploid = as.numeric(str_extract(haploid, "\\d+"))) %>%  # extract numeric part from text like "10-12"
  filter(!is.na(haploid) & haploid > 0)                           # keep only valid positive numeric values

# ---- 5. Filter to tree tips and log-transform ----
# Restrict dataset to clades present in the tree and calculate log10 haploid values for plotting
final_df <- dat_clean %>%
  filter(clade %in% tree$tip.label) %>%                           # drop clades not present in tree tips
  mutate(log_haploid = log10(haploid))                            # log-transform chromosome numbers for visualization

# ---- 6. Add higher classification ----
# Load classification info and standardize clade names
class_df <- read.csv("../data/higher_class.csv", stringsAsFactors = FALSE)
class_df$Clade <- tolower(str_trim(class_df$Clade))              # clean and lowercase clade names for matching

# Merge classification into the main dataset
final_df <- merge(final_df, class_df, by.x = "clade", by.y = "Clade", all.x = TRUE)
names(final_df)[names(final_df) == "Higher.Classification"] <- "higher_class"  # rename column to simpler label

# ---- 7. Count number of species per clade ----
# Summarize how many data points per clade and build label text (e.g., "muridae (n=34)")
counts <- final_df %>%
  count(clade, name = "count") %>%
  mutate(label_text = paste0(clade, " (n=", count, ")"))

# ---- 8. Build tree ----
# Create a ggtree object as the base layer of the figure
p <- ggtree(tree) + theme_tree()

# ---- 9. Align data order with tree tips ----
# Extract tip order from tree plot and align data frame rows with this order
tip_order <- p$data$label[p$data$isTip]                         # ordered list of tip labels from plotted tree
baseline_df <- data.frame(clade = tree$tip.label,              # baseline reference for each clade
                          y_pos = match(tree$tip.label, tip_order))  # map clade name to its vertical position
final_df$y_pos <- match(final_df$clade, tip_order)             # add y-position for chromosome points
counts$y_pos <- match(counts$clade, tip_order)                 # add y-position for labels

# ---- 10. Define x-axis range ----
# Determine minimum and maximum log haploid values for axis scaling
x_min <- min(final_df$log_haploid)
x_max <- max(final_df$log_haploid) + 0.25                       # add padding so points don't hit edge

# ---- 11. Final plot ----
p_final <- p +
  # Draw baseline line across chromosome range for each clade
  geom_facet(
    data = baseline_df,
    panel = "Haploid counts",
    aes(y = y_pos, yend = y_pos, x = x_min, xend = x_max),
    geom = geom_segment,
    linewidth = 0.35,
    color = "grey55",
    alpha = 0.6
  ) +
  # Plot chromosome points (transparent for density visualization)
  geom_facet(
    data = final_df,
    panel = "Haploid counts",
    aes(x = log_haploid, y = y_pos, color = higher_class),
    geom = geom_point,
    position = position_jitter(width = 0.2, height = 0.2),      # jitter points slightly to reduce overlap
    size = 1.2,
    alpha = 0.1                                                # low opacity for overlapping density
  ) +
  # Add clade labels and species counts (non-italic text)
  geom_facet(
    data = transform(counts, x_text = 0),                      # place text at x = 0
    panel = "Clade labels",
    aes(x = x_text, y = y_pos, label = label_text),
    geom = geom_text,
    size = 2.8,
    hjust = 0
  ) +
  # Define custom color palette for higher classification
  scale_color_manual(
    name = "Higher classification",
    values = c(
      "Chondrichthyes" = "#414487",   # deep indigo
      "Actinopterygii" = "#2A788E",   # teal-blue
      "Amphibia"       = "#22A884",   # bright sea green
      "Mammalia"       = "#7AD151",   # saturated green
      "Reptilia"       = "#F8961E",   # rich orange
      "Arachnida"      = "#D1495B",   # deep rose
      "Insecta"        = "#8E3B9E",   # violet-purple
      "Angiosperm"     = "#1F968B",   # turquoise
      "Gymnosperms"    = "#89C2D9",   # sky blue
      "Tracheophytes"  = "#5E4FA2",   # deep purple-blue
      "Bryophytes"     = "#277DA1"    # rich blue
    )
  ) +
  # Make legend points fully opaque (ignore alpha from data points)
  guides(
    color = guide_legend(
      override.aes = list(alpha = 1, size = 2, stroke = 0.5)
    )
  ) +
  # Format x-axis (log-transformed chromosome number)
  scale_x_continuous(
    breaks = log10(c(2, 5, 10, 20, 50, 100, 200)),            # where tick marks appear
    labels = c(2, 5, 10, 20, 50, 100, 200)                    # display unlogged values
  ) +
  labs(x = "Haploid chromosome number (log10)") +             # x-axis label
  theme(
    legend.position = "right",                                # place legend to the right
    panel.grid = element_blank(),                             # remove background grid
    axis.title.x = element_text(size = 12, margin = margin(t = 6)), # style x-axis title
    axis.text.x = element_text(size = 10)                     # style x-axis tick labels
  ) +
  expand_limits(x = x_max + 0.2)                              # extra space on x-axis

# ---- 12. Render final plot ----
print(p_final)                                                # display the figure
