library(tidyverse)
library(FSA)
library(patchwork)
library(scales)
library(knitr)
library(ggbeeswarm)
library(car)
library(FactoMineR)
library(factoextra)
library(vegan)
library(DiagrammeR)

# ======================================================================
# OBJECTIVE 1
# ======================================================================

# Clean Physicochemical Data with global outlier protection (>40°C)
WQP_clean <- read.csv("WQP_R.csv", header = TRUE, sep = ";") %>% 
  distinct() %>% 
  mutate(across(c(dissolved_oxgen_ppm, dissolved_oxygen_percent, ec, ph, temperature), ~ as.numeric(gsub(",", ".", .))), 
         location = str_trim(location)) %>% 
  # GLOBAL QUALITY FILTER: Flag and remove all 7 equipment errors > 40°C
  mutate(temperature = if_else(temperature > 40, as.numeric(NA), temperature))

# FIXED CLARITY DATA PIPE: Removed distinct() to keep all your 95cm observations!
CT_clean <- CT_final %>% 
  mutate(
    location = str_trim(location), 
    clarity_final = as.numeric(gsub(",", ".", clarity_final))
  )

# ======================================================================
# 2. UNIFIED 5-SITE STRUCTURE FUNCTION
# ======================================================================
prepare_plot_data <- function(df) {
  mining <- df %>% filter(location %in% c("Blyde", "Moremela", "Treur")) %>% 
    mutate(
      panel_group = "Illegal mining sites",
      site_label = case_when(
        location == "Blyde" ~ "Blyde",
        location == "Moremela" ~ "Moremela",
        location == "Treur" ~ "Treur (Ref)"
      ),
      disturbance_class = case_when(
        location == "Blyde" ~ "Illegal mining",
        location == "Moremela" ~ "Illegal mining",
        location == "Treur" ~ "Illegal mining reference"
      )
    )
  
  ag <- df %>% filter(location %in% c("Marakalala", "Mohlatsengwana")) %>% 
    mutate(
      panel_group = "Intensive agriculture sites",
      site_label = case_when(
        location == "Marakalala" ~ "Marakalala",
        location == "Mohlatsengwana" ~ "Mohlatsengwana (Ref)"
      ),
      disturbance_class = case_when(
        location == "Marakalala" ~ "Intensive agriculture",
        location == "Mohlatsengwana" ~ "Intensive agriculture reference"
      )
    )
  
  bind_rows(mining, ag) %>% 
    mutate(
      panel_group = factor(panel_group, levels = c("Illegal mining sites", "Intensive agriculture sites")),
      site_label = factor(site_label, levels = c("Blyde", "Moremela", "Treur (Ref)", "Marakalala", "Mohlatsengwana (Ref)")),
      disturbance_class = factor(disturbance_class, levels = c(
        "Illegal mining", "Illegal mining reference", "Intensive agriculture", "Intensive agriculture reference"
      ))
    )
}

# Generate the processed plotting datasets from the cleaned baselines
clarity_plot_df <- prepare_plot_data(CT_clean)
WQP_plot_df <- prepare_plot_data(WQP_clean)

# ======================================================================
# 3. COLOR SYSTEM & GRID PLOT GENERATOR
# ======================================================================
custom_colors <- c(
  "Illegal mining" = "#FFC697",
  "Illegal mining reference" = "#63e5ff",
  "Intensive agriculture" = "#CBD3AD",
  "Intensive agriculture reference" = "#63e5ff"
)

plot_panel <- function(data, y_var, y_label) {
  ggplot(data, aes(x = site_label, y = {{y_var}}, fill = disturbance_class)) + 
    geom_boxplot(alpha = 0.8, outlier.shape = NA, linewidth = 0.5) + 
    geom_jitter(width = 0.15, size = 1.3, alpha = 0.4, color = "black") + 
    facet_wrap(~ panel_group, scales = "free_x") + 
    scale_fill_manual(values = custom_colors) + 
    labs(x = NULL, y = y_label) + 
    theme_bw(base_size = 11) + 
    theme(
      axis.title.y = element_text(face = "bold", size = 11, margin = margin(r = 8)), 
      axis.text.x = element_text(size = 9, color = "black", angle = 25, hjust = 1), 
      axis.text.y = element_text(size = 9, color = "black"), 
      strip.text = element_text(face = "bold", size = 10), 
      strip.background = element_rect(fill = "gray96"), 
      panel.grid.minor = element_blank(), 
      panel.spacing = unit(0.2, "inch"), 
      legend.position = "none"
    )
}

# ======================================================================
# 4. BUILD THE INDIVIDUAL PANELS
# ======================================================================
p_a <- plot_panel(clarity_plot_df, clarity_final, "Clarity (cm)")
p_b <- plot_panel(WQP_plot_df, ec, "EC (uS/cm)")
p_c <- plot_panel(WQP_plot_df, ph, "pH level")
p_d <- plot_panel(WQP_plot_df, temperature, "Temperature (°C)")
p_e <- plot_panel(WQP_plot_df, dissolved_oxgen_ppm, "Dissolved Oxygen (ppm)")

# ======================================================================
# 5. ASSEMBLE WITH PATCHWORK & EXPORT (Figure 4.1)
# ======================================================================
layout_design <- "
AB
CD
EE
"

figure_4_1 <- p_a + p_b + p_c + p_d + p_e + 
  patchwork::plot_layout(design = layout_design) + 
  patchwork::plot_annotation(
    tag_levels = 'a', 
    tag_prefix = '(', 
    tag_suffix = ')', 
    caption = "River site"
  ) & theme(
    plot.tag = element_text(face = "bold", size = 14), 
    plot.caption = element_text(face = "bold", size = 12, hjust = 0.5, margin = margin(t = 10))
  )

ggsave("Figure_4_1_fixed.png", plot = figure_4_1, width = 10, height = 13, dpi = 300)

# ======================================================================
# DUNN'S TEST
# ======================================================================

# 1. Run post-hoc testing for Clarity
cat("--- Pairwise Dunn's Test: Water Clarity ---\n")
dunn_clarity <- dunnTest(clarity_final ~ site_label, data = clarity_plot_df, method = "bh")
print(dunn_clarity$res %>% filter(Comparison %in% c("Blyde - Treur (Ref)", "Moremela - Treur (Ref)", "Marakalala - Mohlatsengwana (Ref)")))

# 2. Run post-hoc testing for Physicochemical Variables (Looping dynamically)
wqp_params <- c("ec", "ph", "temperature", "dissolved_oxgen_ppm")
for (param in wqp_params) {
  cat(paste0("\n--- Pairwise Dunn's Test: ", toupper(param), " ---\n"))
  dunn_wqp <- dunnTest(formula(paste(param, "~ site_label")), data = WQP_plot_df, method = "bh")
  filtered_res <- dunn_wqp$res %>% filter(Comparison %in% c("Blyde - Treur (Ref)", "Moremela - Treur (Ref)", "Marakalala - Mohlatsengwana (Ref)"))
  print(filtered_res)
}

# ======================================================================
# OBJECTIVE 2
# ======================================================================
MS <- read.csv("MS_R.csv", header = TRUE, sep = ";")

MS_final <- MS %>%
  select(location, observation_date, category, species, average_score) %>%
  mutate(
    average_score = as.numeric(str_replace(as.character(average_score), ",", ".")),
    species = as.integer(species),
    location = str_trim(as.character(location)),
    category = str_trim(as.character(category))
  ) %>%
  filter(!is.na(average_score), !is.na(location)) %>%
  mutate(
    Disturbance = case_when(
      location == "Blyde" ~ "Illegal mining",
      location == "Moremela" ~ "Illegal mining",
      location == "Treur" ~ "Illegal mining reference",
      location == "Marakalala" ~ "Intensive agriculture",
      location == "Mohlatsengwana" ~ "Intensive agriculture reference",
      TRUE ~ "Unknown"
    ),
    Catchment = case_when(
      location %in% c("Blyde", "Moremela", "Treur") ~ "Illegal mining sites",
      location %in% c("Marakalala", "Mohlatsengwana") ~ "Intensive agriculture sites",
      TRUE ~ "Unknown"
    ),
    site_label = case_when(
      location == "Blyde" ~ "Blyde",
      location == "Moremela" ~ "Moremela",
      location == "Treur" ~ "Treur (Ref)",
      location == "Marakalala" ~ "Marakalala",
      location == "Mohlatsengwana" ~ "Mohlatsengwana (Ref)"
    )
  ) %>%
  mutate(
    panel_group = factor(Catchment, levels = c("Illegal mining sites", "Intensive agriculture sites")),
    site_label = factor(site_label, levels = c("Blyde", "Moremela", "Treur (Ref)", "Marakalala", "Mohlatsengwana (Ref)")),
    Disturbance = factor(Disturbance, levels = c("Illegal mining", "Illegal mining reference", 
                                                 "Intensive agriculture", "Intensive agriculture reference"))
  )

# --- 1. Illegal Mining Catchment Framework Test ---
mining_sub <- MS_final %>% 
  filter(Catchment == "Illegal mining sites")

# --- 1. ILLEGAL MINING CATCHMENT: IMPACT VS. REFERENCE ---
mining_ttest <- t.test(average_score ~ Disturbance, data = mining_sub, var.equal = FALSE)
print(mining_ttest)

# --- 2. Intensive Agriculture Catchment Framework Test ---
ag_sub <- MS_final %>% 
  filter(Catchment == "Intensive agriculture sites")

# --- 2. INTENSIVE AGRICULTURE CATCHMENT: IMPACT VS. REFERENCE ---\n")
ag_ttest <- t.test(average_score ~ Disturbance, data = ag_sub, var.equal = FALSE)
print(ag_ttest)

# --- Continuous Health Metrics Plot (Figure 4.2) ---
custom_thesis_colors <- c(
  "Illegal mining" = "#FFC697", 
  "Intensive agriculture" = "#CBD3AD", 
  "Illegal mining reference" = "#63e5ff", 
  "Intensive agriculture reference" = "#63e5ff" 
)

ggplot(MS_final, aes(x = site_label, y = average_score, fill = Disturbance)) +
  geom_hline(yintercept = 5.3, linetype = "dashed", color = "#8E44AD", alpha = 0.4, linewidth = 0.5) + 
  geom_hline(yintercept = 5.7, linetype = "dashed", color = "#E67E22", alpha = 0.4, linewidth = 0.5) + 
  geom_hline(yintercept = 6.2, linetype = "dashed", color = "#D4AF37", alpha = 0.4, linewidth = 0.5) + 
  geom_hline(yintercept = 7.2, linetype = "dashed", color = "#27AE60", alpha = 0.4, linewidth = 0.5) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.6, linewidth = 0.6) +
  geom_quasirandom(width = 0.15, size = 1.8, alpha = 0.4, color = "black") +
  facet_wrap(~ panel_group, scales = "free_x") +
  scale_fill_manual(values = custom_thesis_colors) +
  labs(x = "River site", y = "miniSASS biological health score") +
  theme_bw(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "gray96"),
    axis.title.x = element_text(face = "bold", size = 13, margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", size = 13, margin = margin(r = 10)),
    axis.text.x = element_text(size = 11, color = "black", angle = 25, hjust = 1),
    axis.text.y = element_text(size = 11, color = "black"),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )
ggsave("Figure_4_2_continuous.png", width = 10, height = 6, dpi = 300)

##### Categorical Stacked Bar Plot (Figure 4.3) ######
  
  ecological_order <- c("1. Natural", "2. Good", "3. Fair", "4. Poor", "5. Very poor")
MS_categorical_clean <- MS_final %>%
  mutate(category = factor(trimws(category), levels = ecological_order)) %>%
  filter(!is.na(category)) %>%
  count(Catchment, site_label, category) %>%
  group_by(site_label) %>%
  mutate(percentage = n / sum(n)) %>%
  ungroup() %>%
  arrange(site_label, category)
ggplot(MS_categorical_clean, aes(x = site_label, y = percentage, fill = category)) +
  geom_col(color = "black", width = 0.6, linewidth = 0.4) +
  geom_text(
    aes(label = ifelse(percentage >= 0.03, scales::percent(percentage, accuracy = 1), ""), group = category),
    position = position_stack(vjust = 0.5),
    size = 3.5,
    fontface = "bold",
    color = "black"
  ) +
  facet_wrap(~ Catchment, scales = "free_x") +
  scale_fill_manual(
    values = c("1. Natural" = "#1B5E20", "2. Good" = "#388E3C", "3. Fair" = "#FBC02D", "4. Poor" = "#F57C00", "5. Very poor" = "#D32F2F"),
    guide = guide_legend(reverse = FALSE)
  ) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0, 0.05, 0)) +
  labs(
    title = "Figure 4.3: Proportional distribution of MiniSASS ecological-condition categories across river monitoring sites",
    x = "River site", y = "Percentage of samples (%)", fill = "Ecological condition"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5, margin = margin(b = 12)),
    strip.text = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "gray96"),
    axis.title.x = element_text(face = "bold", size = 13, margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", size = 13, margin = margin(r = 10)),
    axis.text.x = element_text(size = 11, color = "black", angle = 25, hjust = 1),
    axis.text.y = element_text(size = 11, color = "black"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )
ggsave("Figure_4_3_categorical_final.png", width = 11, height = 6, dpi = 300)

# ======================================================================
#  OBJECTIVE 3
# ======================================================================

  WQP_plot_df_clean <- WQP_plot_df %>%
  mutate(temperature = if_else(temperature > 40, as.numeric(NA), temperature))
wqp_aggregated_clean <- WQP_plot_df_clean %>%
  group_by(location, observation_date) %>%
  summarise(
    ec = mean(ec, na.rm = TRUE),
    ph = mean(ph, na.rm = TRUE),
    temperature = mean(temperature, na.rm = TRUE),
    dissolved_oxgen_ppm = mean(dissolved_oxgen_ppm, na.rm = TRUE),
    .groups = "drop"
  )
regression_data_clean <- ms_aggregated %>%
  inner_join(wqp_aggregated_clean, by = c("location", "observation_date")) %>%
  inner_join(ct_aggregated, by = "location") %>%
  filter(!is.na(average_score)) %>%
  mutate(
    Disturbance = case_when(
      location == "Blyde" ~ "Illegal mining",
      location == "Moremela" ~ "Illegal mining",
      location == "Treur" ~ "Illegal mining reference",
      location == "Marakalala" ~ "Intensive agriculture",
      location == "Mohlatsengwana" ~ "Intensive agriculture reference",
      TRUE ~ "Unknown"
    )
  )

##### Ordinary Least-Squares Multiple Linear Regression #####
  
  ols_model_final <- lm(average_score ~ clarity_final + ec + ph + temperature + dissolved_oxgen_ppm, data = regression_data_clean)

print(summary(ols_model_final))
print(vif(ols_model_final))

# --- Bivariate Spearman Rank Correlations ---

predictors <- c("clarity_final", "ec", "ph", "temperature", "dissolved_oxgen_ppm")
predictor_labels <- c("Water clarity", "EC", "pH", "Temperature", "DO")
for (i in 1:length(predictors)) {
  var <- predictors[i]
  label <- predictor_labels[i]
  clean_df <- regression_data_clean %>% filter(!is.na(average_score))
  clean_pairs <- clean_df[!is.na(clean_df[[var]]), ]
  cat(paste0("\n--- Spearman Correlation: miniSASS Score vs. ", label, " ---\n"))
  test_res <- cor.test(clean_pairs$average_score, clean_pairs[[var]], method = "spearman", exact = FALSE)
  print(test_res)
}

# ======================================================================
# OBJECTIVE 4: UNIFIED MULTIVARIATE PCA DATA PIPELINE (THESIS READY)
# ======================================================================

pca_matrix_clean <- regression_data_clean %>%
  select(clarity_final, ec, ph, temperature, dissolved_oxgen_ppm, average_score) %>%
  rename(
    `Water Clarity`    = clarity_final,
    `EC`               = ec,
    `pH`               = ph,
    `Temperature`      = temperature,
    `DO` = dissolved_oxgen_ppm,
    `miniSASS Score`   = average_score
  )

pca_engine_clean <- PCA(pca_matrix_clean, scale.unit = TRUE, graph = FALSE)

grouping_factor <- regression_data_clean %>%
  mutate(Broad_Group = case_when(
    location == "Treur"          ~ "Illegal mining reference (Treur)",
    location == "Mohlatsengwana" ~ "Intensive agriculture reference (Mohlatsengwana)",
    location %in% c("Blyde", "Moremela") ~ "Mining-impact grouping (Blyde & Moremela)",
    location == "Marakalala"     ~ "Agricultural-impact grouping (Marakalala)",
    TRUE                         ~ "Unknown"
  )) %>%
  pull(Broad_Group)

custom_group_colors <- c(
  "Illegal mining reference (Treur)"                  = "#2980B9", # Dark Blue
  "Intensive agriculture reference (Mohlatsengwana)" = "#1ABC9C", # Teal Green
  "Mining-impact grouping (Blyde & Moremela)"        = "#E67E22", # Orange
  "Agricultural-impact grouping (Marakalala)"        = "#27AE60"  # Green
)

pca_biplot_final <- fviz_pca_biplot(
  pca_engine_clean,
  repel = TRUE,                        # Keeps labels separated to prevent overlaps
  geom.ind = "point",                  
  pointsize = 2.5,
  habillage = factor(grouping_factor), # Uses the four separate categories
  col.var = "black",                   
  arrowsize = 0.6,
  labelsize = 3.5,                     # Smaller variable labels to prevent crowding
  addEllipses = TRUE,                  
  ellipse.level = 0.95,
  ellipse.type = "convex",
  palette = custom_group_colors    
) +
  scale_fill_manual(values = alpha(custom_group_colors, 0.03)) +
  labs(
    title = "Figure 4.5: PCA ordination of physicochemical conditions and MiniSASS across mining-impacted, agriculture-impacted and reference river sites",
    x = "Principal Component 1 (26.6%)",  # Updated to match console (26.64%)
    y = "Principal Component 2 (23.2%)"   # Updated to match console (23.23%)
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 10, hjust = 0.5, margin = margin(b = 15)),
    axis.title.x = element_text(face = "bold", size = 11, margin = margin(t = 8)),
    axis.title.y = element_text(face = "bold", size = 11, margin = margin(r = 8)),
    legend.position = "right",
    legend.title = element_blank(),       
    legend.text = element_text(size = 9),
    panel.grid.minor = element_blank()
  )

if(!dir.exists("plots")) dir.create("plots")
ggsave("plots/Figure_4_5_PCA_Final.png", plot = pca_biplot_final, width = 11, height = 7, dpi = 300)

# --- Principal component percentages
pca_engine_clean$eig[1:2, 2]

# ======================================================================
# 6. PERMADISP AND PERMANOVA
# ======================================================================
  
  pca_matrix_testing <- regression_data_clean %>%
  select(location, clarity_final, ec, ph, temperature, dissolved_oxgen_ppm, average_score) %>%
  na.omit()
numeric_matrix <- pca_matrix_testing %>%
  select(clarity_final, ec, ph, temperature, dissolved_oxgen_ppm, average_score) %>%
  scale()
grouping_vector <- case_when(
  pca_matrix_testing$location == "Treur" ~ "Illegal mining reference (Treur)",
  pca_matrix_testing$location == "Mohlatsengwana" ~ "Intensive agriculture reference (Mohlatsengwana)",
  pca_matrix_testing$location %in% c("Blyde", "Moremela") ~ "Mining-impact grouping (Blyde & Moremela)",
  pca_matrix_testing$location == "Marakalala" ~ "Agricultural-impact grouping (Marakalala)",
  TRUE ~ "Unknown"
)
grouping_factor <- factor(grouping_vector)
valid_idx <- grouping_factor != "Unknown"
numeric_matrix <- numeric_matrix[valid_idx, ]
grouping_factor <- factor(grouping_factor[valid_idx])
dist_matrix <- vegdist(numeric_matrix, method = "euclidean")

##### STEP A: PERMADISP
  
perm_disp_model <- betadisper(dist_matrix, group = grouping_factor, type = "centroid")
perm_disp_results <- permutest(perm_disp_model, permutations = 999)
print(perm_disp_results)

##### STEP B: PERMANOVA
  
permanova_model <- adonis2(numeric_matrix ~ grouping_factor, method = "euclidean", permutations = 999)
print(permanova_model)

# ======================================================================
# Figure 4.6: CONCEPTUAL SYNTHESIS DIAGRAM
# ======================================================================
grViz("
  digraph conceptual_synthesis {
    
    # Graph setup for a clean, wide layout
    graph [layout = dot, rankdir = LR, splines = ortho, nodesep = 0.4, ranksep = 0.5]
    
    # Global node styling by structural layer
    node [fontname = 'Helvetica', style = 'filled, rounded', penwidth = 2, fixedsize = false, margin = '0.2,0.1']
    
    # --- NODE DEFINITIONS BY LAYER ---
    
    # Layer 1: Anthropogenic Drivers (Grey/Sturdy)
    node [fillcolor = '#F2F4F4', color = '#7F8C8D', fontcolor = '#2C3E50', shape = box]
    D1 [label = 'Mining Disturbance\n(Landscape Level)']
    D2 [label = 'Agricultural Disturbance\n(Landscape Level)']
    
    # Layer 2: Physicochemical Responses (Soft Blue)
    node [fillcolor = '#EBF5FB', color = '#2E86C1', fontcolor = '#1B4F72', shape = box]
    P1 [label = 'Physicochemical Profile:\n• Lower Clarity\n• Higher EC\n(vs. Reference)']
    P2 [label = 'Physicochemical Profile:\n• Altered pH\n• Higher Temperature\n(vs. Reference)']
    
    # Layer 3: Ecological Responses (Soft Green)
    node [fillcolor = '#E8F8F5', color = '#117A65', fontcolor = '#0E6251', shape = box]
    E1 [label = 'Ecological Profile:\n• Lower miniSASS Biotic Condition\n(vs. Mining Reference)']
    E2 [label = 'Ecological Profile:\n• Heterogeneous miniSASS Response\n• No statistically supported continuous diff.']

    # --- EDGE / CONNECTOR DEFINITIONS ---
    # Using clean, thin lines with a descriptive label to explicitly state 'association' instead of 'causation'
    edge [color = '#7F8C8D', arrowsize = 0.8, penwidth = 1.5, fontname = 'Helvetica-Oblique', fontsize = 10, fontcolor = '#7F8C8D']
    
    # Mining Pathway Connections
    D1 -> P1 [label = ' associates with ']
    P1 -> E1 [label = ' tracks with ']
    
    # Agricultural Pathway Connections
    D2 -> P2 [label = ' associates with ']
    P2 -> E2 [label = ' tracks with ']
  }
")

