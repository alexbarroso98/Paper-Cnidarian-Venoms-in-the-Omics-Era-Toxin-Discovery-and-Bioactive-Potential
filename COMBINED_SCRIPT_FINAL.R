# Carregar bibliotecas
library(readxl)
library(dplyr)
library(writexl)
library(ggplot2)
library(tidyr)
library(patchwork)
library(cowplot)

# Organizar dados

df <- read_excel("Tabela_Final_Contagens_Omicas.xlsx")

toxin_cols <- names(df)[
  which(names(df) == "PLA2"):
    (which(names(df) == "T") - 1)
]

family_order <- c(
  
  #------------------------
  # Toxic enzymes
  #------------------------
  "Metalloproteinases",
  "PLA2",
  "PLA1",
  "PLAB",
  "PLAD",
  "Proteases",
  "Hyaluronidase",
  "5'nucleotidase",
  "Chitinase",
  "Other Enzymes",
  
  #------------------------
  # Protease inhibitors
  #------------------------
  "KTx II/Kunitz",
  "Kunitz-like",
  "Cystatin",
  "Others Prot.Inh.",
  
  #------------------------
  # Neurotoxins
  #------------------------
  "NaTx I",
  "NaTx II",
  "NaTx III",
  "KTx I",
  "KTx III/BDS",
  "KTx IV/BBH/St.Cl.9",
  "KTx V",
  "KTx VI",
  "SA8",
  "ShK-like",
  "SCRiP",
  "Acrorhagin",
  "Conotoxins",
  "Latrotoxin",
  "Turripeptides",
  "Three-Finger",
  "Other Neurotoxins",
  
  #------------------------
  # Pore-forming toxins
  #------------------------
  "Actinoporin",
  "JFT",
  "MACPF",
  "Hydralysin",
  "Other PFTs",
  
  #------------------------
  # Venom-associated proteins
  #------------------------
  "Lectin/Snaclec",
  "CRiSP/CAP",
  "EGF-like",
  "PHAB",
  "Venom allergen",
  "Venom coagulation factor",
  "Venom prothrombin activator",
  
  #------------------------
  # Other / Unknown
  #------------------------
  "Others",
  "Unknown"
)

# Contar espécies por família de toxinas e tipo de ómica

summary_table <- data.frame(
  Family = toxin_cols,
  Proteomics = NA,
  Transcriptomics = NA
)

for(i in seq_along(toxin_cols)){
  
  tox <- toxin_cols[i]
  
  summary_table$Proteomics[i] <-
    df %>%
    filter(Omics == "Proteomics") %>%
    summarise(n = sum(.data[[tox]] == 1, na.rm = TRUE)) %>%
    pull(n)
  
  summary_table$Transcriptomics[i] <-
    df %>%
    filter(Omics == "Transcriptomics") %>%
    summarise(n = sum(.data[[tox]] == 1, na.rm = TRUE)) %>%
    pull(n)
  
}

#Remover as familias sem ocurrencias
summary_table <- summary_table %>%
  filter(Proteomics + Transcriptomics > 0)

#GRAFICO 1
# Transformar para formato longo
plot_data <- summary_table %>%
  pivot_longer(
    cols = c(Proteomics, Transcriptomics),
    names_to = "Omics",
    values_to = "Species"
  )

# IMPORTANT: use the same order as Figure B
plot_data$Family <- factor(
  plot_data$Family,
  levels = family_order
)

# Gráfico
p <- ggplot(plot_data,
            aes(x = Species,
                y = Family,
                fill = Omics)) +
  
  annotate("rect", xmin = -3, xmax = -2, ymin = 0.5, ymax = 10.5, fill = "#F4CCCC") +
  annotate("rect", xmin = -3, xmax = -2, ymin = 10.5, ymax = 14.5, fill = "#FFF2CC") +
  annotate("rect", xmin = -3, xmax = -2, ymin = 14.5, ymax = 31.5, fill = "#D9EAD3") +
  annotate("rect", xmin = -3, xmax = -2, ymin = 31.5, ymax = 36.5, fill = "#D0E0E3") +
  annotate("rect", xmin = -3, xmax = -2, ymin = 36.5, ymax = 43.5, fill = "#D9D2E9") +
  
  geom_text(aes(label = Species),
            position = position_dodge(width = 0.75),
            hjust = -0.2,
            size = 4) +
  
  geom_col(position = position_dodge(width = 0.75),
           width = 0.65) +
  
  scale_fill_manual(values = c(
    Proteomics = "#2C6DB2",
    Transcriptomics = "#4DB6AC"
  )) +
  
  labs(
    x = "Number of species",
    y = "Toxin family"
  ) +
  
  coord_cartesian(xlim = c(-3, NA), clip = "off") +
  
  theme_classic(base_size = 14) +
  
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(
      size = 13,
      colour = "black"
    ),
    axis.title.x = element_text(size = 17, face = "bold"),
    axis.title.y = element_text(size = 19, face = "bold"),
    panel.grid.major.y = element_line(
      colour = "grey85",
      linewidth = 0.3
    ),
    panel.grid.minor = element_blank()
  )

#------------------------------------------------------------
# FIGURE B DATA
#------------------------------------------------------------

bubble_data <- df %>%
  pivot_longer(
    cols = all_of(toxin_cols),
    names_to = "Family",
    values_to = "Presence"
  ) %>%
  filter(Presence == 1) %>%
  group_by(Omics, Order, Family) %>%
  summarise(
    N_species = n_distinct(Species),
    .groups = "drop"
  )

#------------------------------------------------------------
# Number of species per order
#------------------------------------------------------------

order_counts <- df %>%
  group_by(Omics, Order) %>%
  summarise(
    N_species_total = n_distinct(Species),
    .groups = "drop"
  )

bubble_data <- bubble_data %>%
  left_join(
    order_counts,
    by = c("Omics", "Order")
  ) %>%
  mutate(
    Order_label = paste0(
      Order, " (", N_species_total, ")"
    )
  )

#------------------------------------------------------------
# Remove empty combinations
#------------------------------------------------------------

bubble_data <- bubble_data %>%
  group_by(Omics, Order_label) %>%
  filter(sum(N_species) > 0) %>%
  ungroup()

#------------------------------------------------------------
# Apply same toxin-family order as Figure A
#------------------------------------------------------------

bubble_data$Family <- factor(
  bubble_data$Family,
  levels = family_order
)

#------------------------------------------------------------
# Order cnidarian orders
#------------------------------------------------------------

proteomics_order <- c(
  "Actiniaria",
  "Corallimorpharia",
  "Scleractinia",
  "Zoantharia",
  "Octocorallia",
  "Cubozoa",
  "Hydrozoa",
  "Scyphozoa",
  "Staurozoa",
  "Myxozoa",
  "Polypodiozoa"
)

transcriptomics_order <- c(
  "Actiniaria",
  "Antipatharia",
  "Ceriantharia",
  "Corallimorpharia",
  "Scleractinia",
  "Zoantharia",
  "Cubozoa",
  "Hydrozoa",
  "Scyphozoa",
  "Staurozoa",
  "Myxozoa",
  "Polypodiozoa"
)

order_levels <- c(
  
  order_counts %>%
    filter(Omics == "Proteomics") %>%
    mutate(
      Order = factor(
        Order,
        levels = proteomics_order
      )
    ) %>%
    arrange(Order) %>%
    mutate(
      Label = paste0(
        Order,
        " (",
        N_species_total,
        ")"
      )
    ) %>%
    pull(Label),
  
  order_counts %>%
    filter(Omics == "Transcriptomics") %>%
    mutate(
      Order = factor(
        Order,
        levels = transcriptomics_order
      )
    ) %>%
    arrange(Order) %>%
    mutate(
      Label = paste0(
        Order,
        " (",
        N_species_total,
        ")"
      )
    ) %>%
    pull(Label)
)

bubble_data$Order_label <- factor(
  bubble_data$Order_label,
  levels = unique(order_levels)
)

#------------------------------------------------------------
# FIGURE B
#------------------------------------------------------------

p2 <- ggplot(
  bubble_data,
  aes(
    x = Order_label,
    y = Family,
    size = N_species,
    colour = Omics
  )
) +
  
  geom_point(alpha = 0.9) +
  
  facet_grid(
    . ~ Omics,
    scales = "free_x",
    space = "free_x"
  ) +
  
  scale_colour_manual(
    values = c(
      Proteomics = "#2C6DB2",
      Transcriptomics = "#4DB6AC"
    )
  ) +
  
  scale_size(
    range = c(2.5, 11),
    name = "Number of Species"
  ) +
  
  labs(
    x = "Cnidarian Group",
    y = NULL
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    
    legend.position = "right",
    
    legend.title = element_text(
      size = 14,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 13
    ),
    
    legend.key.size = unit(0.8, "cm"),
    
    strip.background = element_blank(),
    
    strip.text = element_text(
      face = "bold",
      size = 18
    ),
    
    axis.text.x = element_text(
      size = 16,
      angle = 45,
      hjust = 1,
      vjust = 1,
      colour = "black"
    ),
    
    axis.title.x = element_text(
      size = 18,
      face = "bold"
    ),
    
    # REMOVE y-axis from B
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    
    panel.grid.major.y = element_line(
      colour = "grey85",
      linewidth = 0.3
    ),
    
    panel.grid.minor = element_blank()
  )

p2

#------------------------------------------------------------
# CATEGORY LEGEND
#------------------------------------------------------------

category_legend_data <- data.frame(
  Category = c(
    "Enzymes",
    "Protease Inhibitors",
    "Neurotoxins",
    "Pore-Forming Toxins",
    "Mixed Functions"
  ),
  Color = c(
    "#F4CCCC",
    "#FFF2CC",
    "#D9EAD3",
    "#D0E0E3",
    "#D9D2E9"
  )
)

category_legend <- ggplot(
  category_legend_data,
  aes(x = 1, y = rev(seq_along(Category)))
) +
  geom_tile(
    aes(fill = Color),
    width = 0.20,
    height = 0.65
  ) +
  geom_text(
    aes(label = Category),
    hjust = 0,
    nudge_x = 0.15,
    size = 3.5
  ) +
  scale_fill_identity() +
  xlim(0.85, 3.5) +
  ylim(0.3, 5.7) +
  labs(title = "Category") +
  theme_void() +
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0
    ),
    plot.margin = margin(0, 0, 0, 0)
  )

#------------------------------------------------------------
# EXTRACT EXISTING LEGEND
#------------------------------------------------------------

p2_nolegend <- p2 +
  theme(legend.position = "none")

p2_legend <- cowplot::get_legend(
  p2 +
    theme(
      legend.position = "right"
    )
)

right_panel <- cowplot::plot_grid(
  ggdraw(p2_legend),
  category_legend,
  ncol = 1,
  rel_heights = c(1.4, 1)
)

#------------------------------------------------------------
# COMBINE FIGURE A + FIGURE B
#------------------------------------------------------------

combined <- p + p2_nolegend + right_panel +
  plot_layout(
    widths = c(1.05, 3, 0.9)
  )

combined

ggsave(
  "Figure_AB_combined.svg",
  plot = combined,
  width = 20,
  height = 15,
  bg = "white"
)

ggsave(
  "Figure_AB_combined.png",
  plot = combined,
  width = 20,
  height = 15,
  dpi = 600,
  bg = "white"
)