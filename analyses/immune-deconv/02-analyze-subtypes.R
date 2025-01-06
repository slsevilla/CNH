# Author: Komal S. Rathi updated, 2020-07 Kelsey Keith, 2025-01 S Chill updated
# script to perform immune characterization using R package immunedeconv
## and, when provided filtering params, will analyze immune cell distribution across subtypes

# load libraries
suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(immunedeconv)
  library(tidyverse)
  library(ggplot2)
})

# turn off warnings
options(warn = -1)

# parse parameters
option_list <- list(
  make_option(c("--expr_mat"), type = "character",
              help = "expression data: gene symbol x sample identifiers (.rds)"),
  make_option(c("--clin_file"), type = "character",
              help = "histologies file (.tsv)"),
  make_option(c("--deconv_method"), type = "character",
              help = "deconvolution method"),
  make_option(c("--molecular_subtype_list"), type = "character",
              help = "list of available molecular subtypes, hypen separated MB, WNT"),
  make_option(c("--output_dir"), type = "character", 
              help = "output directory")
)

opt <- parse_args(OptionParser(option_list = option_list))
expr_mat <- opt$expr_mat
clin_file <- opt$clin_file
deconv_method <- tolower(opt$deconv_method)
output_dir <- opt$output_dir

# Parse molecular_subtype_list
if (!is.null(opt$molecular_subtype_list)) {
  molecular_subtype_list <- strsplit(opt$molecular_subtype_list, "-")[[1]]
} else {
  molecular_subtype_list <- NULL
}

# output file
dir.create(output_dir, showWarnings = F, recursive = T)
output_file <- file.path(output_dir, paste0(deconv_method, "_output.rds"))

# method should be one of xcell or quantiseq
if (!(deconv_method %in% c("xcell", "quantiseq")))
  stop("Error: deconv_method must be one of xcell or quantiseq")

# read expression data 
expr_mat <- readRDS(expr_mat)

# read clinical data
# molecular_subtype_list=c("WT, WNT")
clin_file <- readr::read_tsv(clin_file, guess_max = 10000)
clin_file <- clin_file %>% 
  filter(Kids_First_Biospecimen_ID %in% colnames(expr_mat),
         experimental_strategy == "RNA-Seq",
         cohort == "GTEx" | !is.na(cancer_group)) # remove non-annotated samples


# filter by subtype, if there are enough samples associated
if (length(subset(clin_file,molecular_subtype %in% molecular_subtype_list))>0){
  print("Filtering by provided subtype...")
  clin_file <- clin_file %>% 
    filter(molecular_subtype %in% molecular_subtype_list) # add filter for subtype
}

# print if there are subtypes not in the cohort
if(length(setdiff(molecular_subtype_list, unique(clin_file$molecular_subtype)))>0){
  print("There are subtypes provided which are not in the dataset provided.")

  print("--These are included:")
  print(paste("--",sort(unique(clin_file$molecular_subtype))))

  print("--These are missing:")
  print(paste("--",setdiff(molecular_subtype_list, unique(clin_file$molecular_subtype))))
}

# process each group separately because xCell uses the variability among the samples for a linear transformation of the output score
# a group is a combination of cohort + cancer group or gtex group
n_groups <- clin_file %>%
  mutate(group = ifelse(cohort == "GTEx", gtex_group, cancer_group),
        group = paste0(cohort, "_", group)) %>%
        dplyr::select(Kids_First_Biospecimen_ID, group)

full_output <- plyr::dlply(.data = n_groups, .variables = "group", .fun = function(x) {
  # get kf ids for all samples in the group
  # print(x %>% pull(group) %>% unique)
  kf_ids <- x %>%
    pull(Kids_First_Biospecimen_ID)
  
  # at least two samples are needed for processing
  if(length(kf_ids) < 2){
    print("Only one sample available; skip processing")
    return()
  }
  
  # subset input matrix to group
  expr_mat_sub <- expr_mat %>%
    dplyr::select(kf_ids)
  
  # remove features with all 0 values
  expr_mat_sub <- expr_mat_sub[rowSums(expr_mat_sub) > 0,]
  
  # deconvolute using specified method
  print("Starting deconvolution...")
  deconv_output <- deconvolute(gene_expression = as.matrix(expr_mat_sub), 
                               method = deconv_method, arrays = F)
  
  # convert to long format
  deconv_output <- deconv_output %>%
    as.data.frame() %>%
    gather(Kids_First_Biospecimen_ID, fraction, -c(cell_type)) %>%
    as.data.frame()
})

# combine all groups
full_output <- dplyr::bind_rows(full_output)

# merge output with clinical data
full_output <- clin_file %>% 
  dplyr::select(Kids_First_Biospecimen_ID, cohort, sample_type, gtex_group, gtex_subgroup, cancer_group, molecular_subtype) %>%
  inner_join(full_output, by = "Kids_First_Biospecimen_ID") %>%
  mutate(method = deconv_method)

# save output to rds file
print("Writing deconv output to file...")
saveRDS(object = full_output, file = output_file)

# Analysis, if applicable
if (length(subset(clin_file,molecular_subtype %in% molecular_subtype_list))>0){
  print("Performing statistical analysis by subtypes..")
  stats_results <- full_output %>%
    group_by(cell_type) %>%
    summarise(p_value = kruskal.test(fraction ~ molecular_subtype)$p.value) %>%
    mutate(significant = ifelse(p_value < 0.05, "Yes", "No"))

  stats_output_file <- file.path(output_dir, "immune_cell_subtypes_stats.tsv")
  write_tsv(stats_results, stats_output_file)

  # get category subtype
  category_subtype=unique(clin_file$pathology_diagnosis)

  # Create visuals
  print("Generating visualization")
  immune_cell_plot <- full_output %>%
    ggplot(aes(x = molecular_subtype, y = fraction, fill = molecular_subtype)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 0.5, alpha = 0.7) +
    facet_wrap(~ cell_type, scales = "free_y", ncol = 4) +  # Adjust columns for clarity
    theme_bw() +
    labs(
      title = paste0("Immune Cell Fractions Across ", category_subtype, " Subtypes"),
      x = paste0(category_subtype, " Subtype"),
      y = "Immune Cell Fraction"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8), # Adjust text size and alignment
      axis.text.y = element_text(size = 8),
      strip.text = element_text(size = 10, face = "bold"),  # Facet label text
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), # Center title
      plot.margin = margin(1, 1, 2, 1, "cm")  # Increase margins to avoid clipping
    ) +
    coord_cartesian(clip = "off") # Ensure nothing is clipped
  
  # Save Plot
  plot_output_file <- file.path(output_dir, "immune_cell_subtypes_plot.png")
  ggsave(
    filename = plot_output_file,
    plot = immune_cell_plot,
    width = 16, height = 12, dpi = 300
  )
  # Complete
  print(">>> Analysis of immune cell distribution completed successfully!")
}
