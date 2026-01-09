# --------------------------------------------------------------------------------------------------------
# Script: 01_peptide_preprocess.R
# Description: Metaproteomics data preprocessing at peptide level
# Author: Andrés Arroyo Barea
# Date: 2026-01-08
# --------------------------------------------------------------------------------------------------------
source("code/00_setup_environment.R")
source("code/01_config.R")

# -----------------------------
# Objetive
# -----------------------------


# -----------------------------
# Load Peptide Data
# -----------------------------
peptides <- read_metap_data("raw", current_run, "peptide")

# TODO: This should be an external step where names should be properly stablished because this patters could change
# between runs.
# Ex. Samples should follow the same patter from begging: ID_203
colnames(peptides) <- janitor::make_clean_names(colnames(peptides))
colnames(peptides) <- gsub("^x", "ID_", colnames(peptides))

# Print a report to identify var types properly.
glimpse(peptides)

# Missing report

# Sample-metrics filtering according to sample filtering in metadata.
peptides_processed <- filter_samples(peptides, samples_to_include = samples, verbose = T)

# Generación de nuevas variables de interes. Identificación de peptidos que mapean proteínas humanas y de peptidos que mapean a una o varias proteínas. 
peptides_processed <- peptides_processed %>%
  
  mutate(
    
    # Identify human proteins
    human_protein = is_human_protein(protein, "HUMAN"),
    
    # Identify unique peptides
    unique_peptide = is_unique_peptide(mapped_proteins),
    
    # Sum intensity per each peptide. 
    intensity_sum = rowSums(across(matches("^ID_\\d+_intensity")), na.rm = T), 
    
    # Intensity equal to 0.
    intensity_sum_0 = case_when(intensity_sum == 0 ~ "YES", TRUE ~ "NO") ,
    
    
    # CONTINUE HERE.
    
    # Precomputar la información para el 50% de las muestras NDMM .
    cond1_sum = rowSums(across(all_of(cond_list[[1]])) != 0),
    
    # Precomputar la información para el 50% de las muestras RMM.
    cond2_sum = rowSums(across(all_of(cond_list[[2]])) != 0),
    
    # Capturar los peptidos que están en el 50% de los pacientes NDMM y en el 50% de los pacientes RRMM.
    peptides_50_per_groups = case_when(
      cond1_sum >= ceiling(length(cond_list[[1]]) * 0.5) & cond2_sum >= length(cond_list[[2]]) * 0.5 ~ "YES", 
      TRUE ~ "NO"),
    
    # El filtro final para capturar los peptidos relevantes excluye los peptidos humanos, los péptidos
    # que no son únicos, los péptidos con intensidad 0 en todas las muestras y los peptidos que no estén
    # presentes al menos en el 50% de las muestras de cada grupo..
    # SE USARÁ EN ANÁLISIS DE PEPTIDOS.
    filtered_peptides = factor(case_when(
      human_protein == "NO" & unique_peptide == "YES" & intensity_sum > 0 & peptides_50_per_groups == "YES" ~ "YES",
      TRUE ~ "NO")),
    
    # Capturar los peptidos todo o nada en los pacientes NDMM.
    all_nothing_pep_cond1 = factor(case_when(
      human_protein == "NO" & unique_peptide == "YES" & intensity_sum > 0 & peptides_50_per_groups == "YES" & cond1_sum >= ceiling(length(cond_list[[1]])*0.5) & cond2_sum == 0 ~ "YES", 
      TRUE ~ "NO")),
    
    # Capturar los peptidos todo o nada en los pacientes RRMM.
    all_nothing_pep_cond2 = factor(case_when(
      human_protein == "NO" & unique_peptide == "YES" & intensity_sum > 0 & peptides_50_per_groups == "YES" & cond1_sum == 0 & cond2_sum >= length(cond_list[[2]])*0.5 ~ "YES",
      TRUE ~ "NO")),
    
    # Peptidos para seleccionar proteínas de interes posteriormente. No tenemos en cuenta el criterio del 50%.
    # SE USARÁ PARA FILTRAR PROTEINAS
    pep_to_select_proteins = factor(case_when(
      human_protein == "NO" & unique_peptide == "YES" & intensity_sum > 0 ~ "YES",
      TRUE ~ "NO"))
  ) %>%
  select(-cond1_sum, -cond2_sum) %>%
  ungroup()

