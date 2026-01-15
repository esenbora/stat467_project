# STAT 467 - Population Data Correction
# Purpose: Replace inconsistent Population values with World Bank data
# Run this BEFORE 00_data_preprocessing.R
# Input: data.csv | Output: data_population_corrected.csv

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Load required packages
if (!require("WDI")) install.packages("WDI")
if (!require("tidyverse")) install.packages("tidyverse")

library(WDI)
library(tidyverse)

cat("\n=== Population Data Correction ===\n\n")

# --- Section 1: Country Name Mapping ---
# WHO uses long-form names, World Bank uses different conventions
country_mapping <- data.frame(
  who_name = c(
    "Bahamas",
    "Bolivia (Plurinational State of)",
    "Cabo Verde",
    "Congo",
    "Czechia",
    "Côte d'Ivoire",
    "Democratic People's Republic of Korea",
    "Democratic Republic of the Congo",
    "Egypt",
    "Gambia",
    "Iran (Islamic Republic of)",
    "Kyrgyzstan",
    "Lao People's Democratic Republic",
    "Micronesia (Federated States of)",
    "Republic of Korea",
    "Republic of Moldova",
    "Russian Federation",
    "Saint Kitts and Nevis",
    "Saint Lucia",
    "Saint Vincent and the Grenadines",
    "Slovakia",
    "Somalia",
    "Swaziland",
    "Syrian Arab Republic",
    "The former Yugoslav republic of Macedonia",
    "United Kingdom of Great Britain and Northern Ireland",
    "United Republic of Tanzania",
    "United States of America",
    "Venezuela (Bolivarian Republic of)",
    "Viet Nam",
    "Yemen"
  ),
  wb_name = c(
    "Bahamas, The",
    "Bolivia",
    "Cabo Verde",
    "Congo, Rep.",
    "Czechia",
    "Cote d'Ivoire",
    "Korea, Dem. People's Rep.",
    "Congo, Dem. Rep.",
    "Egypt, Arab Rep.",
    "Gambia, The",
    "Iran, Islamic Rep.",
    "Kyrgyz Republic",
    "Lao PDR",
    "Micronesia, Fed. Sts.",
    "Korea, Rep.",
    "Moldova",
    "Russian Federation",
    "St. Kitts and Nevis",
    "St. Lucia",
    "St. Vincent and the Grenadines",
    "Slovak Republic",
    "Somalia, Fed. Rep.",
    "Eswatini",
    "Syrian Arab Republic",
    "North Macedonia",
    "United Kingdom",
    "Tanzania",
    "United States",
    "Venezuela, RB",
    "Viet Nam",
    "Yemen, Rep."
  ),
  stringsAsFactors = FALSE
)

# --- Section 2: Fetch World Bank Data ---
cat("Fetching World Bank population data (2000-2015)...\n")
pop_wb <- WDI(
  indicator = "SP.POP.TOTL",
  country = "all",
  start = 2000,
  end = 2015,
  extra = FALSE
)

cat("Downloaded", nrow(pop_wb), "rows from World Bank\n")

# Clean World Bank data
pop_wb_clean <- pop_wb %>%
  select(country, year, pop_wb = SP.POP.TOTL) %>%
  filter(!is.na(pop_wb))

# --- Section 3: Load WHO Data ---
cat("\nLoading WHO data...\n")
who_raw <- read.csv("data.csv", stringsAsFactors = FALSE, check.names = FALSE)
colnames(who_raw) <- trimws(colnames(who_raw))

cat("WHO data:", nrow(who_raw), "rows,", length(unique(who_raw$Country)), "countries\n")
cat("Original Population coverage:", sum(!is.na(who_raw$Population)), "/", nrow(who_raw), "\n")

# --- Section 4: Apply Country Mapping ---
who_mapped <- who_raw %>%
  left_join(country_mapping, by = c("Country" = "who_name")) %>%
  mutate(Country_match = coalesce(wb_name, Country)) %>%
  select(-wb_name)

# --- Section 5: Merge Population Data ---
cat("\nMerging with World Bank population...\n")
who_corrected <- who_mapped %>%
  left_join(
    pop_wb_clean,
    by = c("Country_match" = "country", "Year" = "year")
  ) %>%
  mutate(
    Population_old = Population,
    Population = coalesce(pop_wb, Population)
  ) %>%
  select(-Country_match, -pop_wb)

# --- Section 6: Validation ---
cat("\n=== Merge Results ===\n")
cat("Total rows:", nrow(who_corrected), "\n")
cat("Population filled (new):", sum(!is.na(who_corrected$Population)), "\n")
cat("Still missing:", sum(is.na(who_corrected$Population)), "\n")

# Show unmatched countries
unmatched <- who_corrected %>%
  filter(is.na(Population)) %>%
  distinct(Country) %>%
  pull(Country)

if (length(unmatched) > 0) {
  cat("\nCountries without World Bank match:\n")
  cat(paste(" -", unmatched), sep = "\n")
}

# Spot check: Show some corrected values
cat("\n=== Sample Corrections ===\n")
sample_check <- who_corrected %>%
  filter(!is.na(Population_old) & !is.na(Population)) %>%
  filter(abs(Population - Population_old) > 1000) %>%
  mutate(ratio = Population / Population_old) %>%
  filter(ratio > 10 | ratio < 0.1) %>%
  select(Country, Year, Population_old, Population, ratio) %>%
  head(10)

if (nrow(sample_check) > 0) {
  print(sample_check)
} else {
  cat("No major corrections detected (values were already consistent)\n")
}

# Spot check: China 2010 should be ~1.34 billion
china_2010 <- who_corrected %>%
  filter(Country == "China", Year == 2010) %>%
  pull(Population)

cat("\nSpot check - China 2010 population:", format(china_2010, big.mark = ","), "\n")

# --- Section 7: Remove countries without Population data ---
# Cook Islands and Niue have no World Bank data - exclude from analysis
who_corrected <- who_corrected %>%
  filter(!is.na(Population))

cat("\nRemoved countries without Population data.\n")
cat("Final dataset:", nrow(who_corrected), "rows\n")

# --- Section 8: Save Output ---
# Remove the old Population column we kept for comparison
who_final <- who_corrected %>%
  select(-Population_old)

write.csv(who_final, "data_population_corrected.csv", row.names = FALSE)
cat("\n=== Saved: data_population_corrected.csv ===\n")
cat("Next step: Run 00_data_preprocessing.R\n")
