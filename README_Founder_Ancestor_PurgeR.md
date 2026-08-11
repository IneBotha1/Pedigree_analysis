# README – Founder and Ancestor Analysis
## South African Boerperd

## Purpose

This analysis was used to describe the **founder and ancestor structure of the South African Boerperd pedigree** and to assess how concentrated the recorded genetic origins of the reference population are.

The analysis was performed in **PurgeR** and was followed by a simple birth-year analysis of the actual founders and ancestors identified in the pedigree.

---

See the scripts `Script_local_PURGER FOUNDER AND ANCESTOR ANALYSIS.R` and `Script_local_visualize_Founder_Ancestor_birthyear.R`
## Input pedigree

The PurgeR analysis used:

`SAH_ped_PurgeR_input_loops_res_2025_excl.csv`

The pedigree had already been cleaned and pedigree loops had been resolved.

The main columns used were:

- `ANI_ID`
- `ANI_ID_SIRE`
- `ANI_ID_DAM`
- `BIRTH_DTM`
- `SEX`
- `STATUS`

Animals born in 2025 were excluded from this analysis.

---

## Reference population

The PurgeR reference population was defined from the `STATUS` column.

Animals whose status began with:

`REG`

were treated as members of the reference population.

This included status values such as:

- `REG`
- `REG (Registered)`

The reference flag was created using:

```r
ped_raw$reference <- (
  !is.na(ped_raw$STATUS) &
    grepl("^REG\\b", trimws(ped_raw$STATUS), ignore.case = TRUE)
)
```

---

## Pedigree preparation

The pedigree was prepared using:

```r
ped <- purgeR::ped_sort(
  ped_raw,
  id = "ANI_ID",
  dam = "ANI_ID_DAM",
  sire = "ANI_ID_SIRE",
  keep_names = TRUE
)
```

`ped_sort()` renumbered animals internally for PurgeR while retaining the original animal IDs.

---

## Founder and ancestor analysis

Founder- and ancestor-based parameters were calculated with:

```r
set.seed(1973)

ancestor_results <- purgeR::pop_Nancestors(
  ped = ped,
  reference = "reference",
  nboot = 10000,
  seed = 1973,
  skip_Ng = FALSE
)
```

The main parameters were:

### Nf – Number of founders

The number of pedigree founders contributing to the selected reference population.

In the final PurgeR-consistent reconstruction:

**Nf = 587 founders**

### Na – Number of ancestors

The number of actual ancestral animals included in the PurgeR ancestor set for the reference population.

In the final reconstruction:

**Na = 3,598 ancestors**

### Nfe – Effective number of founders

The number of equally contributing founders that would produce the same concentration of founder contributions observed in the reference population.

### Nae – Effective number of ancestors

The number of equally contributing ancestors that would produce the same concentration of ancestor contributions observed in the reference population.

### Ng – Founder genome equivalents

A founder-origin diversity measure that additionally reflects the loss of founder genetic material through drift and segregation.

The exact `Nfe`, `Nae` and `Ng` values should be taken directly from the saved PurgeR results rather than typed manually into this README.

---

## Important distinction

`Nf` and `Na` refer to **actual identifiable animals** in the pedigree.

`Nfe`, `Nae` and `Ng` are **effective-number statistics**.

Therefore:

- actual founders can be listed individually;
- actual ancestors can be listed individually;
- birth dates can be attached to those animals where available;
- there is no literal group of individual animals corresponding to `Nfe` or `Nae`.

---

# Founder and Ancestor ID Reconstruction

The actual founder and ancestor IDs were reconstructed from the PurgeR-sorted pedigree so that their original animal IDs and birth dates could be recovered.

A recursive parent trace initially identified **3,574 ancestral animals**.

However, PurgeR reported:

**Na = 3,598**

The difference was **24 animals**.

These 24 animals were members of the registered reference population that were themselves pedigree founders. Because they had no recorded parents, they were not encountered by a parent-only recursive ancestry trace.

They were therefore added to the reconstructed ancestor set.

After doing this:

- reconstructed founders = **587**, matching PurgeR `Nf`;
- reconstructed ancestors = **3,598**, matching PurgeR `Na`.

This cross-check was used to confirm that the extracted individual IDs were consistent with the PurgeR results.

---

# Founder and Ancestor Birth-Year Analysis

After identifying the actual founder and ancestor animals, their recorded birth dates were matched back to the original pedigree.

The purpose was to examine **when the recorded founders and ancestors entered the pedigree** and how much historical birth-year information was available.

---

## Founder birth years

Total founders:

**587**

Birth year available:

**320 founders**

Birth year missing:

**267 founders**

Therefore approximately:

- **54.5%** of founders had a known birth year;
- **45.5%** had an unknown birth year.

Among founders with known birth years, the recorded range was:

**1948–2015**

Because nearly half of the founders have missing birth years, the temporal founder distribution should be interpreted cautiously.

---

## Ancestor birth years

Total ancestors:

**3,598**

Birth year available:

**3,302 ancestors**

Birth year missing:

**296 ancestors**

Therefore approximately:

- **91.8%** of ancestors had a known birth year;
- **8.2%** had an unknown birth year.

Among ancestors with known birth years, the recorded range was:

**1948–2018**

The ancestor birth-year distribution is therefore substantially more complete than the founder birth-year distribution.

---

## Important terminology

The founder animals are also part of the broader ancestor set.

For figures and tables it is therefore clearer to use:

- **Founders**
- **All ancestors**

rather than implying that the two groups are completely separate.

If a non-founder ancestor group is required, it can be defined as:

```r
setdiff(
  ancestor_ped_ids_purger,
  founder_ped_ids
)
```

---

## Birth-year figures

The temporal distributions were summarised using the animals with known birth years.

The planned/main figure used:

- founders;
- all ancestors;
- birth year on the x-axis;
- percentage of animals with known birth years on the y-axis;
- 5-year birth-year bins; and
- historical reference lines at **1973** and **1998**.

A cumulative birth-year figure was also considered to show the proportion of founders and ancestors recorded by successive years.

Unknown birth years were excluded from these temporal distributions but their numbers were reported separately.

---

## Main output files

The analysis produced or was designed to retain:

- `PurgeR_REG_founder_ancestor_results.csv`
- `PurgeR_REG_founders_birthdates.csv`
- `PurgeR_REG_ancestors_birthdates.csv`
- `PurgeR_REG_founder_ancestor_birthdate_summary.csv`
- founder and ancestor birth-year distribution tables
- founder and ancestor birth-year figures
- the R script used for the PurgeR analysis

---

## Short Methods description

Founder- and ancestor-based pedigree parameters were estimated in PurgeR for registered South African Boerperd animals. The pedigree was sorted using `ped_sort()`, and `pop_Nancestors()` was used to estimate the number of founders (`Nf`), effective number of founders (`Nfe`), number of ancestors (`Na`), effective number of ancestors (`Nae`) and founder genome equivalents (`Ng`). The final analysis identified 587 founders and 3,598 ancestors. Original animal IDs were reconstructed from the PurgeR pedigree and linked to recorded birth dates. Birth years were available for 320 of 587 founders and 3,302 of 3,598 ancestors, and temporal distributions were evaluated using animals with known birth years.

---

## Key point

> `Nf` and `Na` describe identifiable pedigree animals, whereas `Nfe`, `Nae` and `Ng` describe the effective concentration and retention of founder/ancestor genetic contributions.

The birth-year analysis was added to place the actual founder and ancestor animals in their historical pedigree context while explicitly retaining the large amount of missing founder birth-year information.
