# README – Pedigree Completeness Analysis
## South African Boerperd

## Purpose

This analysis was used to assess the **quality and depth of the South African Boerperd pedigree** before interpreting pedigree-based measures such as inbreeding, founder contributions and population structure.

Two complementary views of pedigree completeness were used:

1. **Overall pedigree depth and completeness** using `optiSel`.
2. **Change in pedigree completeness through time** using the annual completeness output from POPREP.

---

## Pedigree used

The analysis used the complete, loop-resolved pedigree:

`SAH_basic_ped_input_ungrouped_loopsres.csv`

The complete pedigree was retained so that all available ancestral information could contribute to the completeness calculations.

---

## 1. Overall pedigree completeness with optiSel

The pedigree was prepared in R using `optiSel::prePed()`.

Pedigree depth and completeness were then obtained using `summary.Pedig()`.

The main settings were:

```r
summary.Pedig(
  Pedig,
  maxd = 50,
  d = 4
)
```

### Why these settings were used

`maxd = 50` gives the function a sufficiently large maximum depth to search through the pedigree. It is only a computational tracing limit and does **not** mean that 50 generations were observed.

`d = 4` was used for the pedigree completeness index (PCI), meaning that PCI was calculated across four ancestral generations.

---

## Main completeness measures

### `fullGen`

The number of generations for which **all ancestors are known**.

A higher value indicates that complete pedigree information extends further backwards.

### `maxGen`

The deepest generation in which at least one ancestor can be traced.

This describes the maximum recorded pedigree depth and should not be confused with fully complete generations.

### `equiGen`

Equivalent complete generations.

This combines the proportion of known ancestors across all generations into a single measure. It can be interpreted as the amount of pedigree information equivalent to a certain number of completely known generations.

### `PCI`

Pedigree Completeness Index.

This summarises how complete the pedigree is across the specified ancestral depth. In this analysis, PCI was calculated across four generations.

---

## Overall results

For the complete South African Boerperd pedigree:

- Mean fully traced generations (`fullGen`) = **3.95**
- Deepest recorded pedigree path = **12 generations**
- Mean equivalent complete generations (`equiGen`) = **8.00**
- Mean four-generation PCI = **0.74**

These measures were used to describe the overall quality and depth of the pedigree.

The value of 12 generations refers to the **deepest recorded ancestral path**, not the average pedigree depth.

---

## 2. Completeness by ancestral generation

Pedigree completeness was also examined across successive ancestral generations.

For each generation, the proportion of known ancestors was calculated and plotted.

The purpose of this figure was to show how pedigree information declines as ancestral depth increases.

This is useful because a pedigree can appear deep based on `maxGen` while still becoming incomplete at deeper generations.

---

## 3. Annual pedigree completeness from POPREP

A second analysis examined pedigree completeness by **birth year**.

The POPREP output contained:

- Year
- Number of animals
- Complete generation 1
- Complete generation 2
- Complete generation 3
- Complete generation 4
- Complete generation 5
- Complete generation 6

The annual data were imported from:

`POPREP inbreeding TABLE 1 year number animals pedigree completeness.xlsx`

Only years up to and including **2024** were retained.

---

## Annual completeness graph

The final graph showed:

- six coloured lines representing pedigree completeness from generations 1 to 6;
- pedigree completeness (%) on the left y-axis;
- year on the x-axis; and
- a solid black line showing the number of animals recorded per year on the right y-axis.

This allowed changes in pedigree-recording depth to be viewed together with changes in the number of animals represented in each birth year.

---

## Why both approaches were used

The two analyses answer slightly different questions.

**Overall optiSel metrics** describe the general depth and completeness of the pedigree.

**Annual POPREP completeness** shows how pedigree recording changed through time.

Using both therefore gives a more complete picture than reporting only one pedigree-completeness statistic.

---

## Main outputs to keep

- Overall `fullGen`, `maxGen`, `equiGen` and PCI results
- Pedigree completeness-by-generation figure
- Annual generations 1–6 completeness table
- Annual pedigree completeness + number-of-animals figure
- R scripts used for the analysis

---

## Short Methods description

Pedigree completeness was evaluated using the complete loop-resolved South African Boerperd pedigree. The pedigree was prepared using the R package `optiSel`, and pedigree depth was summarised using the number of fully traced generations (`fullGen`), maximum traced generation (`maxGen`), equivalent complete generations (`equiGen`) and the pedigree completeness index (PCI). PCI was calculated across four ancestral generations. In addition, annual pedigree completeness from one to six ancestral generations was obtained from POPREP output and plotted by birth year together with the number of animals recorded per year.

---

## Key point

> Pedigree completeness was assessed both as **overall pedigree depth** and as **change in recording quality through time**.

These results provide the pedigree-quality context needed before interpreting later pedigree-based genetic diversity analyses.
