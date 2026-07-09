
# Purpose/Overview

The code IDOSE estimates organ-specific inhalation dose coefficients for a selected radionuclide, age, and lung absorption type adjusted for a given aerosol particle size. Total effective dose is calculated as a combination of target organ dose and tissue weighting factors. The total effective dose includes regular target organs and an estimate of dose from a group of remainder tissues. Dose is adjusted for aerosol particle size deposition in multiple regions of the respiratory tract, which affects different target organs. 

## User provided input

- Radionuclide (e.g., Sr-85)
- Lung absorption type (fast/medium/slow absorption: F/M/S)
- Age in days (e.g., 7300 ≈ adult)
- Aerosol aerodynamic diameter (default is 1.0 µm) 

## Standard program inputs

- Tissue weighting factors for organs (`namelist_regular.txt`)
- Organ masses for the “remainder tissues” calculation (`namelist_remainder.txt`)
- Aerosol size factors by respiratory region and age (`DC_PAK3.DEP`)
- Nine “region” HDB files from FGR13 with base per-organ dose coefficients for the selected nuclide/age/absorption type combination, one file per respiratory region:
  - `AI.HDB`
  - `BBE-GEL.HDB`
  - `BBE-SOL.HDB`
  - `BBE-SEQ.HDB`
  - `BBI-GEL.HDB`
  - `BBI-SOL.HDB`
  - `BBI-SEQ.HDB`
  - `ET1.HDB`
  - `ET2.HDB` 

## Outputs

- A table of target organ inhalation dose coefficients
- A table of nuclide-specific inhalation total effective dose for the fast, medium, and slow lung absorption types 

# Inhalation dose calculation approach

### 1. Initialize and select inputs

Set the nuclide, age, aerosol size, and an intake type guess if one has not been set. 

### 2. Read tissue weights and remainder masses

Load per-organ tissue weighting factors, which are used later to compute effective dose.  
Load per-organ masses used when computing the “remainder tissues” dose. 

### 3. Read aerosol size factors and pick the closest particle size bin

From `DC_PAK3.DEP`, read size-dependent factors per respiratory region.  
Find the bin whose aerosol diameter is closest to the input value and extract one size factor for each of the nine regions. 

### 4. Read region-specific organ dose coefficients

For each of the nine respiratory regions, look up the line in the corresponding HDB file that matches the nuclide, age, and intake type, and read the target organ coefficients. 

### 5. Apply aerosol scaling and combine regions

For each target organ, multiply the region’s organ coefficient by that region’s aerosol size factor and sum over all nine regions. This gives one aerosol-adjusted coefficient per organ. 

### 6. Handle gonad weighting

Use the target-organ coefficient for gonads with the appropriate tissue weighting factor contribution in the total effective dose calculation. 

### 7. Compute remainder tissues contribution

For organs classified as remainder tissues, compute the weighted contribution using the remainder tissue masses and the corresponding tissue weighting factors. 

### 8. Sum to total effective dose

Combine the regular-organ contributions and the remainder-tissues contribution to produce the total effective dose. 

### Table 1. Tissue weighting factors by age
<img width="685" height="673" alt="image" src="https://github.com/user-attachments/assets/ca697d2b-4b7a-44ab-9284-eed34c51e5e6" />
