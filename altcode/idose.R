#!/usr/bin/env Rscript

# =======================================================================
# Program: IDOSE version 0.1-alpha
# Purpose: Estimate inhalation dose coefficients for selected
#          nuclide/age/intake-type with aerosol size adjustment
# =======================================================================

# -------------------------------
# Parameters
# -------------------------------
NREGION <- 9
NORG    <- 35
NAGE    <- 6
NTARGET <- 10
NDIAM   <- 50
NNUCS   <- 1252

CONVFACTOR <- 3.1536e7

indices <- c(1,2,3,4,5,6,7,8,9,10,11,20,21,22,23,24,25,26,27,28,29,30,31,32,33)

# -------------------------------
# User selections and defaults
# -------------------------------
CHEMFORM <- "INORGANIC"

# -------------------------------
# Helper functions
# -------------------------------
upcase <- function(str) {
  toupper(str)
}

read_region_hdb <- function(fname, kreg, XNUC, XAGE, XINTAKE, xnum ) {
  if (!file.exists(fname)) {
    warning(paste("Cannot open", fname))
    return(regcoef)
  }

  lines <- readLines(fname, warn = FALSE)
  if (length(lines) < 3) return(regcoef)

  # Skip first two header lines
  data_lines <- lines[-c(1,2)]

  for (ln in data_lines) {
    # Fixed-width parsing approximating:
    # nuc A7, age I5, itype A2, aerodiam E8.1, 35E10.2
    nuc      <- trimws(substr(ln, 1, 7))
    age_str  <- trimws(substr(ln, 8, 12))
    itype    <- trimws(substr(ln, 13, 14))
    # aerodiam is not used in the matching logic
    rest     <- substr(ln, 15, nchar(ln))

    age <- suppressWarnings(as.integer(age_str))
    if (is.na(age)) next
    if (age > 7300) age <- 7300

#    nums <- suppressWarnings(as.numeric(strsplit(trimws(rest), "\\s+")[[1]]))
    nums <- unlist(strsplit(trimws(rest),split=" "))

    #if (length(nums) < NORG + 1) next

    if (nuc == XNUC && age == XAGE && itype == XINTAKE) {
      print("test")
      nums <- c(nums[2:34],0.0,0.0)
      xnums <- as.numeric(unlist(nums))
      regcoef[kreg, ] <- xnums
      all <- rbind(all,xnums)
      return(regcoef)
      break
    }
  }

  #regcoef
}

# -------------------------------
# Read in nuclide list
# -------------------------------
allnucs_file <- "/home/kbaker/idose/allnucslist.txt"
if (!file.exists(allnucs_file)) stop("Cannot open allnucslist.txt")

allnucs_lines <- readLines(allnucs_file, warn = FALSE)
if (length(allnucs_lines) < 1) stop("allnucslist.txt is empty")

NUMBER <- as.integer(trimws(allnucs_lines[1]))
if (is.na(NUMBER)) stop("Could not read NUMBER from allnucslist.txt")

if (length(allnucs_lines) < (NUMBER + 1)) {
  stop("allnucslist.txt does not contain enough nuclide records")
}

# Fixed-width parsing based on the Fortran format (A7,f9.2,A2,I5)
NUCNAME <- character(NNUCS)
QAERO   <- numeric(NNUCS)
QINTAKE <- character(NNUCS)
QAGE    <- integer(NNUCS)

for (i in seq_len(NUMBER)) {
  ln <- allnucs_lines[i + 1]
  NUCNAME[i] <- trimws(substr(ln, 1, 7))
  QAERO[i]   <- suppressWarnings(as.numeric(trimws(substr(ln, 8, 16))))
  QINTAKE[i] <- trimws(substr(ln, 17, 18))
  QAGE[i]    <- suppressWarnings(as.integer(trimws(substr(ln, 19, 23))))
}

# -------------------------------
# Output files
# -------------------------------
out_comply <- "/home/kbaker/idose/output.complyinhalationdc.csv"
out_cap88  <- "/home/kbaker/idose/output.cap88v4inhalationdc.csv"

cap88_con <- file(out_cap88, open = "wt")
comp_con  <- file(out_comply, open = "wt")

on.exit({
  close(cap88_con)
  close(comp_con)
}, add = TRUE)

writeLines("IDOSE v0.1-alpha : dose coefficients provided as Sv/Bq", cap88_con)
writeLines("Nuclide,Default_Type,Age,AerosolDiam,F,M,S,QAFLAG", comp_con)

# -------------------------------
# Main loop over nuclides
# -------------------------------
for (n in seq_len(NUMBER)) {

  XNUC <- NUCNAME[n]

  # Default intake type by nuclide
  if (substr(XNUC, 1, 2) %in% c("Cs", "H-", "Se", "I-")) {
    DINTAKE <- "F"
  } else if (substr(XNUC, 1, 2) == "Th") {
    DINTAKE <- "S"
  } else {
    DINTAKE <- "M"
  }

  # Special cases
  if (substr(XNUC, 1, 2) == "Po" && CHEMFORM == "ORGANIC") {
    DINTAKE <- "F"
  } else if (substr(XNUC, 1, 2) == "Po" && CHEMFORM != "INORGANIC") {
    DINTAKE <- "M"
  }

  if (substr(XNUC, 1, 2) == "Hg" && CHEMFORM == "METHYL") {
    DINTAKE <- "M"
  }

  REMQA <- ""
  TEDE  <- numeric(3)

  # -------------------------------
  # Loop over lung absorption types
  # -------------------------------
  for (m in 1:3) {

    XINTAKE <- c("F", "M", "S")[m]

    # Initialize arrays
    all <- NULL
    REGIONCOEFFS <- matrix(0.0, nrow = NREGION, ncol = NORG)
    regcoef <- matrix(0.0, nrow = NREGION, ncol = 35)
    SIZECOEFFS   <- numeric(NREGION)
    REGWEIGHTS   <- numeric(NORG)
    REMMASS      <- numeric(NORG)
    ORGCOEFFS    <- numeric(NORG)
    XREGIONCOEFFS      <- matrix(0.0, nrow=1,ncol=NORG)
    ADIAMFACTOR  <- array(0.0, dim = c(NAGE, NTARGET, NDIAM))
    ARRAYWEIGHT  <- matrix(0.0, nrow = NORG, ncol = 7)
    ORGANNAME    <- rep(" ", NORG)

    # -------------------------------
    # Read tissue weights (regular)
    # -------------------------------
    if (!file.exists("/home/kbaker/idose/namelist_regular.txt")) stop("Cannot open namelist_regular.txt")
    reg_lines <- readLines("/home/kbaker/idose/namelist_regular.txt", warn = FALSE)
    if (length(reg_lines) < 34) stop("namelist_regular.txt too short")

    for (k in 1:33) {
      #vals <- suppressWarnings(as.numeric(strsplit(trimws(reg_lines[k + 1]), "\\s+")[[1]]))
      vals <- unlist(strsplit(reg_lines[k+1],split=","))
      if (length(vals) > 8) stop("Read error in namelist_regular.txt body")
      ORGANNAME[k]   <- trimws(substr(reg_lines[k + 1], 1, 9))
      #ARRAYWEIGHT[k,] <- vals[2:7]
      REGWEIGHTS[k]   <- as.numeric(as.character(vals[7]))
    }

    ORGANNAME[34] <- "REM"
    ORGANNAME[35] <- "E50"
    REGWEIGHTS[34:35] <- 0.0

    if (m == 1 && n == 1) {
      header_names <- c(ORGANNAME[indices], "E_50")
      writeLines(paste(header_names, collapse = ","), cap88_con)
    }

    # -------------------------------
    # Read remainder masses
    # -------------------------------
    if (!file.exists("/home/kbaker/idose/namelist_remainder.txt")) stop("Cannot open namelist_remainder.txt")
    rem_lines <- readLines("/home/kbaker/idose/namelist_remainder.txt", warn = FALSE)
    if (length(rem_lines) < 34) stop("namelist_remainder.txt too short")

    for (k in 1:33) {
      #vals <- suppressWarnings(as.numeric(strsplit(trimws(rem_lines[k + 1]), "\\s+")[[1]]))
      vals <- unlist(strsplit(rem_lines[k+1],split=","))
      if (length(vals) > 7) stop("Read error in namelist_remainder.txt body")
      REMMASS[k] <- as.numeric(as.character(vals[7]))
    }
    REMMASS[34:35] <- 0.0

    # -------------------------------
    # Read aerosol size factors
    # -------------------------------
    if (!file.exists("/home/kbaker/idose/DC_PAK3.DEP")) stop("Cannot open DC_PAK3.DEP")
    dep_lines <- readLines("/home/kbaker/idose/DC_PAK3.DEP", warn = FALSE)

    line_idx <- 1
    for (i in 1:NAGE) {
      if (line_idx > length(dep_lines)) stop("Read error: AGEGROUPNAME in DC_PAK3.DEP")
      AGEGROUPNAME <- trimws(dep_lines[line_idx])
      line_idx <- line_idx + 1

      for (j in 1:NTARGET) {
        if (line_idx > length(dep_lines)) stop("Read error: target lines in DC_PAK3.DEP")
        vals <- suppressWarnings(as.numeric(strsplit(trimws(dep_lines[line_idx]), "\\s+")[[1]]))
        line_idx <- line_idx + 1
        if (length(vals) < NDIAM + 1) stop("Read error: target lines in DC_PAK3.DEP")
        TARGET <- trimws(substr(dep_lines[line_idx - 1], 1, 9))
        ADIAMFACTOR[i, j, ] <- vals[2:(NDIAM + 1)]
      }
    }

    # -------------------------------
    # Select nearest aerosol bin
    # -------------------------------
    best <- Inf
    xref <- NA_integer_
    for (k in 1:NDIAM) {
      diff <- abs(ADIAMFACTOR[6, 1, k] - QAERO[n])
      if (diff < best) {
        best <- diff
        xref <- k
      }
    }
    if (is.na(xref)) stop("No matching aerosol diameter bin found")

    # -------------------------------
    # Build size factors for 9 regions
    # -------------------------------
    for (k in 1:NREGION) {
      SIZECOEFFS[k] <- ADIAMFACTOR[6, k + 1, xref]
    }

    # -------------------------------
    # Read 9 region HDB files
    # -------------------------------
    
    REGIONCOEFFS <- read_region_hdb("/home/kbaker/idose/AI.HDB",       1, XNUC, QAGE[n], XINTAKE ,xnums)
    all <- rbind(all,REGIONCOEFFS[1,])
    REGIONCOEFFS <- read_region_hdb("/home/kbaker/idose/BBE-GEL.HDB",  2, XNUC, QAGE[n], XINTAKE, xnums)
    all <- rbind(all,REGIONCOEFFS[2,])
    REGIONCOEFFS <- read_region_hdb("/home/kbaker/idose/BBE-SOL.HDB",  3, XNUC, QAGE[n], XINTAKE, xnums)
    all <- rbind(all,REGIONCOEFFS[3,])
    REGIONCOEFFS <- read_region_hdb("/home/kbaker/idose/BBE-SEQ.HDB",  4, XNUC, QAGE[n], XINTAKE, xnums)
    all <- rbind(all,REGIONCOEFFS[4,])
    REGIONCOEFFS <- read_region_hdb("/home/kbaker/idose/BBI-GEL.HDB",  5, XNUC, QAGE[n], XINTAKE, xnums)
    all <- rbind(all,REGIONCOEFFS[5,])
    REGIONCOEFFS <- read_region_hdb("/home/kbaker/idose/BBI-SOL.HDB",  6, XNUC, QAGE[n], XINTAKE, xnums)
    all <- rbind(all,REGIONCOEFFS[6,])
    REGIONCOEFFS <- read_region_hdb("/home/kbaker/idose/BBI-SEQ.HDB",  7, XNUC, QAGE[n], XINTAKE, xnums)
    all <- rbind(all,REGIONCOEFFS[7,])
    REGIONCOEFFS <- read_region_hdb("/home/kbaker/idose/ET1.HDB",      8, XNUC, QAGE[n], XINTAKE, xnums)
    all <- rbind(all,REGIONCOEFFS[8,])
    REGIONCOEFFS <- read_region_hdb("/home/kbaker/idose/ET2.HDB",      9, XNUC, QAGE[n], XINTAKE, xnums)
    all <- rbind(all,REGIONCOEFFS[9,])

    REGIONCOEFFS <- all
    
    # -------------------------------
    # Aggregate per-organ coefficients
    # -------------------------------
    for (i in 1:NORG) {
      ORGCOEFFS[i] <- sum(REGIONCOEFFS[, i] * SIZECOEFFS)
    }

    # -------------------------------
    # Gonad weighting selection
    # -------------------------------
    idxov <- 21
    idxte <- 26
    WGONOV <- 0.0
    WGONTE <- 0.0
    if (ORGCOEFFS[idxov] >= ORGCOEFFS[idxte]) {
      WGONOV <- 0.2
    } else {
      WGONTE <- 0.2
    }

    # -------------------------------
    # Sum regular component
    # -------------------------------
    SUMREG <- 0.0
    for (i in 1:NORG) {
      if (i == idxov) {
        SUMREG <- SUMREG + ORGCOEFFS[i] * WGONOV
      } else if (i == idxte) {
        SUMREG <- SUMREG + ORGCOEFFS[i] * WGONTE
      } else {
        SUMREG <- SUMREG + ORGCOEFFS[i] * REGWEIGHTS[i]
      }
    }

    # Max regular coefficient
    MAXREGC <- 0.0
    for (i in 1:NORG) {
      if (REGWEIGHTS[i] > 0.0 && ORGCOEFFS[i] > MAXREGC) {
        MAXREGC <- ORGCOEFFS[i]
      }
    }
    if (WGONOV > 0.0 && ORGCOEFFS[idxov] > MAXREGC) MAXREGC <- ORGCOEFFS[idxov]
    if (WGONTE > 0.0 && ORGCOEFFS[idxte] > MAXREGC) MAXREGC <- ORGCOEFFS[idxte]

    # -------------------------------
    # Remainder dose (HREM)
    # -------------------------------
    SUMREM   <- 0.0
    SUMREMW  <- 0.0
    MAXREMC  <- 0.0
    idxmaxrem <- 0

    for (i in 1:NORG) {
      if (REMMASS[i] > 0.0) {
        SUMREM  <- SUMREM  + ORGCOEFFS[i] * REMMASS[i]
        SUMREMW <- SUMREMW + REMMASS[i]
        if (ORGCOEFFS[i] > MAXREMC) {
          MAXREMC <- ORGCOEFFS[i]
          idxmaxrem <- i
        }
      }
    }

    if (SUMREMW <= 0.0) {
      HREM <- 0.0
      REMQA <- paste0(REMQA, "Z")
    } else {
      if (MAXREGC > MAXREMC) {
        HREM <- SUMREM / SUMREMW
        REMQA <- paste0(REMQA, "A")
      } else {
        SUMREM2  <- 0.0
        SUMREMW2 <- 0.0
        for (i in 1:NORG) {
          if (REMMASS[i] > 0.0 && i != idxmaxrem) {
            SUMREM2  <- SUMREM2  + ORGCOEFFS[i] * REMMASS[i]
            SUMREMW2 <- SUMREMW2 + REMMASS[i]
          }
        }
        if (SUMREMW2 > 0.0) {
          HREM <- 0.5 * ((SUMREM2 / SUMREMW2) + MAXREMC)
          REMQA <- paste0(REMQA, "B")
        } else {
          HREM <- MAXREMC
          REMQA <- paste0(REMQA, "C")
        }
      }
    }

    # -------------------------------
    # Effective dose
    # -------------------------------
    HEFF <- SUMREG + 0.05 * HREM
    ORGCOEFFS[34] <- HREM
    ORGCOEFFS[35] <- HEFF
    TEDE[m] <- HEFF

    # CAP88 output only when intake types match
    if (XINTAKE == QINTAKE[n]) {
      row_vals <- c(NUCNAME[n], QINTAKE[n], as.character(QAGE[n]), sprintf("%.2f", QAERO[n]))
      # write line 1
      cat(paste(row_vals, collapse = ","), file = cap88_con, sep = "\n", append = TRUE)

      coeff_text <- sprintf("%.3E", ORGCOEFFS)
      cap88_row <- paste(c(coeff_text[indices], sprintf("%.3E", HEFF)), collapse = ",")
      cat(cap88_row, file = cap88_con, sep = "\n", append = TRUE)
    }
  }

  # -------------------------------
  # COMPLY output
  # -------------------------------
  comply_line <- paste(
    NUCNAME[n],
    DINTAKE,
    QAGE[n],
    sprintf("%.2f", QAERO[n]),
    sprintf("%.4E", TEDE[1]),
    sprintf("%.4E", TEDE[2]),
    sprintf("%.4E", TEDE[3]),
    REMQA,
    sep = ","
  )
  cat(comply_line, file = comp_con, sep = "\n", append = TRUE)
}

message("Done.")

