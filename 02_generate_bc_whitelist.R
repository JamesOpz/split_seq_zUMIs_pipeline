### Generate all BC combinations from for split_3 analysis
### JO 12/04/2021

# Libraries
library(stringr)
library(readxl)
library(readr)

# Load bc files that are provided to IDT
rt_bc_file <- "../rt_primers_idt_order.xlsx"
r2_bc_file <- "../r1_ligation_idt_order.xlsx"
r3_bc_file <- "../r2_ligation_idt_order.xlsx"

# Read files
rt_oligo_seq <- read_xlsx(rt_bc_file)
r2_oligo_seq <- read_xlsx(r2_bc_file)
r3_oligo_seq <- read_xlsx(r3_bc_file)

# Define bc positions within the IDT formatted sequence for each file
rt_bc_pos <- c(23,30)
r2_bc_pos <- c(23,30)
r3_bc_pos <- c(41,48)

bc_length <- 8
bc_number <- 96

# Extract BCs from IDT provided sequence
# RT
rt_bc_list_split <- str_split(rt_oligo_seq$Sequence, "")
rt_bc <-  sapply(rt_bc_list_split,
                 FUN = function(x) paste(x[rt_bc_pos[1]:rt_bc_pos[2]], collapse = ""))

# R2 (ligation)
r2_bc_spit <- str_split(r2_oligo_seq$Sequence, "")
r2_bc <-  sapply(r2_bc_spit,
                 FUN = function(x) paste(x[r2_bc_pos[1]:r2_bc_pos[2]], collapse = ""))

# R3 (ligation)
r3_bc_spit <- str_split(r3_oligo_seq$Sequence, "")
r3_bc <-  sapply(r3_bc_spit,
                 FUN = function(x) paste(x[r3_bc_pos[1]:r3_bc_pos[2]], collapse = ""))

# Checks for correct extraction
# RT
length(rt_bc) == bc_number
nchar(rt_bc[1]) == bc_length
# R2
length(r2_bc) == bc_number
nchar(r2_bc[1]) == bc_length
#R3
length(r3_bc) == bc_number
nchar(r3_bc[1]) == bc_length

# Write barcodes to file
# Reconstruct map
r2_bc_map <- data.frame(well_position = r2_oligo_seq$WellPosition, bc = r2_bc)
r3_bc_map <- data.frame(well_position = r3_oligo_seq$WellPosition, bc = r3_bc)

# Write
write_csv(r2_bc_map, file = "../R2_bc_map.csv")
write_csv(r3_bc_map, file = "../R3_bc_map.csv")

# Generate the whitelist of all BC combinations
# The reads are sequenced outside in, so barcodes are (UMI) 3 -> 2 -> 1 PloyA seq
bc_grid = expand.grid(r3_bc, r2_bc, rt_bc)

# cocatenate
bc_wlist = paste(bc_grid$Var1, bc_grid$Var2, bc_grid$Var3, sep="")

# Write the expected barcode whitelist out to file
write_lines(bc_wlist, file = "../split_seq_barcode_wlist.txt")

