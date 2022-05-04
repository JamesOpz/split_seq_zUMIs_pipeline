## Generate list of shared barcodes for split_3 analysis
# JO 2021/04/12

# Libraries
library(stringr)
library(readxl)
library(readr)

# Load bc files that are provided to IDT
rt_bc_file <- "../rt_primers_idt_order.xlsx"

# Read files
rt_oligo_seq <- read_xlsx(rt_bc_file)

# Define bc positions within the IDT formatted sequence for each file
rt_bc_pos <- c(23,30)

bc_length <- 8
bc_number <- 96

# Extract BCs from IDT provided sequence
# RT
rt_bc_list_split <- str_split(rt_oligo_seq$Sequence, "")
rt_bc <-  sapply(rt_bc_list_split,
                 FUN = function(x) paste(x[rt_bc_pos[1]:rt_bc_pos[2]], collapse = ""))

# Create the sharing datafram rows A1 (1) + E1 (49), etc.
rt_sharing_df <- data.frame(polyA = rt_bc[1:48], rand_hex = rt_bc[49:96])

# Write map to file
rt_sharing_map <- data.frame(well_position = rt_oligo_seq$WellPosition[1:48], 
                             polya = rt_sharing_df$polyA, 
                             rand_hex = rt_sharing_df$rand_hex)

write_csv(rt_sharing_map, file = "../RT_bc_map.csv")

# cocatenate with a tab
rt_bc_sharing = paste(rt_sharing_df$polyA, rt_sharing_df$rand_hex, sep="\t")

# Write to file - Remember that the #BC position needs to be inserted into the top manually, 
# check previous file or zUMI wiki
write_lines(rt_bc_sharing, file = "../rt_bc_sharing.txt")
