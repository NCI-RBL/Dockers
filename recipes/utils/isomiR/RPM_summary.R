library(dplyr)
args = commandArgs(trailingOnly=TRUE)

#read group_output.isomir.tsv from quagmiR
group_output <- read.csv(file = args[1], sep = '\t')

#calculate RPM (Reads per million total reads) and transpose table
group_output <- group_output %>% 
  mutate(RPM=TOTAL_READS/TOTAL_READS_IN_SAMPLE*1000000) %>%
  select(SAMPLE, MIRNA, RPM) %>%
  pivot_wider(names_from = SAMPLE, values_from = RPM, values_fill = 0.01)

write.csv(group_output, args[2], row.names = F)
