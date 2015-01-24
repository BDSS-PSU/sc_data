### Merge CTC_2004-2011 files

# Merge inflow
indir <- "data/CTC_Mig_2004_2011/county_in_flow/"
files <- list.files(indir)
for(fname in flist) {
  dat <- read.table(paste0(indir, fname), sep = ",", header = T)
  
  gsub()
  dat$year <- 
}