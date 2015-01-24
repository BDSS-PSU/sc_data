rm(list = ls())
### Merge CTC_2004-2011 files

# Merge inflow
indir <- "data/CTC_Mig_2004_2011/county_in_flow/"
<<<<<<< HEAD
flist <- list.files(indir)
out <- vector(mode = "list", length = length(flist))
i <- 1
for(fname in flist) {
  print(fname)
  dat <- read.table(paste0(indir, fname), sep = ",", header = T)
  gsub("[a-z]|\\.", "", fname)
  dat$year <- 2000 + as.numeric(substring(gsub("[a-z]|\\.", "", fname), 1, 2))
  dat$flow <- "in"
  out[[i]] <- dat
  i <- i + 1
}
inflows <- do.call("rbind", out)

# Merge outflow
indir <- "data/CTC_Mig_2004_2011/county_out_flow/"
flist <- list.files(indir)
out <- vector(mode = "list", length = length(flist))
i <- 1
for(fname in flist) {
  cat(i, " ")
  print(fname)
  dat <- read.table(paste0(indir, fname), sep = ",", header = T)
  gsub("[a-z]|\\.", "", fname)
  dat$year <- 2000 + as.numeric(substring(gsub("[a-z]|\\.", "", fname), 1, 2))
  dat$flow <- "out"
  out[[i]] <- dat
  i <- i + 1
}
outflows <- do.call("rbind", out)

# Merge them both and write on file
out <- rbind(inflows, outflows)
write.table(out, file = "data/ctc_04_11.csv", sep = ",", row.names = F)