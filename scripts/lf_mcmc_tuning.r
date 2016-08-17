library(zoo)
par(mfrow=c(5,2) ,oma = c(0, 0, 2, 0) )
par(mar=c(4,2,2,4))


for(subfolder in list.dirs(full.names = FALSE))
{
  if (subfolder != "")
  {
  setwd(sprintf("C:/cygwin64/home/Gregg/mcmc/output/%s",subfolder))
  list_of_files = list.files(pattern="LangevinTimeSeries_.*.txt")
  len = length(list_of_files)
  a = numeric(len)
  b = numeric(len)
  c = numeric(len)
  
  for(i in 1:len)
  {
    mcmc_data = read.table(list_of_files[i])
    x <- head( readLines(list_of_files[i]))
    meta_data_string <- grep( pattern=subfolder, x, value=TRUE)
    real_parts = mcmc_data[7:7]
    imag_parts = mcmc_data[8:8]
    a[i] = var(as.numeric(unlist(real_parts)))
    b[i] = var(as.numeric(unlist(imag_parts)))
    c[i] = as.numeric( gsub("[^0-9.]", "", meta_data_string) )
  }
  
  plot(c, a, log='y', xlab=subfolder, ylab="var", main="Re")
  plot(c, b, log='y', xlab=subfolder, ylab="var", main="Im")
  #mtext("Parameter proposal variance vs. posterior variance.", outer=TRUE, cex=1.5)
  }
}