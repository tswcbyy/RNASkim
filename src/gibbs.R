library("locfit");
dat <- read.table("rsgibbs_locfit.dat", header=TRUE);
fit <- locfit::locfit(w~m, data=dat);
result <- predict(fit, dat[,"m"]);
dat$r <- result;
dat$w <- NULL;
write.table(dat, "rsgibbs_locfit.dat", append=FALSE, sep="\t", row.names=FALSE, col.names=FALSE); 
