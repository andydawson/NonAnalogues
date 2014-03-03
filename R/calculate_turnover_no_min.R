#  Calculate the dissimilarities for all sites in the large compiled.pollen data.frame.
#  We want to calculate several metrics in two directions, forward and backwards,
#  and we ideally want to keep track of where the closest landscape points are.
#  Landscape similarity gives us a sense of how close a site is to an ecosystem that
#  has previously existed on the landscape.  We also look forward to see

library(snowfall)
sfStop()
sfInit(parallel = TRUE, cpus = 4)

#  If 'Other' accounts for more than 10% of the total pollen sum we're going 
#  to exclude the sample.  This excludes fully a quarter of the sites in the dataset!
#  n = 5441 of 21397

data(pollen.equiv)

pol.types <- c('Other', unique(as.character(pollen.equiv$WhitmoreSmall)))
pol.types <- colnames(compiled.pollen)[colnames(compiled.pollen)%in%pol.types[!is.na(pol.types)]]

cp.pct <- compiled.pollen[,pol.types]/rowSums(compiled.pollen[,pol.types])

foo = sqrt(sum(cp.pct[1,]^2 - cp.pct[2,]^2))

no.others <- cp.pct$Other > 0.10

rep.frame <- data.frame(site = compiled.pollen$sitename,
                        dataset = compiled.pollen$dataset,
                        age = compiled.pollen$age,
                        self.min.past = rep(NA, nrow(compiled.pollen)),
                        self.min.futu = rep(NA, nrow(compiled.pollen)),
                        land.min.past = rep(NA, nrow(compiled.pollen)),
                        land.min.futu = rep(NA, nrow(compiled.pollen)),
                        land.pts.past = rep(NA, nrow(compiled.pollen)),
                        land.pts.futu = rep(NA, nrow(compiled.pollen)),
                        sample.size = rep(NA, nrow(compiled.pollen)),
                        self.size = rep(NA, nrow(compiled.pollen)),
                        matrix(nrow=nrow(compiled.pollen), ncol=100))

landscape.dat <- data.frame(site = compiled.pollen$sitename,
                            dataset = compiled.pollen$dataset,
                            age = compiled.pollen$age)

# # first make some rank plots
# i=1
# right.site <- compiled.pollen$sitename == rep.frame$site[i]
# site.ages <-  compiled.pollen$age[right.site]
# 
# # min dis for the right site
# site.samples <- cp.pct[right.site, ]#[right.age[-i] & right.site[-i],]      
# 
# 
# init_sorted = sort(site.samples[1,], decreasing=TRUE)
# counts_sorted = site.samples[,names(init_sorted)]
# 
# for (j in 1:nrow(counts_sorted)){
#   counts_yr = counts_sorted[j,]
#   
#   #   counts_yr = sort(counts_yr, decreasing=TRUE)
#   x = seq(1:length(counts_yr))
#   
#   plot(x, log(counts_yr), main=j)
# }
# 
# 
# self.minus <- apply(self.samples, 1, function(x) (arrow.mat - x)^2)
# self.vals  <- min(c(1000, sqrt(colSums(self.minus, na.rm=TRUE))))


#  Calculate the dissimilarities for all sites in the large compiled.pollen data.frame.







for(i in 10:30){#nrow(rep.frame)){
  
  if(any(is.na(rep.frame[i, 4:103])) & (no.others[i] == FALSE)){
    #  For each sample in the dataset we need to find it, and then check if it
    #  has any samples that are between 250 and 750 years older than it.
    right.site <- compiled.pollen$sitename == rep.frame$site[i]
    site.ages <-  compiled.pollen$age[right.site]
    site.lat  <-  unique(compiled.pollen$lat[right.site])
    site.long  <-  unique(compiled.pollen$long[right.site])
    site.coords <- c(site.long, site.lat)
    
    right.age.past <- ((compiled.pollen$age > (rep.frame$age[i] + 250)) & 
                    (compiled.pollen$age < (rep.frame$age[i] + 750)))
    right.age.futu <- ((compiled.pollen$age > (rep.frame$age[i] - 250)) & 
                     (compiled.pollen$age < (rep.frame$age[i] - 750)))
                
    right.age = right.age.past
    
#     if(any((right.age.past & right.site) & sum(right.age.past) > 5) &
#        any((right.age.futu & right.site) & sum(right.age.futu) > 5)){
      
      #  Now we know that the sample has something to compare to (that is 
      #  between 250 and 750 years older), we can create a vector for the sample.
      #  We exclude 'Other', because it's super big some times.
      
      arrow <- cp.pct[i,]
      arrow.mat <- as.matrix(arrow)
      
      # what are the calib?
      calib.samples <- cp.pct[-i, ][right.age[-i],]
      calib.sites <- compiled.pollen[-i, ][right.age[-i],][,1]
      calib.minus <- apply(calib.samples, 1, function(x) (arrow.mat - x)^2)
      calib.vals <- sqrt(colSums(calib.minus, na.rm=TRUE))
      
      # all sites (excluding right site) with samples in the 250 - 750 years 
      # prior to the right age
      landscape.samples <- cp.pct[-i, ][right.age[-i] & !right.site[-i],]
      landscape.minus <- apply(landscape.samples, 1, function(x) (arrow.mat - x)^2)
      landscape.vals  <- sqrt(colSums(landscape.minus, na.rm=TRUE))
      
      landscape.meta   <- compiled.pollen[-i, ][right.age[-i] & !right.site[-i],]
      landscape.age    <- landscape.meta$age
      landscape.coords <- cbind(landscape.meta$long, landscape.meta$lat) 
      
      landscape.dist   <- rdist.earth(matrix(site.coords, nrow=1), matrix(landscape.coords, ncol=2))
      
      # min dis for the right site
      self.samples <- cp.pct[-i, ][right.age[-i] & right.site[-i],]      
      self.minus <- apply(self.samples, 1, function(x) (arrow.mat - x)^2)
      self.vals  <- min(c(1000, sqrt(colSums(self.minus, na.rm=TRUE))))
      
      par(ask=TRUE)
      plot(landscape.dist, landscape.vals, ylim=c(0,max(landscape.vals)), pch=20, 
           col=ifelse(abs(landscape.age - rep.frame$age[i]) < 300, "blue", "gray23"))#col='gray23')
      points(0, self.vals, col='red', pch=8)
      abline(h=self.vals, col='black', lty=2)
      
      plot(landscape.age, landscape.vals, ylim=c(0,max(landscape.vals)), pch=20, col='gray23')
      
      landscape.vals.min  <- min(c(1000, sqrt(colSums(landscape.minus, na.rm=TRUE))))
      

      
      sfExport("calib.sites")
      sfExport("calib.samples")
      sfExport('calib.vals')
      
      min.dist <- function(x){
        resampled <- sample(nrow(calib.samples), replace=TRUE)
        dist.test <- calib.vals[resampled][!duplicated(calib.sites[resampled])]
        min(dist.test)
      }
      
      rep.frame[i,12:ncol(rep.frame)] <- unlist(sfLapply(1:100, min.dist))
      rep.frame$self.min.past[i] <- self.vals
      rep.frame$land.min.past[i] <- landscape.vals.min
#       rep.frame$min.dist[i] <- min(calib.vals)
      rep.frame$self.size[i] <- sum(right.age & right.site)
      rep.frame$sample.size[i] <- length(calib.vals)
      
      cat(paste(as.character(rep.frame[i, 1]), 
                        rep.frame[i,2], 
                        round(i/nrow(rep.frame), 4)*100, sep=', '), '\n')
#     }
  }
}

save(rep.frame, file='data/rep.frame.RData')
load('data/rep.frame.RData')
