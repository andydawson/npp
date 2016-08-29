#!/usr/bin/Rscript
library(plotrix)

source("config20_v1")
vers = 'v2'

setwd(projectDir)
centroids <- read.csv('centroids.csv')
nPlots <- nrow(centroids)

ftPerMeter <- 3.2808399

lastYear <- 2012
firstYear <- 1960#1905
years <- firstYear:lastYear
nT <- length(years)

growDays <- 30+31+30+31+31+30+31

library(dplR)
library(fields)
library(reshape2)

# setwd(dataDir)
# setwd('lyford')

wd <- paste0(dataDir, '/lyford')
wd20 <- paste0(dataDir, '/lyford20')


# encapsulate firstYear in getTimeIndex closure
getTimeIndex_gen <- function(firstYear) {
    function(year) {
        year - firstYear + 1
    }
}
getTimeIndex <- getTimeIndex_gen(firstYear)


## read data in

censusFull <- read.csv(paste0(wd, '/', censusFile), stringsAsFactors = FALSE)
ringMeta <- read.csv(paste0(wd20, '/', ringMetaFile), stringsAsFactors = FALSE)
treeMeta <- read.csv(paste0(wd20, '/', treeMetaFile), skip = 2, stringsAsFactors = FALSE)
censusDates <- read.csv(paste0(wd, '/', censusSampleDatesFile), stringsAsFactors = FALSE)
rwSampleDates <- read.csv(paste0(wd, '/', rwPlotSampleDatesFile), stringsAsFactors = FALSE)

# # read in plot data to get live dead status
# rwStatus <- list(length=nPlots)
# for (j in 1:nPlots){
#   rwStatus[[j]] <- read.csv(paste0(wd20, '/', 'Lyford_Data_Full20m_Final/LyFordPlot', j, '.csv'), 
#                             skip=7,header=TRUE,stringsAsFactors = FALSE)
# }
    
# remove second tree that has Audrey id of 1863; this one is not consistent with Audrey's data
treeMeta <- treeMeta[!(treeMeta$Site == "LF2" & treeMeta$Tree.Number == 27), ]
ringMeta <- ringMeta[!(ringMeta$SITE == "LF2" & ringMeta$TREE == 27), ]

treeMeta <- treeMeta[!(treeMeta$Site == "LF2" & treeMeta$Tree.Number == 45), ]

rwFiles <- list.files(paste0(wd20, '/', "RW"))
rwFiles <- rwFiles[grep(".rwl$", rwFiles)]
rwData <- list()
for(fn in rwFiles) {
    id <- gsub(".rw", "", fn)
    rwData[[id]] <- t(read.tucson(file.path(wd20, "RW", fn)))  # rows are tree, cols are times
}

## initial processing of input datasets
names(censusFull)[names(censusFull) == "treeid"] <- "census_id"

ringMeta$SITE <- as.numeric(gsub("LF", "", ringMeta$SITE))
ringMeta$X <-  NULL
names(ringMeta) <- tolower(names(ringMeta))

names(treeMeta) <- tolower(names(treeMeta))
treeMeta$site <- as.numeric(gsub("LF", "", treeMeta$site))
wh <- which(names(treeMeta) == "tag")
names(treeMeta)[wh] <- "census_id" # to match census file
treeMeta$id <- treeMeta$site*100 + treeMeta$tree.number


idx_prob = c()
treeMeta$dbh_year = NA
lastTreeIdFirstSampling <- c(21, 31, 37)
for (i in 1:nrow(treeMeta)){
  if (treeMeta$status[i] == 'Li'){
  if (treeMeta$tree.number[i]<lastTreeIdFirstSampling[treeMeta$site[i]]){
    treeMeta$dbh_year[i] = 2013
  } else {
    treeMeta$dbh_year[i] = 2014
  }
  } else {
    if (any((ringMeta$site==treeMeta$site[i]) & 
                           (ringMeta$tree==treeMeta$tree.number[i]))){
    treeMeta$dbh_year[i] = max(ringMeta[which((ringMeta$site==treeMeta$site[i]) & 
                           (ringMeta$tree==treeMeta$tree.number[i])),'meas_outer'])
    } else {
      print(i)
      idx_prob = c(idx_prob, i)
    }
  }
}


## initial processing of census sample dates
threshold <- function(vals, minDate, lower = TRUE) {
    if(!is(minDate, "Date")) minDate <- as.Date(minDate)
    if(length(minDate) == 1) minDate <- rep(minDate, length(vals))
    if(lower) {
        vals[julian(vals) < julian(minDate)] <- minDate[julian(vals) < julian(minDate)]
    } else vals[julian(vals) > julian(minDate)] <- minDate[julian(vals) > julian(minDate)]
    return(vals)
}

inds <- which(censusDates$X1987.1992_start == "")
censusDates$X1987.1992_start[inds] <-  censusDates$X1987.1992_end[inds]
censusDates$yr1991 <- as.numeric(gsub(".*(19[89][0-9])$", "\\1", censusDates$X1987.1992_start))

start <- as.Date("1969-3-31")
end <- as.Date("1969-11-1")
tmp1 <- as.Date(censusDates$X1969_start, format = "%m/%d/%Y")
tmp1 <- threshold(tmp1, start)
tmp1 <- threshold(tmp1, end, lower = FALSE)
tmp2 <- as.Date(censusDates$X1969_end, format = "%m/%d/%Y")
tmp2 <- threshold(tmp2, start)
tmp2 <- threshold(tmp2, end, lower = FALSE)
censusDates$day69 <- (julian(tmp1) + julian(tmp2)) / 2 - julian(start)

start <- as.Date("1975-3-31")
end <- as.Date("1975-11-1")
censusDates$X1975[is.na(censusDates$X1975)] <- "5/15/1975" # rough average of other blocks
tmp1 <- as.Date(censusDates$X1975, format = "%m/%d/%Y")
tmp1 <- threshold(tmp1, start)
tmp1 <- threshold(tmp1, end, lower = FALSE)
censusDates$day75 <- julian(tmp1) - julian(start)

start <- as.Date(paste0(censusDates$yr1991, "-3-31"))
end <- as.Date(paste0(censusDates$yr1991, "-11-1"))
tmp1 <- as.Date(censusDates$X1987.1992_start, format = "%m/%d/%Y")
tmp1 <- threshold(tmp1, start)
tmp1 <- threshold(tmp1, end, lower = FALSE)
tmp2 <- as.Date(censusDates$X1987.1992_end, format = "%m/%d/%Y")
tmp2 <- threshold(tmp2, start)
tmp2 <- threshold(tmp2, end, lower = FALSE)
censusDates$day91 <- (julian(tmp1) + julian(tmp2)) / 2 - julian(start)

start <- as.Date("2001-3-31")
end <- as.Date("2001-11-1")
tmp1 <- as.Date(censusDates$X2001_start, format = "%m/%d/%Y")
tmp1 <- threshold(tmp1, start)
tmp1 <- threshold(tmp1, end, lower = FALSE)
tmp2 <- as.Date(censusDates$X2001_end, format = "%m/%d/%Y")
tmp2 <- threshold(tmp2, start)
tmp2 <- threshold(tmp2, end, lower = FALSE)
censusDates$day01 <- (julian(tmp1) + julian(tmp2)) / 2 - julian(start)

start <- as.Date("2011-3-31")
end <- as.Date("2011-11-1")
tmp1 <- as.Date(censusDates$X2011_start, format = "%m/%d/%Y")
tmp1 <- threshold(tmp1, start)
tmp1 <- threshold(tmp1, end, lower = FALSE)
tmp2 <- as.Date(censusDates$X2011_end, format = "%m/%d/%Y")
tmp2 <- threshold(tmp2, start)
tmp2 <- threshold(tmp2, end, lower = FALSE)
censusDates$day11 <- (julian(tmp1) + julian(tmp2)) / 2 - julian(start)

start <- as.Date("1962-3-31")
end <- as.Date("1962-11-1")
tmp1 <- as.Date(censusDates$X1962, format = "%m/%d/%Y")
censusDates$day62 <- julian(tmp1) - julian(start)




## get census trees only in plots
census <- list(); length(census) <- nPlots
for(i in seq_len(nPlots)) {
    dist <- rdist(centroids[i, 2:3], censusFull[ , c('xsite', 'ysite')])
    census[[i]] <- censusFull[dist < plotRadius*ftPerMeter, ]
    census[[i]]$site <- i
    census[[i]]$dist_census <- dist[dist < plotRadius*ftPerMeter]
}

# check that we are getting the right trees
par(mfrow=c(1,1))
plot(censusFull$xsite, censusFull$ysite, asp=1, col="lightgrey")#, xlim=c(-500,500))
points(centroids$x, centroids$y, col="red", pch='X')

for (i in 1:nPlots){
  points(census[[i]]$xsite, census[[i]]$ysite, pch=19, col='blue')
  draw.circle(centroids$x[i], centroids$y[i], plotRadius*ftPerMeter)
}
# looks right to me!


census <- do.call(rbind, census)
dbhCols <- grep("dbh", names(census))
condCols <- grep("cond", names(census))

existFun <- function(x) {
    dbhs <- x[ , dbhCols]
    conds <- x[ , condCols]
    numAlive <- sum(conds == "L", na.rm = TRUE)
    numAlive[is.na(numAlive)] <- 0
    numAlive[]
    return(numAlive > 0)
}

# remove trees that we deem as out
out_trees <- c(3829, 1866, 1890, 2111, 995, 1089)
census = census[which(!(census$census_id %in% out_trees)),]

numAlive <- apply(census, 1, function(x) sum(x[condCols] == "L", na.rm = TRUE))
noDbh <- apply(census[ , dbhCols], 1, function(x) sum(is.na(x)) == length(dbhCols))
numAlive[noDbh] <- 0
census <- census[numAlive > 0, ]
# for now throw out trees never alive in a census - assume these were never big so limited impact on biomass increment

censusYears <- strsplit(censusYears, ",")[[1]]
censusYears2digit <- substring(censusYears, 3, 4)
nCensus <- length(censusYears)


census <- merge(census, treeMeta, by=c('census_id', 'site'), all.x = TRUE, all.y = FALSE)

census <- merge(census, censusDates[ , c('Block', 'yr1991', paste0('day', censusYears2digit))], by.x = 'block', by.y = 'Block', all.x = TRUE, all.y = FALSE)

census <- census[order(census$site, census$id), ]

census$stat_id <- 1:nrow(census)

# points(census[which(census$census_id==1890),'xsite'], census[which(census$census_id==1890),'ysite'], pch=19, col='red')
# points(census[which(census$census_id==2111),'xsite'], census[which(census$census_id==2111),'ysite'], pch=19, col='red')

## extract dbh values for living trees

dbhCols <- paste0('dbh', censusYears2digit)
condCols <- paste0('cond', censusYears2digit)
dayCols <- paste0('day', censusYears2digit)

colnames(census)[colnames(census) == "species.x"] = "species"

# get the dbh values from the RW plots
cols <- c("census_id", "site", "id", "stat_id", "species", "dbh", "dbh_year")
tmp <- census[ , cols]
dbhRW <- melt(tmp, id.vars = c("species","census_id", "site", "stat_id", "id","dbh_year"))
dbhRW <- dbhRW[!is.na(dbhRW$value),]
dbhRW$yr = substring(dbhRW$dbh_year,3,4)
dbhRW$day  = rwSampleDates[match(dbhRW$dbh_year, rwSampleDates$year),'date']

# FIXME
for (i in 1:nrow(dbhRW)){
  if (!is.na(dbhRW$dbh_year[i])){
    print(i)
    start <- as.Date(paste0(dbhRW$dbh_year[i], "-3-31"))
    end   <- as.Date(paste0(dbhRW$dbh_year[i], "-11-1"))
    if (dbhRW$dbh_year[i] %in% c(2013,2014)) {
      tmp1 <- as.Date(dbhRW$day[i])
      tmp1 <- threshold(tmp1, start)
      tmp1 <- threshold(tmp1, end, lower = FALSE)
    } else{
      tmp1   <- as.Date(paste0(dbhRW$dbh_year[i], "-11-1"))
    }
      dbhRW$day[i] <- julian(tmp1) - julian(start)
  } else {
    next
  }
}

dbhRW = dbhRW[!is.na(dbhRW$day),]

# get the dbh values from the census
cols <- c("census_id", "site", "id", "stat_id", "species", dbhCols, condCols)
tmp <- census[ , cols]
dbh <- melt(tmp, id.vars = c("species","census_id", "site", "stat_id", "id", condCols))
dbh$yr=substring(dbh$variable,4,5)

cols <- c("census_id", "site", "id", "stat_id", "species", dayCols, condCols)
tmp <- census[ , cols]
day <- melt(tmp, id.vars = c("species","census_id", "site", "stat_id", "id", condCols))
day$yr=substring(day$variable,4,5)

names(day)[names(day) == 'value'] <- 'day'

dbh <- merge(dbh, day[ , c('stat_id', 'day', 'yr')], all.x = TRUE, all.y = FALSE)
# get in terms of 4-digit year, not 2-digit
match <- data.frame(yr = censusYears2digit, year = as.numeric(censusYears)) 
dbh <- merge(dbh, match, all.x = TRUE, all.y = FALSE)

colMatch = data.frame(censusYears2digit, which(names(dbh) %in% condCols))
names(colMatch) = c('yr','col')
dbh2 = merge(dbh, colMatch, by.x= 'yr', by.y = 'yr')
dbh2$status <- dbh2[cbind(1:nrow(dbh2), dbh2$col)]
dbh2 <- merge(dbh2, census[ , c('stat_id', 'yr1991')], all.x = TRUE, all.y = FALSE)
dbh2$year[dbh2$year == 1991] <- dbh2$yr1991[dbh2$year == 1991]

dbh <- subset(dbh2, status == "L", c('yr', 'year', 'species', 'day', 'census_id', 'id', 'stat_id', 'site', 'value'))
dbh <- subset(dbh, !is.na(value))


dbh = rbind(data.frame(dbh, type=rep('census')), 
            data.frame(yr=dbhRW$yr,
                       year=dbhRW$dbh_year, 
                       species=dbhRW$species,
                       day=dbhRW$day,
                       census_id=dbhRW$census_id,
                       id=dbhRW$id,
                       stat_id=dbhRW$stat_id,
                       site=dbhRW$site,
                       value=dbhRW$value,
                       type=rep('rw')))#[,('yr', 'year','species', 'day', 'census_id', 'id', 'stat_id','site','value')])


dbh_sub = dbh[which(dbh$year>=2011),]

dbh_sub = dcast(dbh_sub, stat_id~type)
dbh_sub = transform(dbh_sub, census=as.numeric(census), rw=as.numeric(rw))


dbh_sub[which(abs(dbh_sub$census - dbh_sub$rw) > 8),'stat_id']
dbh[which((dbh$stat_id %in% dbh_sub[which(abs(dbh_sub$census - dbh_sub$rw) > 8),'stat_id']) & (dbh$year >= 2011)), ]


ggplot(dbh_sub) + geom_point(aes(x=census, y=rw))
ggsave(file='figures/dbh_census_vs_rw.pdf')


dbh_sub$diff = dbh_sub$rw - dbh_sub$census
ggplot(dbh_sub) + geom_point(aes(x=census, y=diff)) + xlab('Census DBH') + ylab('RW - Census DBH')
ggsave(file='figures/dbh_diff_vs_census.pdf')

## dbh is a core data object going into the stat model

saplings <- subset(dbh2, (is.na(status) & (yr != 62)))
saplings <- rbind(saplings, subset(dbh2, (is.na(status) & (yr == 62) & (site == 1))))
## saplings is a core data object going into the stat model

# sum over multiple cores to get tree-level increment
# insert 0 for NA when match non-NA in another core

combineCores <- function(rw) {
    zeroGrowthFlag <- -1
    id <- as.numeric(substring(dimnames(rw)[[1]], 3, 5))
    rw <- rw[!is.na(id), ]
    id <- id[!is.na(id)]
    orient <- substring(dimnames(rw)[[1]], 6, 7)
    trees <- unique(id)
    plot <- unique(as.numeric(substring(dimnames(rw)[[1]], 3, 3)))
    if(length(plot) > 1) stop("multiple plots in object")
    vals <- matrix(zeroGrowthFlag, length(trees), ncol(rw)) # -1 is code for no ring 
    dimnames(vals)[[2]] <- dimnames(rw)[[2]]
    for(i in seq_along(trees)) {
        tree_rw <- rw[id == trees[i], , drop = FALSE]
        cat("processing ", plot, trees[i], substring(dimnames(tree_rw)[[1]], 6, 6), "\n")
        anyNonNA <- apply(tree_rw, 2, function(x) sum(!is.na(x)) > 0)
        # vals[i, anyNonNA] <- (colSums(tree_rw, na.rm = TRUE) / (0.5*nrow(tree_rw)))[anyNonNA]
        vals[i, anyNonNA] <- (colSums(tree_rw, na.rm = TRUE) / (nrow(tree_rw)))[anyNonNA]
        if(!anyNonNA[1]) {
            lastZero <- max(which(cumsum(vals[i,]) == -(1:ncol(rw))))
            vals[i, 1:lastZero] <- NA  # NA is for before tree existed
        }
    }
    dimnames(vals)[[1]] <- trees
    return(vals)
}
# 
# # organize cores in the same way as combineCores, without cobining cores for individuals trees
# orgCores <- function(rw) {
#   zeroGrowthFlag <- -1
#   id <- as.numeric(substring(dimnames(rw)[[1]], 3, 5))
#   rw <- rw[!is.na(id), ]
#   id <- id[!is.na(id)]
#   orient <- substring(dimnames(rw)[[1]], 6, 7)
#   trees <- unique(id)
#   plot <- unique(as.numeric(substring(dimnames(rw)[[1]], 3, 3)))
#   if(length(plot) > 1) stop("multiple plots in object")
#   vals <- matrix(zeroGrowthFlag, length(trees), ncol(rw)) # -1 is code for no ring 
#   dimnames(vals)[[2]] <- dimnames(rw)[[2]]
#   for(i in seq_along(trees)) {
#     tree_rw <- rw[id == trees[i], , drop = FALSE]
#     cat("processing ", plot, trees[i], substring(dimnames(tree_rw)[[1]], 6, 6), "\n")
#     anyNonNA <- apply(tree_rw, 2, function(x) sum(!is.na(x)) > 0)
#     vals[i, anyNonNA] <- (colSums(tree_rw, na.rm = TRUE) / (0.5*nrow(tree_rw)))[anyNonNA]
#     if(!anyNonNA[1]) {
#       lastZero <- max(which(cumsum(vals[i,]) == -(1:ncol(rw))))
#       vals[i, 1:lastZero] <- NA  # NA is for before tree existed
#     }
#   }
#   dimnames(vals)[[1]] <- trees
#   return(vals)
# }


# organize cores in the same way as combineCores, without cobining cores for individuals trees
orgCores <- function(rw) {
  zeroGrowthFlag <- -1
#   id <- as.numeric(substring(dimnames(rw)[[1]], 3, 5))
  id <- substring(dimnames(rw)[[1]], 3, 6)
  rw <- rw[!is.na(id), ]
  id <- id[!is.na(id)]
  orient <- substring(dimnames(rw)[[1]], 6, 7)
  trees <- unique(id)
  plot <- unique(as.numeric(substring(dimnames(rw)[[1]], 3, 3)))
  if(length(plot) > 1) stop("multiple plots in object")
  vals <- matrix(zeroGrowthFlag, length(id), ncol(rw)) # -1 is code for no ring 
  dimnames(vals)[[2]] <- dimnames(rw)[[2]]
  for(i in seq_along(trees)) {
    tree_rw <- rw[id == trees[i], , drop = FALSE]
    cat("processing ", plot, trees[i], substring(dimnames(tree_rw)[[1]], 6, 6), "\n")
    anyNonNA <- apply(tree_rw, 2, function(x) sum(!is.na(x)) > 0)
    # vals[i, anyNonNA] <- (colSums(tree_rw, na.rm = TRUE) / (0.5*nrow(tree_rw)))[anyNonNA]
    vals[i, anyNonNA] <- (colSums(tree_rw, na.rm = TRUE) / (nrow(tree_rw)))[anyNonNA]
    if(!anyNonNA[1]) {
      lastZero <- max(which(cumsum(vals[i,]) == -(1:ncol(rw))))
      vals[i, 1:lastZero] <- NA  # NA is for before tree existed
    }
  }
  dimnames(vals)[[1]] <- trees
#   return(vals)
  # vals = rw
  dimnames(vals)[[2]] <- dimnames(rw)[[2]]
  dimnames(vals)[[1]] <- id

  return(vals)
}

# # single measurement per tree (average over cores)
# treeRW <- lapply(rwData, combineCores)

# multiple measurements per tree (no averaging)
treeRW <- lapply(rwData, orgCores)

# combine rw data to get increment dataset
ringYears <- sapply(treeRW, function(x) as.numeric(dimnames(x)[[2]]))
nTreesWithCores <- sum(sapply(treeRW, function(x) nrow(x)))
incr <- matrix(NA, nTreesWithCores, nT)
dimnames(incr)[[2]] <- years

start <- 1
for(i in seq_along(treeRW)) {
    end <-  start -1 + nrow(treeRW[[i]])
    include <- which(ringYears[[i]] <= lastYear & ringYears[[i]] >= firstYear)
    incr[start:end, as.character(ringYears[[i]][include])] <- treeRW[[i]][ , include]
 
    start <- end + 1
}
dimnames(incr)[[1]] <- unlist(sapply(treeRW, function(x) dimnames(x)[[1]]))

# need to figure out lastDate and put in census
census$lastYear <- 2012
census$lastYear[census$cond11 != "L" & !is.na(census$cond11)] <- 2010
census$lastYear[census$cond01 != "L" & !is.na(census$cond01)] <- 2000
census$lastYear[census$cond91 != "L" & !is.na(census$cond91)] <- census$yr1991[census$cond91 != "L" & !is.na(census$cond91)] - 1
census$lastYear[census$cond75 != "L" & !is.na(census$cond75)] <- 1974
census$lastYear[census$cond69 != "L" & !is.na(census$cond69)] <- 1968
census$lastYear[census$cond62 != "L" & !is.na(census$cond62)] <- 1962

# now walk through trees with rings
firstNoGrowth <- function(x) {
    tmp <- which(x == -1)
    if(length(tmp)) return(as.numeric(names(x)[min(tmp)])) else return(NA)
}
zeroRing <- apply(incr, 1, firstNoGrowth)
zeroRing[is.na(zeroRing)] <- lastYear + 1 # kludge so that 2012 is given as last year in code below for trees with ring in 2012

censusMatches <-  which(census$id %in% substr(names(zeroRing),1,3))

census$lastYear[censusMatches] <- zeroRing[as.character(census$id[censusMatches])] - 1


# temporary restrict to one site to quicken calcs
if(F){
census <- subset(census, site == 1)
dbh <- dbh[substring(dbh$id, 1,1)=="1", ]
wh <- substring(dimnames(incr)[[1]], 1 ,1) == "1"
incr <- incr[wh,]
}

tbl <- table(census$species)
taxaMatch <- data.frame(tbl); names(taxaMatch) <- c('species', 'count')
taxaMatch$taxon <- 1:nrow(taxaMatch)

# get the site number
site_data    <- data.frame(stat_id=dbh$stat_id, site=dbh$site)
site_data    <- site_data[!duplicated(dbh$stat_id),]
tree_site_id <- site_data[order(site_data$stat_id), 'site']
nSites       <- length(unique(tree_site_id))

# figure out plot number business!
plot_data <- dbh
plot_data <- plot_data[!duplicated(plot_data$stat_id),]


dbh$value <- as.numeric(dbh$value)

nDBH <- nrow(dbh)
logDobs <- log(dbh$value)
dbh_tree_id <- dbh$stat_id
dbh_day_id <- as.numeric(dbh$day) / growDays
dbh_year_id <- getTimeIndex(dbh$year)

nSaplings <- nrow(saplings)
sapling_tree_id <- saplings$stat_id
sapling_year_id <- getTimeIndex(saplings$year)
max_size <- ifelse(saplings$year == 1969, 4.5, 5)

incrMelted <- melt(incr)
names(incrMelted) <- c('id', 'year', 'incr')

incrMelted$orient = tolower(substr(incrMelted$id, 4, 4))

incrMelted$id     = substr(incrMelted$id, 1, 3)

incrData <- merge(incrMelted, census[ , c('id', 'stat_id')])
incrData <- subset(incrData, !is.na(incr))
incrData <- subset(incrData, incr != -1)

nWidths <- nrow(incrData)
incr_tree_id <- incrData$stat_id
incr_year_id <- getTimeIndex(incrData$year)
logXobs <- log(incrData$incr)
# fudge 0 incr for now
logXobs[logXobs < log(.02)] <- log(.02)

tmp <- merge(census[ , c('stat_id', 'species')], taxaMatch[ , c('species', 'taxon')], all.x = TRUE, all.y = FALSE)
tmp <- tmp[order(tmp$stat_id), ]
taxon <- tmp$taxon
nTaxa <- nrow(taxaMatch)

n <- nrow(census)
dead <- census$lastYear < lastYear

last_time <- getTimeIndex(census$lastYear)

last_time = vector(length=n)
for (i in 1:n){
  last_time[i] = max(c(dbh$year[which(dbh$stat_id == i)], incrData$year[which(incrData$stat_id == i)]))
  print(last_time[i])
}

last_ti = match(years, last_time)

# save(n, nT, nDBH, nWidths, nSites, dbh_tree_id, dbh_day_id, dbh_year_id, incr_tree_id, incr_year_id, tree_site_id,
#      logDobs, logXobs, last_time, census, dbh, incrData, incr, treeMeta, ringMeta,
#      nSaplings, sapling_tree_id, sapling_year_id, max_size, taxon, nTaxa,
#      file = 'nimble/test/data.Rda')
save(incrData, incr, treeMeta, ringMeta,
     file = paste0('/home/adawson/Documents/projects/npp/data/meas/incr_data_', vers,'.Rda'))

######################################################################################################################################
## make nimble data
######################################################################################################################################

incrDataOrig = incrData

year_start = 1960
year_end   = max(incrData$year)

# reset years for now
years = seq(year_start, year_end)
incrData  = incrData[which(incrData$year %in% years),]
trees_inc = sort(unique(incrData$stat_id))

# order by tree and year
incrData = incrData[order(incrData$stat_id, incrData$year),]

# how many measurements for each tree for each year?
ncores = vector(length=nrow(incrData))
for (i in 1:length(trees_inc)){
  for (t in 1:length(years)){
    idx = which((incrData$stat_id == trees_inc[i]) & (incrData$year == years[t]))
    ncores[idx] = rep(length(idx), length(idx))
  }
}

incrData = cbind(incrData, ncores)
incrData = data.frame(incrData, measno = seq(1, nrow(incrData)))

head(incrData)

N_inc   = nrow(incrData) # number of measurement 
m2t     = incrData$year
m2tree  = incrData$stat_id
m2treecode = incrData$id
m2nc = incrData$ncores
m2ti    = match(m2t, years)
m2orient = incrData$orient
Xobs    = incrData$incr
Xobs[Xobs==0] = 0.0001
logXobs = log(Xobs)

dbh = dbh[which(dbh$year %in% years),]
dbh = dbh[order(dbh$stat_id, dbh$year),]
N_dbh   = nrow(dbh)
logDobs = log(dbh$value)
dbh_tree_id = dbh$stat_id
dbh_day_id  = dbh$day / growDays
dbh_year_id = dbh$year - year_start + 1#getTimeIndex(dbh$year)
dbh_tree_code = dbh$id

ids_table = dbh[,c(5,6,7)]
ids_table = ids_table[!duplicated(ids_table$stat_id),]

trees   = sort(unique(dbh$stat_id))
N_trees = length(unique(dbh$stat_id))
N_years = length(years)

# get the site number
site_data    <- data.frame(stat_id=dbh$stat_id, site=dbh$site, census_id=dbh$census_id)
site_data    <- site_data[!duplicated(dbh$stat_id),]
tree_site_id <- site_data[order(site_data$stat_id), 'site']
tree_census_id <- site_data[order(site_data$stat_id), 'census_id']
N_sites      <- length(unique(tree_site_id))


census_years = c(1969, 1975, 1991, 2001, 2011)

last_time_data = vector(length=N_trees)
last_time = vector(length=N_trees)
for (i in 1:N_trees){
  tree = trees[i]
  print(tree)
  last_time[i] = max(c(dbh$year[which(dbh$stat_id == tree)], incrData$year[which(incrData$stat_id == tree)]), na.rm=TRUE)
  last_time_data[i] = last_time[i]
  if (last_time[i] %in% census_years) {
    if (which(census_years == last_time[i]) == length(census_years)) {
      last_time[i] = max(years)
    } else if (last_time_data[i] == 1991) {
      last_time[i] = 2001
    } else {
      last_time[i] = census_years[which(census_years == last_time[i]) + 1]
    }
  } else if (last_time_data[i] == 1992) {
    last_time[i] = 2001
  }
  print(last_time[i])
}
last_time = as.numeric(last_time)
last_time_data = as.numeric(last_time_data)

X_ord = data.frame(meas=numeric(0), tree_id=numeric(0), year=numeric(0))
n = 1
for (i in 1:N_trees){
  print(i)
  print(last_time[i])
  year = seq(year_start, last_time[i])
  meas = seq(n, n+length(year)-1)
  n = n + length(year)
  
  X_ord = rbind(X_ord, data.frame(meas=meas, tree_id=rep(trees[i], length(year)), year=year))
}

x2tree  = X_ord$tree_id
x2year  = match(X_ord$year, years) 
last_ti = last_time-year_start +1

# X_ord = aggregate(orient~year + stat_id, incrData, function(x) length(unique(x)))
# X_ord = X_ord[order(X_ord$stat_id, X_ord$year),]
N_X   = nrow(X_ord)
N_D   = N_X

meas2x = vector(length=N_inc)
for (i in 1:N_inc) {
  stat_id = incrData$stat_id[i]
  year    = incrData$year[i]
  
  meas2x[i] = which((X_ord$tree_id == stat_id) & (X_ord$year == year))
}

meas2d = vector(length=N_dbh)
for (i in 1:N_dbh) {
  stat_id = dbh$stat_id[i]
  year    = dbh$year[i]
  
  meas2d[i] = which((X_ord$tree_id == stat_id) & (X_ord$year == year))
}


# ncores_years = ncores$year
# ncores_tree  = ncores$stat_id
# ncores = ncores$orient
N_ncores = length(ncores)

i1core2m = vector(mode="integer")
i2core2m = vector(mode="integer")
i3core2m = vector(mode="integer")
i4core2m = vector(mode="integer")
n = 1
while (n <= length(m2nc)) {
  if (m2nc[n] == 1) {
    i1core2m = c(i1core2m, n)
    n = n + 1
  } else if (m2nc[n] == 2) {
    i2core2m = c(i2core2m, n)
    n = n + 2
  } else if (m2nc[n] == 3) {
    i3core2m = c(i3core2m, n)
    n = n + 3
  } else if (m2nc[n] == 4) {
    i4core2m = c(i4core2m, n)
    n = n + 4
  } else {
    print("WTF")
  }
}

n1cores = length(i1core2m)
n2cores = length(i2core2m)
n3cores = length(i3core2m)
n4cores = length(i4core2m)

ones = rep(1, 4)
open_dbh =25
N_taxa=nTaxa

cs_last_ti = cumsum(last_ti)
first_ti   = c(1,cs_last_ti[1:(length(cs_last_ti)-1)]+1)

N_saplings = nSaplings 
sap2x = vector(length=N_saplings)
for (i in 1:N_saplings) {
  sap2x[i] = which((X_ord$tree_id == sapling_tree_id[i]) & (X_ord$year == years[sapling_year_id[i]]))
}

not_sap2x = which(!(seq(1,N_X) %in% sap2x))

X = rep(0.1, N_X)
D = rep(0.1, N_X)
logX = rep(log(0.1), N_X)
D0 = rep(3, N_trees)

beta = rep(0.1, N_trees)
beta_t = rep(0.1, N_years)
beta0 = 0.1
sig_x_obs = 0.7
sig_d_obs = 0.01
sig_d = 0.1
sig_d_sap = 0.1
sig_x = 0.1
beta_sd = 0.1
beta_t_sd = 0.1
beta_spp_sd = 0.1
beta_spp = rep(0.1, N_taxa)
beta_slope = 0.1

b0 = 0.1
b1 = 10

tau2 = 0
tau3 = 0
tau4 = 0

nu = 0.1
rho = 0.1

D_saplings = rep(3, N_saplings)
D_pre = rep(3, N_X-N_saplings)

dump(c('X', 'logX', 'D0', 'beta', 'beta_t', 'beta0', 'sig_x_obs', 'sig_d_obs', 'sig_d', 'sig_d_sap',
       'sig_x', 'beta_sd', 'beta_t_sd', 'beta_spp', 'beta_spp_sd', 'beta_slope', 
       'tau2', 'tau3', 'tau4', 'b0', 'b1', 'D_saplings', 'D_pre', 'nu', 'rho'), 
     file=paste0('data/dump/tree_full_20_', vers, '_inits.dump'))


# average rw values
incrData$incr[incrData$incr == 0] = 0.0001
if (any(incrData$incr==0)){
  print('Zero ring-widths!')
}

logOFmean=FALSE

incrData = data.frame(incrData, incr_log = log(incrData$incr))
if (logOFmean){
  incrAvg = aggregate(incr ~ stat_id + year, incrData, mean)
} else {
  incrAvg = aggregate(incr_log ~ stat_id + year, incrData, mean)
}
# incrAvg = aggregate(incr_log ~ stat_id + year, incrData, mean)
incrAvg = incrAvg[order(incrAvg$stat_id, incrAvg$year),]
incrAvg = data.frame(incrAvg, measno = seq(1, nrow(incrAvg)))

if (any(is.infinite(incrData$incr_log))){
  print('Infinitely negative log ring-widths!')
}

N_inc_a   = nrow(incrAvg) # number of measurement 
m2t_a     = incrAvg$year
m2tree_a  = incrAvg$stat_id
m2treecode_a = incrAvg$id
m2ti_a    = match(m2t_a, years)
# Xobs_a    = incrAvg$incr
# Xobs_a[Xobs_a==0] = 0.0001
# logXobs_a = log(Xobs_a)
if (logOFmean){
  logXobs_a    = log(incrAvg$incr)
} else {
  logXobs_a    = incrAvg$incr_log
}
  
# X_ord the same; still estimate same times and trees!
meas2x_a = vector(length=N_inc_a)
for (i in 1:N_inc_a) {
  stat_id = incrAvg$stat_id[i]
  year    = incrAvg$year[i]
  
  meas2x_a[i] = which((X_ord$tree_id == stat_id) & (X_ord$year == year))
}

meas2d = vector(length=N_dbh)
for (i in 1:N_dbh) {
  stat_id = dbh$stat_id[i]
  year    = dbh$year[i]
  
  meas2d[i] = which((X_ord$tree_id == stat_id) & (X_ord$year == year))
}

# census dates
aggregate(year~site, dbh, function(x) sort(unique(x)))


dump(c('N_inc', 'N_dbh', 'N_X', 'N_D', 'N_trees', 'N_years', 'N_sites',
       'logXobs', 'm2t', 'm2ti', 'm2nc', 'm2tree', 'ncores',
       'N_inc_a', 'logXobs_a', 'm2t_a', 'm2ti_a', 'm2tree_a', 
       'logDobs', 'dbh_tree_id', 'dbh_day_id', 'dbh_year_id',
       'tree_site_id', 'tree_census_id',
       'n1cores', 'n2cores', 'n3cores', 'n4cores', 
       'i1core2m', 'i2core2m', 'i3core2m', 'i4core2m',
       'meas2x',  'x2year', 'x2tree', 
       'meas2d',
       'last_ti',
       'ones',
       'year_start', 'year_end',
       'taxon', 'N_taxa', 'open_dbh',
       'N_saplings', 'sap2x', 'not_sap2x','max_size', 'sapling_tree_id', 'sapling_year_id',
       'first_ti', 'cs_last_ti'),
     file=paste0('data/dump/tree_full_20_', vers, '.dump'))

save(N_inc, N_dbh, N_X, N_D, N_trees, N_years, N_sites,
     logXobs, m2t, m2ti, m2nc, m2treecode, m2tree, ncores,
     N_inc_a, logXobs_a, m2t_a, m2ti_a, m2tree_a, 
     logDobs, dbh_tree_id, dbh_day_id, dbh_year_id,
     tree_site_id, tree_census_id,
     n1cores, n2cores, n3cores, n4cores, 
     i1core2m, i2core2m, i3core2m, i4core2m,
     meas2x, x2year, x2tree, 
     meas2d,
     last_ti, last_time, last_time_data,
     ones,
     year_start, year_end,
     trees, years, m2orient,
     taxon, N_taxa, taxaMatch, open_dbh,
     N_saplings, sap2x, not_sap2x, max_size, sapling_tree_id, sapling_year_id,
     first_ti, cs_last_ti, dbh,
     file=paste0('data/dump/tree_full_20_', vers, '.rdata'))



census0=readRDS(file='r/census0.RDS')

foo0 = census0[census0$site==1,]

foo1 = census[census$site==1,]

miss = !(foo1$census_id %in% foo0$census_id)

miss1 = foo1[miss,]

plot(miss1$xsite, miss1$ysite)

####################################################################################################################################
## plot dead tree rws
####################################################################################################################################

dead_tree_ids = census$stat_id[census$census_id %in% c(1211, 4797, 4793,3844)]
N_rwdead = length(dead_tree_ids)

idx_dead = which(m2tree %in% dead_tree_ids)
years_dead = m2t[idx_dead]
rw_dead = exp(logXobs[m2tree %in% dead_tree_ids])

for (i in 1:N_rwdead){
  tree = dead_tree_ids[i]
  idx_dead = which(m2tree == tree)
  if (length(idx_dead) == 0 ) {
    print(paste0('No ring-widths for tree ', tree))
    next
  }
  
  years_dead = m2t[idx_dead]
  rw_dead = exp(logXobs[idx_dead])
  
  core_nums = unique(m2orient[idx_dead])
  
  par(mfrow=c(1,1))
  plot(c(0,0), xlim=c(min(years),max(years)), ylim=c(min(rw_dead), max(rw_dead)), xlab='Year', ylab='ring-width',
       type='n', main=paste0('Tree ', tree))
  
  # lines(years_dead, , col='blue')
  
  for (core in core_nums){
    idx = which((m2tree == tree) & (m2orient == core))
    yrs = years[m2ti[idx] ]
    lines(yrs, exp(logXobs[idx]), col='black', lty=2)
  }

}


####################################################################################################################################
## plot all rws
####################################################################################################################################

dead_tree_ids = census$stat_id[census$census_id %in% c(1211, 4797, 4793,3844)]
live_tree_ids = unique(m2tree[!(m2tree %in% dead_tree_ids)])
N_live = length(live_tree_ids)
N_dead = length(dead_tree_ids)

pdf(file='figures/live_vs_dead_rws.pdf')
par(mfrow=c(2,1))
plot(c(0,0), xlim=c(min(years),max(years)), ylim=c(min(exp(logXobs)), max(exp(logXobs))), xlab='Year', ylab='ring-width',
     type='n', main='Live tree ring-widths')
for (i in 1:N_live){
  
  # tree = trees[i]
  tree = live_tree_ids[i]
  
  if (!(any(m2tree_a == tree))) {
    print(paste0('No ring-widths for tree ', tree))
    next
  }
  
  idx_tree = which(m2tree_a == tree)
  
  years_tree = m2t_a[idx_tree]
  rw_tree = exp(logXobs_a[idx_tree])
  
  # core_nums = unique(m2orient[idx_tree])
  
  # lines(years_dead, , col='blue')
  
  if (tree %in% dead_tree_ids) {
    col='blue'
    lty = 1
  } else {
    col = 'black'
    lty = 2
  }
  
  lines(years_tree, exp(logXobs_a[idx_tree]), col=col, lty=2)
  
#   for (core in core_nums){
#     idx = which((m2tree == tree) & (m2orient == core))
#     yrs = years[m2ti[idx] ]
#     lines(yrs, exp(logXobs[idx]), col=col, lty=2)
#   }
  
}

plot(c(0,0), xlim=c(min(years),max(years)), ylim=c(min(exp(logXobs)), max(exp(logXobs))), xlab='Year', ylab='ring-width',
     type='n', main='Dead tree ring-widths')
for (i in 1:N_dead){
  
  # tree = trees[i]
  tree = dead_tree_ids[i]
  
  if (!(any(m2tree_a == tree))) {
    print(paste0('No ring-widths for tree ', tree))
    next
  }
  
  idx_tree = which(m2tree_a == tree)
  
  years_tree = m2t_a[idx_tree]
  rw_tree = exp(logXobs_a[idx_tree])
  
  # core_nums = unique(m2orient[idx_tree])
  
  # lines(years_dead, , col='blue')
  
  if (tree %in% dead_tree_ids) {
    col='blue'
    lty = 1
  } else {
    col = 'black'
    lty = 2
  }
  
  lines(years_tree, exp(logXobs_a[idx_tree]), col=col, lty=2)
  
  #   for (core in core_nums){
  #     idx = which((m2tree == tree) & (m2orient == core))
  #     yrs = years[m2ti[idx] ]
  #     lines(yrs, exp(logXobs[idx]), col=col, lty=2)
  #   }
  
}
dev.off()

####################################################################################################################################
## check on last_ti
####################################################################################################################################

dead_census = data.frame(census_id=census$census_id, stat_id=census$stat_id, dead=rep(NA, N_trees))

dead_census$dead[which(census$cond69 == 'D')] = '69'
dead_census$dead[which(census$cond75 == 'D')] = '75'
dead_census$dead[which(census$cond91 == 'D')] = '91'
dead_census$dead[which(census$cond01 == 'D')] = '01'
dead_census$dead[which(census$cond11 == 'D')] = '11'

dead_census$dead[which(dead_census$census_id==3844)] = 'D'
dead_census$dead[which(dead_census$census_id==1863)] = 'D'
dead_census$dead[which(dead_census$census_id==4797)] = 'D'

# check that they died since last census date
idx_prob = dead_census[which(dead_census$dead == 'D'),]
last_ti[idx_prob$stat_id] # good

# check if we think others are dead
dead_meas = which(last_ti<52)
dead_obs  = which(!is.na(dead_census$dead))

dead_both = which((last_ti < 52) & !is.na(dead_census$dead))

dead_rw  = which((last_ti < 52) & is.na(dead_census$dead)) 
dead_one = which((!(last_ti < 52)) & (!is.na(dead_census$dead)) | ((last_ti < 52) & is.na(dead_census$dead)) )

dead_census[dead_rw,]

dead_open = data.frame(ids_table[which(ids_table$stat_id %in% dead_rw),], last_ti=years[last_ti[dead_rw]])

idx_check = 

# dump(c('N_inc_a', 'N_dbh', 'N_X', 'N_D', 'N_trees', 'N_years',
#        'logXobs_a', 'm2t_a', 'm2ti_a', 'm2tree_a', 
#        'logDobs', 'dbh_tree_id', 'dbh_day_id', 'dbh_year_id',
#        'last_ti', 'year_start', 'year_end',
#        'meas2x', 'x2tree', 'x2year', 'meas2d'),
#      file=paste0('data/dump/tree_full_avgrws.dump'))
# 
# save(N_inc_a, N_dbh, N_X, N_D, N_trees, N_years,
#      logXobs_a, m2t_a, m2ti_a, m2treecode_a, m2tree_a,
#      logDobs, dbh_tree_id, dbh_day_id, dbh_year_id,
#      last_ti, year_start, year_end,
#      trees, years,
#      meas2x, x2tree, x2year, meas2d,
#      N_inc, logXobs, m2t, m2ti, m2nc, m2treecode, m2tree, ncores,
#      file=paste0('data/dump/tree_full_avgrws.rdata'))
# 
# 
# 














# # need y2 and y3 ???=
# 
# incr[incr == 0] = 0.002
# 
# incr.na = incr
# 
# N_trees = length(unique(substr(rownames(incr),1,3)))
# T       = ncol(incr)
# N_tot   = nrow(incr)
# N_cores = as.vector(unname(table(substr(rownames(incr),1,3))))
# years   = as.numeric(colnames(incr))
# 
# # y = unname(incr.na.rm)
# y = unname(incr)
# logy = log(y)
# 
# dump(c('N_trees', 'T', 'N_tot', 'N_cores', 'y', 'logy'),
#      file=paste0('data/dump/rw_data', suff, '.dump'))
# 
# save(N_trees, T, N_tot, N_cores, y, logy,
#      years, incr.na,
#      file=paste0('data/dump/rw_data', suff, '.rdata'))
# 
# # figure out what to do with na
# na.rows = apply(y, 1, function(x) any(is.na(x)))
# rownames(incr)[na.rows]
# 
# # # what to do with zero rows
# # zero.rows = apply(y, 1, function(x) any(x==0))
# # rownames(incr.na.rm)[zero.rows]
# # incr.na.rm[zero.rows,]
# 
# rw = melt(incr.na)
# rw = cbind(rw[,1], rw)
# colnames(rw) = c('tree', 'core', 'year', 'rw')
# 
# rw$tree = substr(rw$tree,1,3)
# 
# rw$core = tolower(substr(rw$core, 4,4))
# angles = data.frame(core=c('n', 'e', 's', 'w', 'a', 'b', 'c'), theta=c(0, 90, 180, 270, 0, 90, 180))
# 
# rw$core[rw$core ]
# 
# rw$angle = angles$theta[match(rw$core, angles$core)]
# 
# trees = unique(rw$tree)
# 
# for (i in 1:length(trees)){
#  tree = trees[i]
#  rw_sub = rw[which(rw$tree == tree),]
#  
#  orients = unique(rw_sub$core)
#   c(dist(unique(rw_sub$angle)))
#  # rw_sub[which(rw_sub$core == 'n'), 'year']
#  
#  c(dist(v))
#  
#  
#  years = 
#  
#    
# }
# 
# 
# 
# cores = unique(rw$core)
# trees = unique(rw$tree)
# 
# #core_dict = list()
# #for (c in 1:length(cores)) {
# #  core_dict[[levels(cores)[c]]] = c
# #}
# 
# core_num = vector(length=length(cores))
# for (i in 1:length(trees)){
#   tree = trees[i]
#   idx = which(substr(cores,1,3) == tree)
#   core_num[idx] = seq(1, length(idx))
# }
# 
# rw$core = core_num[match(rw$core, cores)]
# 
# rw = rw[!is.na(rw$rw),]
# rw = data.frame(rw, measno = seq(1, nrow(rw)))
# 
# head(rw)