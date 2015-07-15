theta1 <- c(-4.583810004,0,1,0.1)
tw_groupin <- wdcMerged$tw_group[1]
wdcMergedday <- subset(wdcMerged, tw_group==tw_groupin)
deltain_tw <- rep(-3, nrow(wdcMergedday))
density_ridership_col <<- 3
density_metro_col <<- 4
