# script for testing SPDE model

# fit model ----
out = fitSPDEsimDat(wellTestDat, seismicTestDat)

names(out)

# plot predictions, seismic data, truth ----

quilt.plot()








