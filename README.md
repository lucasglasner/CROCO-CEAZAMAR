### CROCO-CEAZAMAR Forecast/hindcast system

This repository contains the scripts and preprocessing tools for running the croco model in forecast mode. The development is part of CEAZA efforts
to provide monitoring and forecast for the Coquimbo Region (Chile) coastal zone.

From a technical point of view, the methodology consist of a daily update of two croco simulations: (i) the croco hindcast and (ii) the croco forecast.
The hindcast is a simulation forced by the near real time (NRT) data of Mercator Analysis and ERA5 without any kind of data assimilation. On the other hand,
the forecast simulation maintains the ocean forcing data but changes the atmospheric forcing to GFS analysis and forecast dataset. The following scheme shows 
the data flow of the hindcast and forecast simulation:

```
<<<<<< HINDCAST RUN >>>>>>>                  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< FORECAST RUN >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                                              <<<<<<<< HINDCAST* >>>>>>>  < NOWCAST >   <<<<<<<<<<<<<<<<<<<<< FORECAST >>>>>>>>>>>>>>>>>
<--------- ERA5 ---------- | -------------- | ------ GFS ANALYSIS ----- | ----------- | ------ GFS FORECAST ----- | -------------------- |
<------------------------- |   6 days ago   | ------------------------- | --- NOW --- | ------------------------- | --- 10 lead days --- |
<--- MERCATOR ANALYSIS --- | -------------- | --- MERCATOR ANALYSIS --- | ----------- | --- MERCATOR FORECAST --- | -------------------- |
                             RESTART FILE       
     DAILY AVG OUTPUTS                                                             HOURLY AVG OUTPUTS
```

In simple words, the idea is to run a dynamical downscaling of ERA5 and Mercator to build the initial conditions for a forecast run. Since ERA5 has an approx. 6 days delay
this downscaling is run every day giving an unbiased restart file that is 6 days apart from real time. The forecast run consist of a small hindcast that takes account
this delay and uses the GFS analysis for filling the missing atmospheric forcing, that allows to perform a croco simulation to the current date. The second part of the forecast 
run is just a free run of croco forcing with the GFS forecast and the Mercator forecast. The idea is to postprocess the forecast run for operational uses because it will save
hourly 3D fields everyday from 6 days past to 10 days to the future. The hindcast simulation with save everyday a single daily avg output that can be used for research purposes.  