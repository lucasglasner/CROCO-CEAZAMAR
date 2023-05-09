### CROCO-CEAZAMAR Forecast/hindcast system

This repository contains the scripts and preprocessing tools for running the croco model in forecast mode. The development is part of CEAZA efforts\
to provide monitoring and forecast for the Coquimbo Region (Chile) coastal zone.

From a technical point of view, the methodology consist of a daily update of two croco simulations: (i) the croco hindcast and (ii) the croco forecast.\
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