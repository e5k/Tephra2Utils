# TephraFits

This is a collection of Matlab functions that were developed around [Tephra2](https://github.com/ljc-geo/tephra2). 

## Content

### plotT2.m 
Plots the output of Tephra2 on a map. It can be used with the functions *plot_google_map* and *utm2ll* (both available on the Matlab File Exchange) to plot the output on a Google Map background. This is described [here](https://e5k.github.io/codes/utilities/2017/12/10/plot_tephra2/).

### runT2.m 
Runs Tephra2 and allows sensitivity analyses of parameters. *This is still experimental*.

### plotBetaPlume.m 
Plots the mass distribution in the column using a Beta PDF, which is the model used in Tephra2.

### plotTGSD.m 
Plots a gaussian TGSD.

### processT2Inversion.m
Process the output of inversion runs described [here](https://e5k.github.io/codes/utilities/2018/06/06/inversion/).

### plotT2Inversion.m
A script to pliot the results of inversion runs.

### plotWind.m
Plots 3-columns, tab-delimited ascii wind profiles.

## Versions

#### v1.0 (pre-Nov 18)
In Novembre 2018, the core of the inversion method slightly changed. To retrieve the functions required to post-process runs earlier than than, download [this version](https://github.com/e5k/Tephra2Utils/archive/v1.0.zip) instead. 