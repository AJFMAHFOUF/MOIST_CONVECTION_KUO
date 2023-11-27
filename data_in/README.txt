Documentation of the file: profiles.dat
---------------------------------------

Jean-François Mahfouf (23/11/2023)


This file contains 1535 vertical profiles from the ECMWF model in a configuration with
60 vertical levels and a TL511 truncation (probably CY23R4) that I must have created
in 2001 just before leaving for Canada (see Mahfouf et al. (2005) at QJRMS). 
As it's been more than 20 years, various uncertainties remain on the exact content of this file.

These are probably short-range forecasts (6 h). The profiles were sampled
on 8 latitude bands around the equator between 1.23°N and -1.23°S (every 40 km 
corresponding to the Gaussian grid spacing). On each latitude circle, approximately
200 longitudes (out of the 1023 in the Gaussian grid). This sampling is not regular and
must therefore correspond to points where the deep convection pattern has been activated. 

The time step of this model was 900 s (to estimate evolved variables from trends).

Latitude 1: number of points = 207
Latitude 2: number of points = 214
Latitude 3: number of points = 195
Latitude 4: number of points = 188
Latitude 5: number of points = 191
Latitude 6: number of points = 186
Latitude 7: number of points = 187
Latitude 8: number of points = 167

The file is structured as follows for each profile:

1) Header: 20-character string indicating the profile number (1 to 1535)

2) 5 real numbers: Latitude (degrees), Longitude (degrees), Land/sea index (1 for land and 0 for sea),
Surface temperature (K), Surface pressure (Pa)

3) 4 real numbers whose meaning remains unknow to me, all profiles have the same values:
(2.71) - (-39.72) - (-0.039) - (0.007)

4) 7 tables of 60 values (vertical profiles) ordered from the top of the model to the surface
- Column 1: Temperature (K)
- Column 2: Specific humidity (kg/kg)
- Column 3: Pressure at full levels (Pa)
- Column 4: Pressure at half levels (Pa) - surface pressure is at half level values
- Column 5: Vertical wind speed (m/s)
- Column 6: Horizontal wind speed - zonal component (m/s)
- Column 7: Horizontal wind speed - meridian component (m/s)


Documentation of the file : input_model_mean_values.dat
-------------------------------------------------------

This file contains mean convective profiles on 60 levels of:
- Temperature (K) 
- Specific humidity (kg/kg)
- Temperature tendency (K/s)
- Specific humidity tendency ((kg/kg)/s) 

They have been obtained from an average of 53 model profiles
provided by P. Lopez (ECMWF) in the 1D-Var system for the assimilation
of surface rainrates (devised in 2006 for HIRLAM)


