This is intended to be a R packages to be published along the manuscript. It will be a all-in-one package for GxE. 
All the necessary subfunctions are already tested in-house. Wrapping with necessry R packages will be the next step.
The main functions including:
a). Automatically download weather data , i.e. temperature, from NOAA website https://www.ncdc.noaa.gov/cdo-web/datatools/findstation,
    and daylength from http://aa.usno.navy.mil/data/docs/Dur_OneYear.php based on the locations of tested fields (latitude and longitude)
b). Combine the downloaded data (temperature and daylength) into a master environment parameters according to planting date.
c). Trait coorelations among tested fields.
d). Find the highest coorelation bewteen trait and environment environment parameter.
e). Cross validation for using environment parameter to predict trais in untested environment.
