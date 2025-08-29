using SatelliteToolbox
SpaceIndices.init()

tstart = (2016, 6, 14, 00, 0, 0)
tend =  (2016, 6, 15, 00, 0, 0)
interval = 10
dr = DateTime(tstart...):Minute(interval):DateTime(tend...)
jdate = date_to_jd.(dr)
heights = [80e3, 90e3, 100e3, 150e3]
lon = 19.23*pi/180
lat = 19.23*pi/180

atm = [AtmosphericModels.nrlmsise00(jd,                  # Julian Day
                                    h,                   # Altitude [m]
                                    la,         # Latitude [rad]
                                    lo         # Longitude [rad]
                                    ) for jd in jdate, h in heights, lo in lon, la in lat]

using Plots
n2dens = [a.N2_number_density for a in atm[:, 2]]
plot(diff(n2dens)./n2dens[1:end-1])

function getmsis_h(time, heights, long, lat)
    jdate = date_to_jd.(time)
    atm = [AtmosphericModels.nrlmsise00(jdate,                  # Julian Day
                                    heigths,                   # Altitude [m]
                                    69.58*pi/180,         # Latitude [rad]
                                    19.23*pi/180,         # Longitude [rad]
                                    ) for h in heights]

end




