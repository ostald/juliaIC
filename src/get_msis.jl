using SatelliteToolbox
SpaceIndices.init()

function msis(times, heights, lats, longs)
    jdate = date_to_jd.(times)
    atm = [AtmosphericModels.nrlmsise00(
            jd,         # Julian Day
            h,          # Altitude [m]
            lat,         # Latitude [rad]
            long          # Longitude [rad]
            ) for jd in jdate, h in heights, lat in lats, long in longs];
    return atm
end





