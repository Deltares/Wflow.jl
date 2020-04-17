
function lattometres(lat::Float64)
    m1 = 111132.92     #latitude calculation term 1
    m2 = -559.82       #latitude calculation term 2
    m3 = 1.175         #latitude calculation term 3
    m4 = -0.0023       #latitude calculation term 4
    p1 = 111412.84     #longitude calculation term 1
    p2 = -93.5         #longitude calculation term 2
    p3 = 0.118         #longitude calculation term 3

    #Calculate the length of a degree of latitude and longitude in meters
    latlen =
        m1 +
        (m2 * cosd(2.0 * lat)) +
        (m3 * cosd(4.0 * lat)) +
        (m4 * cosd(6.0 * lat))
    longlen = (p1 * cosd(lat)) + (p2 * cosd(3.0 * lat)) + (p3 * cosd(5.0 * lat))

    return longlen, latlen
end
