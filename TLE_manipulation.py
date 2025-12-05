"""
Process TLE lines to get the initial state vector 
to propagate whithin URANIO. 
"""

import sgp4
from sgp4.api import WGS72, Satrec
from sgp4.api import jday
from sgp4.api import days2mdhms
from datetime import datetime


# Example: NewSat - Satellogic 04/12/2025
tle_line1 = '1 59122U 24043AA  25338.43666971  .00033378  00000-0  11161-2 0  9991'
tle_line2 = '2 59122  97.5260 111.0397 0008144 190.8999 169.2066 15.30890055 97499'

newsat = Satrec.twoline2rv(tle_line1, tle_line2, WGS72)
yr = newsat.epochyr
days = newsat.epochdays
month, day, hour, minute, second = days2mdhms(yr, days)

jd, fr = jday(2000+int(yr), month, day, hour, minute, second)
e, r, v = newsat.sgp4(jd, fr)

print(r)
print(v)