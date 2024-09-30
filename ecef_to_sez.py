# ecef_to_sez.py
#
# Usage: python3 script_name.py arg1 arg2 ...
#  Text explaining script usage
# Parameters:
# 
#  ...
# Output:
#  A description of the script output
#
# Written by Brad Denby
# Other contributors: Josh Smith
#
# Optional license statement, e.g., See the LICENSE file for the license.

# import Python modules
# e.g., import math # math module
import sys # argv
import math

# "constants"
R_E_KM = 6378.1368
E_E=0.081819221456

# helper functions

## function description
def calc_denom(ecc, lat_rad):
  return math.sqrt(1.0-(ecc**2)*(math.sin(lat_rad)**2))
#   pass

# initialize script arguments
o_x_km=float('nan')
o_y_km=float('nan')
o_z_km=float('nan')
x_km=float('nan')
y_km=float('nan')
z_km=float('nan')

# parse script arguments
if len(sys.argv)==7:
    o_x_km = float(sys.argv[1])
    o_y_km = float(sys.argv[2])
    o_z_km = float(sys.argv[3])
    x_km = float(sys.argv[4])
    y_km = float(sys.argv[5]) 
    zz_km = float(sys.argv[6])
    ...
else:
    print(\
        'Usage: '\
            'python3 o_x_km o_y_km o_z_km x_km y_km zz_km ...'\
                )
    exit()

# write script below this line

# First step, determine ECEF vector from the station to the object
ecef_x_km=x_km-o_x_km
ecef_y_km=y_km-o_y_km
ecef_z_km=zz_km-o_z_km
# plug origin values into ecef_to_llh.py to get lat lon and hae
lon_rad = math.atan2(o_y_km,o_x_km)
lon_deg = lon_rad*180.0/math.pi

lat_rad = math.asin(o_z_km/math.sqrt(o_x_km**2+o_y_km**2+o_z_km**2))
r_lon_km = math.sqrt(o_x_km**2+o_y_km**2)
prev_lat_rad = float('nan')

# iteratively find latitude
c_E = float('nan')
count = 0
while (math.isnan(prev_lat_rad) or abs(lat_rad-prev_lat_rad)>10e-7) and count<5:
  denom = calc_denom(E_E,lat_rad)
  c_E = R_E_KM/denom
  prev_lat_rad = lat_rad
  lat_rad = math.atan((o_z_km+c_E*(E_E**2)*math.sin(lat_rad))/r_lon_km)
  count = count+1
  
# calculate hae
hae_km = r_lon_km/math.cos(lat_rad)-c_E
lat_deg=lat_rad*(180/math.pi)

#initialize cos and sin for lat and lon
c_phi=math.cos(lat_rad)
s_phi=math.sin(lat_rad)
c_th=math.cos(lon_rad)
s_th=math.sin(lon_rad)
# first, calculate the first matrix multiplication with Rz(theta) and Ry(90-phi)
R=[[s_phi*c_th+0+0,s_phi*s_th+0+0,0+0-c_phi],
   [0-s_th+0,0+c_th+0,0+0+0],
   [c_phi*c_th+0+0,c_phi*s_th+0+0,0+0+s_phi]]
#initialize dummy variables for the matrix entries
[[x1,x2,x3],[y1,y2,y3],[z1,z2,z3]]=R
r_sez=[ecef_x_km*x1+ecef_y_km*x2+ecef_z_km*x3,
       ecef_x_km*y1+ecef_y_km*y2+ecef_z_km*y3,
       ecef_x_km*z1+ecef_y_km*z2+ecef_z_km*z3]

s_km,e_km,zz_km=r_sez

print(s_km)
print(e_km)
print(z_km)
