import numpy as np

float_type = 'f8'

# # Get number of vertical levels and size from .ini file
# with open('moser180.ini') as f:
#     for line in f:
#         if (line.split('=')[0]=='ktot'):
#             kmax = int(line.split('=')[1])
#         if (line.split('=')[0]=='zsize'):

kmax=128
zsize=61.72

# define the variables
z = np.zeros(kmax+1)
u = np.zeros(kmax+1)
s = np.zeros(kmax+1)

# create non-equidistant grid
alpha = 0.85
for k in range(0,kmax+1):
    #eta  = -1. + 2.*((k+1)-0.5) / kmax
    eta  = -1. + 2.*((k+0)-0.0) / kmax
    z[k] = zsize / (2.*alpha) * np.tanh(eta*0.5*(np.log(1.+alpha) - np.log(1.-alpha))) + 0.5*zsize
    s[k] = z[k]

    print(z[k])

# # create initial parabolic shape
# dpdxls = -1.5e-6
# visc   =  1.0e-5
# for k in range(kmax):
#     u[k] = 1./(2.*visc)*dpdxls*(z[k]**2. - zsize*z[k])

# nc_file = nc.Dataset("moser180_input.nc", mode="w", datamodel="NETCDF4", clobber=False)

# nc_file.createDimension("z", kmax)
# nc_z  = nc_file.createVariable("z" , float_type, ("z"))

# nc_group_init = nc_file.createGroup("init");
# nc_u = nc_group_init.createVariable("u", float_type, ("z"))
# nc_s = nc_group_init.createVariable("s", float_type, ("z"))

# nc_z[:] = z[:]
# nc_u[:] = u[:]
# nc_s[:] = s[:]

# nc_file.close()

