import h5py
from astropy import wcs
import numpy as np
import math
import copy


# generate input data and input coords
def setup_data(mapcenter, mapsize, beamsize_fwhm, num_samples, num_sources, channels):
    lon_scale = np.cos(np.radians(mapcenter[1]))
    scale_size = 1.0
    map_l, map_r = (
        mapcenter[0] - scale_size * mapsize[0] / 2. / lon_scale,
        mapcenter[0] + scale_size * mapsize[0] / 2. / lon_scale
    )
    map_b, map_t = mapcenter[1] - scale_size * mapsize[1] / 2., mapcenter[1] + scale_size * mapsize[1] / 2.

    # coordinates are drawn from a uniform distribution
    xcoords = np.random.uniform(map_l, map_r, num_samples).astype(np.float32)
    ycoords = np.random.uniform(map_b, map_t, num_samples).astype(np.float32)

    # Create a dataset under group coords.
    c1 = g1.create_dataset("xcoords", data=xcoords)
    c2 = g1.create_dataset("ycoords", data=ycoords)
    #==========================Generate signal==============================
    # add Gaussian noise
    mock_data = np.random.normal(0., 1., len(xcoords)).astype(np.float32)

    beamsize_sigma = beamsize_fwhm / np.sqrt(8 * np.log(2))

    # put in artifical point source, with random amplitudes
    # we'll assume a Gaussian-shaped PSF

    def gauss2d(x, y, x0, y0, A, s):
        return A * np.exp(-((x - x0) ** 2 + (y - y0) ** 2) / 2. / s ** 2)

    
    for i in range(channels):
        data = copy.deepcopy(mock_data)
        for _ in range(num_sources):
            sou_x = np.random.uniform(map_l, map_r, 1).astype(np.float32)
            sou_y = np.random.uniform(map_b, map_t, 1).astype(np.float32)
            A = np.random.uniform(0, 100, 1).astype(np.float32)
            point_source = gauss2d(xcoords, ycoords, sou_x, sou_y, A, beamsize_sigma)
            data += point_source
        signal = data[:, np.newaxis]  # need dummy spectral axis
        # Create a dataset under group signal.
        dset = g2.create_dataset('signal_%d' %i, data=signal)
        # Create a dataset under group signal.
        n1 = dset.attrs["num_samples"] = num_samples

    return xcoords, ycoords


#==========================Creat HDF5 File==============================
file_group = h5py.File('input2.hdf5', 'w')

# Create a group and a dataset under group "coords".
g1 = file_group.create_group("coords")
g2 = file_group.create_group("signal")

#==========================Set Related Parameter=========================
# set input data and input coords parameters
mapcenter = 60., 30.  # all in degrees
map_size = 5.
beam_size = 300
# pixelsize = 0.01
mapsize = map_size, map_size
beamsize_fwhm = 2 * beam_size / 3600.
# a good pixel size is a third of the FWHM of the PSF (avoids aliasing)
# pixsize = beamsize_fwhm / 3.
num_samples = 10000000
num_sources = 10
channels = 4

#==========================Create Channels files in HDF5 File==============================
# get data and coords
xcoords, ycoords = setup_data(
    mapcenter, mapsize, beamsize_fwhm, num_samples, num_sources, channels
)

# Save and exit the file.
file_group.close()
