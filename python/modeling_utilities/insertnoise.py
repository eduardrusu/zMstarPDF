# Given an image containing flux from a perfect model, as well as a value for the original (unsubtracted) sky, create an image containing pure Poisson noise, with null sky level

import numpy as np
from astropy.io import fits

model = "ilens_out_model.fits"
sky = 250000

image = fits.open(model)
data = image[0].data
data = data + sky
noise = np.random.poisson(data)
noise = noise - sky
image[0].data = noise

image.writeto(model[:-5]+"_noise.fits",clobber=True)
