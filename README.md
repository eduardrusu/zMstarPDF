# zMstarPDF

Calculating joint photo-z and Mstar PDFs, while coping with various systematics.

Motivated primarily by understanding massive structures along the lines of sight to time delay strong gravitational lenses so that we can attempt to make accurate time delay distance measurements, we are trying to infer photometric redshifts simultaneously (or at least self-consistently) with stellar masses, from optical and near infra-red survey photometry.

Our main test data are the 5 H0LiCOW lens fields, observed in _ugriJK_. Our calibration data of choice are the CFHTLS object catalogs, generated with `SExtractor` and `LePhare` (for the photo-zs and stellar masses). We have two options for inferring z and Mstar from these datasets:

* Follow the same `LePhare` procedure as the CFHTLenS team, so that we can simply re-use their data products.
* Infer z and Mstar from the CFHTLS and H0LiCOW photometry afresh, generating MCMC samples from Pr(z,Mstar|data) directly. This is possible using the `stellarpops` code (Auger et al 2010).

This repository contains scripts and notes from our investigations of these options.

### People

* Edi Rusu (UC Davis)
* Chris Fassnacht (UC Davis)
* Phil Marshall (KIPAC)

### Contacts, License etc.

This is astronomy research in progress: while the contents of this repository are publically visible, they are Copyright 2015 the authors, and not available for re-use. Please cite _(Rusu et al, in preparation)_ if you need to refer to this work, and feel free to get in touch via [this repo's issues](https://github.com/eduardrusu/zMstarPDF/issues).
