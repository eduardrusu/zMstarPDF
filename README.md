# zMstarPDF

Calculating joint photo-z and Mstar PDFs, while coping with various systematics.

Motivated primarily by understanding massive structures along the lines of sight to time delay strong gravitational lenses so that we can attempt to make accurate time delay distance measurements, we are trying to infer photometric redshifts simultaneously (or at least self-consistently) with stellar masses, from optical and near infrared survey photometry.

Our main test data are the 5 H0LiCOW lens fields, observed in _ugriJHK_ as well as 4 _IRAC_ bands. We use a weighted counts approach to translate measured over/under-densities in these fields, with respect to a large calibration survey, into a probability distribution for the external convergence (see [_Rusu et al, 2017_] (http://adsabs.harvard.edu/abs/2016arXiv160701047R) for details). Our calibration data of choice are the CFHTLenS object catalogs, generated with `SExtractor`, with photo-zs and stellar mass esstimates from `BPZ` and `LePhare`, respectively. We have two options for inferring z and Mstar from these datasets:

* Follow the same procedure as the CFHTLenS team, so that we can simply re-use their data products. This is the approach we used in [_Rusu et al, 2017_] (http://adsabs.harvard.edu/abs/2016arXiv160701047R).
* Infer z and Mstar from the CFHTLS and H0LiCOW photometry afresh, generating MCMC samples from Pr(z,Mstar|data) directly. This is possible using the `stellarpops` code (Auger et al 2010).

This repository contains scripts and notes, as well as code used in our investigations of these options. In particular, it contains the complete code used for the analysis presented in [_Rusu et al, 2017_] (http://adsabs.harvard.edu/abs/2016arXiv160701047R). **The porting of the code is not yet complete.** A summary of our results is available on this [webpage] (http://shsuyu.github.io/H0LiCOW/site/paperIII.html).

### People

* Cristian Eduard Rusu (UC Davis)
* Chris Fassnacht (UC Davis)
* Phil Marshall (KIPAC)

### Contacts, License etc.

This is astronomy research in progress: while the contents of this repository are publically visible, they are Copyright 2015 the authors, and not available for re-use. Please cite _([Rusu et al, 2017] (http://adsabs.harvard.edu/abs/2016arXiv160701047R))_, _(Rusu et al, in preparation)_ if you need to refer to this work, and feel free to get in touch via [this repo's issues](https://github.com/eduardrusu/zMstarPDF/issues).
