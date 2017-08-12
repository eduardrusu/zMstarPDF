import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as so

def find_confidence_interval(x, pdf, confidence_level):
    return pdf[pdf > x].sum() - confidence_level

def density_contour(xdata, ydata, nbins_x, nbins_y, ax=None, **contour_kwargs):
    """ Create a density contour plot.

    Parameters
    ----------
    xdata : numpy.ndarray
    ydata : numpy.ndarray
    nbins_x : int
        Number of bins along x dimension
    nbins_y : int
        Number of bins along y dimension
    ax : matplotlib.Axes (optional)
        If supplied, plot the contour to this axis. Otherwise, open a new figure
    contour_kwargs : dict
        kwargs to be passed to pyplot.contour()
    """

    H, xedges, yedges = np.histogram2d(xdata, ydata, bins=(nbins_x,nbins_y), normed=True)
    x_bin_sizes = (xedges[1:] - xedges[:-1]).reshape((1,nbins_x))
    y_bin_sizes = (yedges[1:] - yedges[:-1]).reshape((nbins_y,1))

    pdf = (H*(x_bin_sizes*y_bin_sizes))

    one_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.68))
    two_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.95))
    three_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.99))
    levels = [one_sigma, two_sigma, three_sigma]

    X, Y = 0.5*(xedges[1:]+xedges[:-1]), 0.5*(yedges[1:]+yedges[:-1])
    Z = pdf.T

    if ax == None:
        #plt.rcParams['contour.negative_linestyle'] = 'solid'
        contour = plt.contour(X, Y, Z, levels=levels, origin="lower", colors=('red','magenta','blue'), linewidths = 0.7)
    else:
        #plt.rcParams['contour.negative_linestyle'] = 'solid'
        contour = ax.contour(X, Y, Z, levels=levels, origin="lower", colors=('red','magenta','blue'), linewidths = 0.7)

    return contour

def density_contourdash(xdata, ydata, nbins_x, nbins_y, ax=None, **contour_kwargs):
    """ Create a density contour plot.

    Parameters
    ----------
    xdata : numpy.ndarray
    ydata : numpy.ndarray
    nbins_x : int
        Number of bins along x dimension
    nbins_y : int
        Number of bins along y dimension
    ax : matplotlib.Axes (optional)
        If supplied, plot the contour to this axis. Otherwise, open a new figure
    contour_kwargs : dict
        kwargs to be passed to pyplot.contour()
    """

    H, xedges, yedges = np.histogram2d(xdata, ydata, bins=(nbins_x,nbins_y), normed=True)
    x_bin_sizes = (xedges[1:] - xedges[:-1]).reshape((1,nbins_x))
    y_bin_sizes = (yedges[1:] - yedges[:-1]).reshape((nbins_y,1))

    pdf = (H*(x_bin_sizes*y_bin_sizes))

    one_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.68))
    two_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.95))
    three_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.99))
    levels1 = [one_sigma]
    levels2 = [two_sigma]
    levels3 = [three_sigma]

    X, Y = 0.5*(xedges[1:]+xedges[:-1]), 0.5*(yedges[1:]+yedges[:-1])
    Z = pdf.T

    if ax == None:
        #plt.rcParams['contour.negative_linestyle'] = 'solid'
        contour1 = plt.contour(X, Y, Z, levels=levels1, origin="lower", colors=('black'), linewidths = 0.7, linestyles = 'solid')
    else:
        #plt.rcParams['contour.negative_linestyle'] = 'solid'
        contour1 = ax.contour(X, Y, Z, levels=levels1, origin="lower", colors=('black'), linewidths = 0.7, linestyles = 'solid')

    if ax == None:
        #plt.rcParams['contour.negative_linestyle'] = 'solid'
        contour2 = plt.contour(X, Y, Z, levels=levels2, origin="lower", colors=('black'), linewidths = 0.7, linestyles = 'dashed')
    else:
        #plt.rcParams['contour.negative_linestyle'] = 'solid'
        contour2 = ax.contour(X, Y, Z, levels=levels2, origin="lower", colors=('black'), linewidths = 0.7, linestyles = 'dashed')

    if ax == None:
        #plt.rcParams['contour.negative_linestyle'] = 'solid'
        contour3 = plt.contour(X, Y, Z, levels=levels3, origin="lower", colors=('black'), linewidths = 0.7, linestyles = 'dotted')
    else:
        #plt.rcParams['contour.negative_linestyle'] = 'solid'
        contour3 = ax.contour(X, Y, Z, levels=levels3, origin="lower", colors=('black'), linewidths = 0.7, linestyles = 'dotted')

    return contour1,contour2,contour3

def density_contourdashmagenta(xdata, ydata, nbins_x, nbins_y, ax=None, **contour_kwargs):
    """ Create a density contour plot.
        
        Parameters
        ----------
        xdata : numpy.ndarray
        ydata : numpy.ndarray
        nbins_x : int
        Number of bins along x dimension
        nbins_y : int
        Number of bins along y dimension
        ax : matplotlib.Axes (optional)
        If supplied, plot the contour to this axis. Otherwise, open a new figure
        contour_kwargs : dict
        kwargs to be passed to pyplot.contour()
        """
    
    H, xedges, yedges = np.histogram2d(xdata, ydata, bins=(nbins_x,nbins_y), normed=True)
    x_bin_sizes = (xedges[1:] - xedges[:-1]).reshape((1,nbins_x))
    y_bin_sizes = (yedges[1:] - yedges[:-1]).reshape((nbins_y,1))
    
    pdf = (H*(x_bin_sizes*y_bin_sizes))
    
    one_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.68))
    two_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.95))
    three_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.99))
    levels1 = [one_sigma]
    levels2 = [two_sigma]
    levels3 = [three_sigma]
    
    X, Y = 0.5*(xedges[1:]+xedges[:-1]), 0.5*(yedges[1:]+yedges[:-1])
    Z = pdf.T
    
    if ax == None:
        #plt.rcParams['contour.negative_linestyle'] = 'solid'
        contour1 = plt.contour(X, Y, Z, levels=levels1, origin="lower", colors=('magenta'), linewidths = 0.7, linestyles = 'solid')
    else:
        #plt.rcParams['contour.negative_linestyle'] = 'solid'
        contour1 = ax.contour(X, Y, Z, levels=levels1, origin="lower", colors=('magenta'), linewidths = 0.7, linestyles = 'solid')
    
    if ax == None:
        #plt.rcParams['contour.negative_linestyle'] = 'solid'
        contour2 = plt.contour(X, Y, Z, levels=levels2, origin="lower", colors=('magenta'), linewidths = 0.7, linestyles = 'dashed')
    else:
        #plt.rcParams['contour.negative_linestyle'] = 'solid'
        contour2 = ax.contour(X, Y, Z, levels=levels2, origin="lower", colors=('magenta'), linewidths = 0.7, linestyles = 'dashed')
    
    if ax == None:
        #plt.rcParams['contour.negative_linestyle'] = 'solid'
        contour3 = plt.contour(X, Y, Z, levels=levels3, origin="lower", colors=('magenta'), linewidths = 0.7, linestyles = 'dotted')
    else:
        #plt.rcParams['contour.negative_linestyle'] = 'solid'
        contour3 = ax.contour(X, Y, Z, levels=levels3, origin="lower", colors=('magenta'), linewidths = 0.7, linestyles = 'dotted')
    
    return contour1,contour2,contour3

def density_contourdashblue(xdata, ydata, nbins_x, nbins_y, ax=None, **contour_kwargs):
    """ Create a density contour plot.
        
        Parameters
        ----------
        xdata : numpy.ndarray
        ydata : numpy.ndarray
        nbins_x : int
        Number of bins along x dimension
        nbins_y : int
        Number of bins along y dimension
        ax : matplotlib.Axes (optional)
        If supplied, plot the contour to this axis. Otherwise, open a new figure
        contour_kwargs : dict
        kwargs to be passed to pyplot.contour()
        """
    
    H, xedges, yedges = np.histogram2d(xdata, ydata, bins=(nbins_x,nbins_y), normed=True)
    x_bin_sizes = (xedges[1:] - xedges[:-1]).reshape((1,nbins_x))
    y_bin_sizes = (yedges[1:] - yedges[:-1]).reshape((nbins_y,1))
    
    pdf = (H*(x_bin_sizes*y_bin_sizes))
    
    one_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.68))
    two_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.95))
    three_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.99))
    levels1 = [one_sigma]
    levels2 = [two_sigma]
    levels3 = [three_sigma]
    
    X, Y = 0.5*(xedges[1:]+xedges[:-1]), 0.5*(yedges[1:]+yedges[:-1])
    Z = pdf.T
    
    if ax == None:
        #plt.rcParams['contour.negative_linestyle'] = 'solid'
        contour1 = plt.contour(X, Y, Z, levels=levels1, origin="lower", colors=('blue'), linewidths = 0.7, linestyles = 'solid')
    else:
        #plt.rcParams['contour.negative_linestyle'] = 'solid'
        contour1 = ax.contour(X, Y, Z, levels=levels1, origin="lower", colors=('blue'), linewidths = 0.7, linestyles = 'solid')
    
    if ax == None:
        #plt.rcParams['contour.negative_linestyle'] = 'solid'
        contour2 = plt.contour(X, Y, Z, levels=levels2, origin="lower", colors=('blue'), linewidths = 0.7, linestyles = 'dashed')
    else:
        #plt.rcParams['contour.negative_linestyle'] = 'solid'
        contour2 = ax.contour(X, Y, Z, levels=levels2, origin="lower", colors=('blue'), linewidths = 0.7, linestyles = 'dashed')
    
    if ax == None:
        #plt.rcParams['contour.negative_linestyle'] = 'solid'
        contour3 = plt.contour(X, Y, Z, levels=levels3, origin="lower", colors=('blue'), linewidths = 0.7, linestyles = 'dotted')
    else:
        #plt.rcParams['contour.negative_linestyle'] = 'solid'
        contour3 = ax.contour(X, Y, Z, levels=levels3, origin="lower", colors=('blue'), linewidths = 0.7, linestyles = 'dotted')
    
    return contour1,contour2,contour3

def density_contourfill(xdata, ydata, nbins_x, nbins_y, ax=None, **contour_kwargs):
    """ Create a density contour plot.

    Parameters
    ----------
    xdata : numpy.ndarray
    ydata : numpy.ndarray
    nbins_x : int
        Number of bins along x dimension
    nbins_y : int
        Number of bins along y dimension
    ax : matplotlib.Axes (optional)
        If supplied, plot the contour to this axis. Otherwise, open a new figure
    contour_kwargs : dict
        kwargs to be passed to pyplot.contour()
    """

    H, xedges, yedges = np.histogram2d(xdata, ydata, bins=(nbins_x,nbins_y), normed=True)
    x_bin_sizes = (xedges[1:] - xedges[:-1]).reshape((1,nbins_x))
    y_bin_sizes = (yedges[1:] - yedges[:-1]).reshape((nbins_y,1))

    pdf = (H*(x_bin_sizes*y_bin_sizes))

    one_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.68))
    two_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.95))
    three_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.99))
    levels1 = [one_sigma,0]
    levels2 = [two_sigma,one_sigma]
    levels3 = [three_sigma,two_sigma]

    X, Y = 0.5*(xedges[1:]+xedges[:-1]), 0.5*(yedges[1:]+yedges[:-1])
    Z = pdf.T

    if ax == None:
        #plt.rcParams['contour.negative_linestyle'] = 'solid'
        contour1 = plt.contourf(X, Y, Z, levels=levels1, origin="lower", colors=('black'))
    else:
        #plt.rcParams['contour.negative_linestyle'] = 'solid'
        contour1 = ax.contourf(X, Y, Z, levels=levels1, origin="lower", colors=('black'))

    if ax == None:
        #plt.rcParams['contour.negative_linestyle'] = 'solid'
        contour2 = plt.contourf(X, Y, Z, levels=levels2, origin="lower", colors=('gray'))
    else:
        #plt.rcParams['contour.negative_linestyle'] = 'solid'
        contour2 = ax.contourf(X, Y, Z, levels=levels2, origin="lower", colors=('gray'))

    if ax == None:
        #plt.rcParams['contour.negative_linestyle'] = 'solid'
        contour3 = plt.contourf(X, Y, Z, levels=levels3, alpha=0.5, origin="lower", colors=('gray'))
    else:
        #plt.rcParams['contour.negative_linestyle'] = 'solid'
        contour3 = ax.contourf(X, Y, Z, levels=levels3, alpha=0.5, origin="lower", colors=('gray'))

    return contour1,contour2,contour3



#def test_density_contour():
#    norm = np.random.normal(10., 15., size=(12540035, 2))
#    density_contour(norm[:,0], norm[:,1], 100, 100)
#    plt.show()

#test_density_contour()