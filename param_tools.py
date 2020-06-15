import numpy as np
from scipy.interpolate import interp1d, interp2d


def arc_cumulator(t, coords):
    """
    Parameters
    ----------
    t : array or None
        Parameter values associated with data coordinates. Should be sorted
        already. Pass None to use np.linspace(0, 1, n).
    coords : array
        Values of the curve at each t. 
        
    Returns
    -------
    t : array
        As above.
    cum_s : array
        Cumulative arc length on [t[0], t].
    
    Evaluates the cumulative arc length at each coordinate in sequence.
    """
    
    if np.all(t)==None:
        t = np.linspace(0, 1, coords.shape[-1])
    
    assert t.shape == coords.shape[1:], "Need same number of parameters as coordinates"
    delta_s = np.linalg.norm(np.diff(coords), axis=0)
    cum_s = np.concatenate([[0], np.cumsum(delta_s)])
    
    return t, cum_s


def r_arc_from_data(n, t, coords, interp=True, kind='linear'):
    """
    Parameters
    ----------
    n : int
        Number of points to generate.
    t : array or None
        Parameter values associated with data coordinates. Should be sorted
        already. Pass None to use np.linspace(0, 1, n).
    coords : array
        Values of the curve at each t. 
    interp : boolean
        Whether to generate random function values or not. Set to False if
        you have another way of evaluating the function.
    kind : str
        Interpolation method to be passed to scipy.interpolate.
        
    Returns
    -------
    (rand_coords : array)
        Random coordinates, if interp=True was passed.
    rand_t : array
        Parameters associated with the coordinates.
    rand_s : array
        Cumulative arc length at each coordinate.
        
    If the parameterizing function is known, use r_arc() instead, as the
    coordinates will be exactly computed from t instead of interpolating.
    """
    
    t, cum_s = arc_cumulator(t, coords)
    
    # Get random values on [0, length]
    # To test uniformity, you can use:
    # rand_s = np.linspace(0,cum_s[-1], n)
    rand_s = np.random.rand(n) * cum_s[-1]
    
    # Find corresponding t-values by interpolation
    rand_t = interp1d(cum_s, t)(rand_s)
    
    if interp==True:
        # Interpolate coordinates, e.g. if func unknown
        rand_coords = interp1d(t, coords, kind=kind)(rand_t)
        return rand_coords, rand_t, rand_s

    else:
        return rand_t, rand_s
    

def r_arc(n, func, t0, t1, precision=225):
    """
    Parameters
    ----------
    n : int
        Number of points to generate.
    func : function
        Parametric function describing the curve on which points should
        be generated.
    t0, t1 : ints or floats
        Range over which func is evaluated.
    precision : int
        Number of t-values at which func is evaluated when computing
        arc length. 
        
    Returns
    -------
    rand_coords : array
        Random coordinates.
    rand_t : array
        Parameters associated with the coordinates.
    rand_s : array
        Cumulative arc length at each coordinate.
        
    Generates random points distributed uniformly along a parametric curve.
    """
    
    t = np.linspace(t0, t1, precision)
    coords = func(t)
    
    rand_t, rand_s = r_arc_from_data(n, t, coords, interp=False)
    rand_coords = func(rand_t)
    
    return rand_coords, rand_t, rand_s


def arc_length(func, t0, t1, precision=225):
    """
    Parameters
    ----------
    func : function
        Parametric function describing the curve on which points should
        be generated.
    t0, t1: ints or floats
        Range over which func is evaluated.
    precision : int
        Number of t-values at which func is evaluated.
        
    Returns
    -------
    length : float
        Estimate of surface area.
    
    Convenience function to evaluate total arc length over given range.
    """
    
    t = np.linspace(t0, t1, precision)
    coords = func(t)
    
    length = arc_cumulator(t, coords)[1][-1]
    
    return length


def surface_cumulator(t, u, coords):
    """
    Parameters
    ----------
    t : array or None
        Parameter values associated with data coordinates. Should be sorted
        already. Pass None to use np.linspace(0, 1, n).
    u : array or None
        Parameter values associated with data coordinates. Should be sorted
        already. Pass None to use np.linspace(0, 1, n).
    coords : array
        Values of the curve at each t, u pair. 
        
    Returns
    -------
    t, u : arrays
        As above.
    cum_S_t : array
        Cumulative surface area on [t[0], t], all u.
    cum_S_u : array
        Cumulative surface area on all t, [u[0], u].
    
    Evaluates the cumulative surface area at each coordinate.
    """
    
    if np.all(t)==None:
        t, _ = np.meshgrid(np.linspace(0, 1, coords.shape[-2]), np.linspace(0, 1, coords.shape[-1]))
    if np.all(u)==None:
        _, u = np.meshgrid(np.linspace(0, 1, coords.shape[-2]), np.linspace(0, 1, coords.shape[-1]))
    
    assert t.shape == u.shape == coords.shape[1:], "Need same number of parameters as coordinates"
    delta_t_temp = np.diff(coords, axis=2)
    delta_u_temp = np.diff(coords, axis=1)

    # Pad with zeros so that small rand_S can still be interpd
    delta_t = np.zeros(coords.shape)
    delta_u = np.zeros(coords.shape)

    delta_t[:coords.shape[0], :coords.shape[1], 1:coords.shape[2]] = delta_t_temp
    delta_u[:coords.shape[0], 1:coords.shape[1], :coords.shape[2]] = delta_u_temp

    # Area of each parallelogram
    delta_S = np.linalg.norm(np.cross(delta_t, delta_u, 0, 0), axis=2)

    cum_S_t = np.cumsum(delta_S.sum(axis=0))
    cum_S_u = np.cumsum(delta_S.sum(axis=1))
    
    return t, u, cum_S_t, cum_S_u


def r_surface_from_data(n, t, u, coords, interp=True, kind='linear'):   
    """
    Parameters
    ----------
    n : int
        Number of points to generate.
    t : array or None
        Parameter values associated with data coordinates. Should be sorted
        already. Pass None to use np.linspace(0, 1, n).
    u : array or None
        Parameter values associated with data coordinates. Should be sorted
        already. Pass None to use np.linspace(0, 1, n).
    coords : array
        Values of the curve at each t, u pair. 
    interp : boolean
        Whether to generate random function values or not. Set to false if you
        have another way of evaluating the function.
    kind : str
        Interpolation method to be passed to scipy.interpolate.
        
    Returns
    -------
    (rand_coords : array)
        Random coordinates, if interp=True was passed.
    rand_t : array
        t-values associated with the coordinates.
    rand_u : array
        u_values associated with the coordinates.
    rand_S_t : array
        Cumulative area at each t-value over full range of u.
    rand_S_u : array
        Cumulative area at each u-value over full range of t.
        
    If the parameterizing function is known, use r_surface() instead, as the
    coordinates will be exactly computed from t instead of interpolating.
    """
    
    t, u, cum_S_t, cum_S_u = surface_cumulator(t, u, coords)

    # Random values
    rand_S_t = np.random.rand(n) * cum_S_t[-1]
    rand_S_u = np.random.rand(n) * cum_S_u[-1]

    # Find corresponding t-values by interpolation
    rand_t = interp1d(cum_S_t, t[0,:])(rand_S_t)
    rand_u = interp1d(cum_S_u, u[:,0])(rand_S_u)

    if interp==True:
        # Interpolate coordinates, e.g. if func unknown
        
        rand_coords = np.empty([coords.shape[0], n])
        
        # One axis at a time, or else scipy throws dim mismatch
        for i in range(coords.shape[0]):
            f = interp2d(t, u, coords[i], kind=kind)
            
            # One point at time, or else scipy does a meshgrid
            for j in range(n):
                rand_coords[i,j] = f(rand_t[j], rand_u[j])
                
        return rand_coords, rand_t, rand_u, rand_S_t, rand_S_u

    else:
        return rand_t, rand_u, rand_S_t, rand_S_u


def r_surface(n, func, t0, t1, u0, u1, t_precision=25, u_precision=25):
    """
    Parameters
    ----------
    n : int
        Number of points to generate.
    func : function
        Parametric function describing the curve on which points should
        be generated.
    t0, t1, u0, u1 : ints or floats
        Range over which func is evaluated.
    t_precision, u_precision : ints
        Number of t-values at which func is evaluated when computing
        surface area. 
       
    Returns
    -------
    rand_coords : array
        Random coordinates.
    rand_t : array
        t-values associated with the coordinates.
    rand_u : array
        u_values associated with the coordinates.
    rand_S_t : array
        Cumulative area at each t-value over full range of u.
    rand_S_u : array
        Cumulative area at each u-value over full range of t.
        
    Generates random points distributed uniformly over a parametric surface.
    """
    
    t, u = np.meshgrid(np.linspace(t0, t1, t_precision), np.linspace(u0, u1, u_precision))
    coords = func(t, u)
    
    rand_t, rand_u, rand_S_t, rand_S_u = r_surface_from_data(n, t, u, coords, interp=False)
    rand_coords = func(rand_t, rand_u)
    
    return rand_coords, rand_t, rand_u, rand_S_t, rand_S_u


def surface_area(func, t0, t1, u0, u1, t_precision=25, u_precision=25):
    """
    Parameters
    ----------
    func : function
        Parametric function describing the curve on which points should
        be generated.
    t0, t1, u0, u1 : ints or floats
        Range over which func is evaluated.
    t_precision, u_precision : ints
        Number of t- and u-values at which func is evaluated.
        
    Returns
    -------
    area : float
        Estimate of surface area.
    
    Convenience function to evaluate total surface area over given range.    
    """
    
    t, u = np.meshgrid(np.linspace(t0, t1, t_precision), np.linspace(u0, u1, u_precision))
    coords = func(t, u)
    
    area = surface_cumulator(t, u, coords)[3][-1]
    
    return area