### Created March 2022, I.Zelko


##########################
### Distribution functions
##########################

import numpy as np

def normalize_and_center_the_dist(m,p):
    """
    in: 
        m : mass array
        p : unnormalized distribution
    out:
        delta_m: the intervals for integration
        c_m: the centered mass points
        c_p: the centered unnormalized prob
        p_norm: uncentered normalized prob
        c_p_norm: centered normalized prob
        
    
    """
    #getting the centers of the bins
    c_m = (m[1:]+m[0:-1])/2
    ## this part could be made probably better with interpolation
    c_p = (p[1:]+p[0:-1])/2
    ## getting the intervals of integration
    delta_m = m[1:]-m[0:-1]
    ## getting the entire area
    p_total = np.sum(c_p*delta_m)
    ## normalize the distribution
    p_norm=p/p_total
    ## Normalize the centered distribution
    c_p_norm=c_p/p_total
    return delta_m, c_m, p_norm, c_p_norm
def normalize_the_dist(m,p):
    """
    in: 
        m : mass array
        p : unnormalized distribution
        p_norm: uncentered normalized prob
    """
    
    #getting the centers of the bins
    c_m = (m[1:]+m[0:-1])/2
    ## this part could be made probably better with interpolation
    c_p = (p[1:]+p[0:-1])/2
    ## getting the intervals of integration
    delta_m = m[1:]-m[0:-1]
    ## getting the entire area
    p_total = np.sum(c_p*delta_m)
    ## normalize the distribution
    p_norm=p/p_total
    return p_norm
    
def integrate_lower_95_limit(c_m,c_p_norm,delta_m):
    i=0
    integrand=0
    ## it is 0.954 because of that is how lower 2sigma comes down to a Gaussian
    while integrand<0.954:
        integrand+=c_p_norm[-i]*delta_m[-i]
        i+=1
    print("Lower 2 sigma lower limit on m is",c_m[-i], "keV")
    return c_m[-i]
def integrate_upper_95_limit(c_m,c_p_norm,delta_m):
    i=0
    integrand=0
    ## it is 0.954 because of that is how lower 2sigma comes down to a Gaussian
    while integrand<0.954:
        integrand+=c_p_norm[i]*delta_m[i]
        i+=1
    print("Upper 2 sigma upper limit on m is",c_m[i])
    return c_m[i]
def norm_and_limit(x, p):
    delta_x, c_X, p_norm, c_p_norm = normalize_and_center_the_dist(x,p )
    upper_95_limit = integrate_upper_95_limit(c_X,c_p_norm,delta_x)
    lower_95_limit =integrate_lower_95_limit(c_X,c_p_norm,delta_x)
    return lower_95_limit, upper_95_limit, p_norm

