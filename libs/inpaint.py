"""Algorithms for filling gaps in ND arrays"""
import numpy as np
from scipy import ndimage as nd

def fill_ndimage(data,invalid=None):
    """Replace the value of invalid 'data' cells (indicated by 'invalid')
    by the value of the nearest valid data cell
    Parameters
    ----------
    data: numpy array of any dimension
    invalid: a binary array of same shape as 'data'. True cells set where data
    value should be replaced.
    If None (default), use: invalid = np.isnan(data)
    Returns
    -------
    Return a filled array.
    Credits
    -------
    http://stackoverflow.com/a/9262129
    """
    if invalid is None: invalid = np.isnan(data)
    ind = nd.distance_transform_edt(invalid, return_distances=False, return_indices=True)
    return data[tuple(ind)]

def fill_crystalgrowth(data,invalid):
    """Crystal-growth filler
    Credits
    -------
    http://stackoverflow.com/a/5556426
    """
    # invert invalid so that masked values are True
    invalid = ~invalid
    shape = data.shape
    dim = len(shape)
    slcs = [slice(None)]*dim
    while np.any(~invalid): # as long as there are any False's in invalid
       for i in xrange(dim): # do each axis
          # make slices to shift view one element along the axis
          slcs1 = slcs[:]
          slcs2 = slcs[:]
          slcs1[i] = slice(0, -1)
          slcs2[i] = slice(1, None)
          # replace from the right
          repmask = np.logical_and(~invalid[slcs1], invalid[slcs2])
          data[slcs1][repmask] = data[slcs2][repmask]
          invalid[slcs1][repmask] = True
          # replace from the left
          repmask = np.logical_and(~invalid[slcs2], invalid[slcs1])
          data[slcs2][repmask] = data[slcs1][repmask]
          invalid[slcs2][repmask] = True
    return data

def fill_inpaint(array,invalid=None,max_iter=5,tol=0.5,kernel_size=1,method='localmean'):
    """Replace NaN elements in an array using an iterative image inpainting algorithm.
    The algorithm is the following:
    1) For each element in the input array, replace it by a weighted average
    of the neighbouring elements which are not NaN themselves. The weights depends
    of the method type. If ``method=localmean`` weight are equal to 1/( (2*kernel_size+1)**2 -1 )
    2) Several iterations are needed if there are adjacent NaN elements.
    If this is the case, information is "spread" from the edges of the missing
    regions iteratively, until the variation is below a certain threshold.
    Parameters
    ----------
    array : 2d np.ndarray
    invalid : mask for array (masked elements True), optional
    an array containing NaN elements that have to be replaced
    max_iter : int
    the number of iterations
    kernel_size : int
    the size of the kernel, default is 1
    method : str
    the method used to replace invalid values. Valid options are
    `localmean`, 'idw'.
    Returns
    -------
    filled : 2d np.ndarray
    a copy of the input array, where NaN elements have been replaced.
    Credits
    -------
    http://stackoverflow.com/a/17125125/512111
    """
    array = np.asarray(array,np.float64)
    if invalid is not None:
       array[np.where(invalid)] = np.nan
 
    # convert masked array to array filled with nans
    array = np.ma.filled(array,np.nan)
    filled = np.empty( [array.shape[0], array.shape[1]], dtype=np.float64)
    kernel = np.empty( (2*kernel_size+1, 2*kernel_size+1), dtype=np.float64 )
    # indices where array is NaN
    inans, jnans = np.nonzero(np.isnan(array))
    # number of NaN elements
    n_nans = len(inans)
    # arrays which contain replaced values to check for convergence
    replaced_new = np.zeros( n_nans, dtype=np.float64)
    replaced_old = np.zeros( n_nans, dtype=np.float64)
    # depending on kernel type, fill kernel array
    if method == 'localmean':
        print('kernel_size', kernel_size)
        for i in xrange(2*kernel_size+1):
            for j in xrange(2*kernel_size+1):
                kernel[i,j] = 1
        print(kernel, 'kernel')
    
    elif method == 'idw':
        kernel = np.array([[0, 0.5, 0.5, 0.5,0],
            [0.5,0.75,0.75,0.75,0.5],
            [0.5,0.75,1,0.75,0.5],
            [0.5,0.75,0.75,0.5,1],
            [0, 0.5, 0.5 ,0.5 ,0]])
        print(kernel, 'kernel')
    
    else:
        raise NotImplementedError('Method not valid. Should be one of [\'localmean\',\'idw\'].')
    
    # fill new array with input elements
    for i in xrange(array.shape[0]):
        for j in xrange(array.shape[1]):
            filled[i,j] = array[i,j]
    
    # make several passes
    # until we reach convergence
    for it in xrange(max_iter):
        print('iteration', it)
        # for each NaN element
        for k in xrange(n_nans):
            i = inans[k]
            j = jnans[k]

            # initialize to zero
            filled[i,j] = 0.0
            n = 0

            # loop over the kernel
            for I in xrange(2*kernel_size+1):
                for J in xrange(2*kernel_size+1):

                    # if we are not out of the boundaries
                    if i+I-kernel_size < array.shape[0] and i+I-kernel_size >= 0:
                        if j+J-kernel_size < array.shape[1] and j+J-kernel_size >= 0:

                            # if the neighbour element is not NaN itself.
                            if filled[i+I-kernel_size, j+J-kernel_size] == filled[i+I-kernel_size, j+J-kernel_size] :

                                # do not sum itself
                                if I-kernel_size != 0 and J-kernel_size != 0:

                                    # convolve kernel with original array
                                    filled[i,j] = filled[i,j] + filled[i+I-kernel_size, j+J-kernel_size]*kernel[I, J]
                                    n = n + 1*kernel[I,J]

            # divide value by effective number of added elements
            if n != 0:
                filled[i,j] = filled[i,j] / n
                replaced_new[k] = filled[i,j]
            else:
                filled[i,j] = np.nan

        # check if mean square difference between values of replaced
        #elements is below a certain tolerance
        print('tolerance', np.mean( (replaced_new-replaced_old)**2 ))
        if np.mean( (replaced_new-replaced_old)**2 ) < tol:
            break
        else:
            for l in xrange(n_nans):
                replaced_old[l] = replaced_new[l]
    
    return filled

def test_fill(s,d,fill=fill_ndimage):

    import matplotlib.pyplot as plt
    # s is size of one dimension, d is the number of dimension
    data = np.arange(s**d).reshape((s,)*d)
    seed = np.zeros(data.shape,dtype=bool)
    seed.flat[np.random.randint(0,seed.size,int(data.size/20**d))] = True

    data_filled = fill(data,-seed)

    # draw (dilated) seeds in black
    data_filled[nd.binary_dilation(seed,iterations=2)] = 0

    plt.figure()
    plt.imshow(np.mod(data_filled,42))
    plt.show()
    
if __name__ == '__main__':

    test_fill(500,2,fill_ndimage)
