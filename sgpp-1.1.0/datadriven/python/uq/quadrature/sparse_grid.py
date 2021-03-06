from pysgpp import (createOperationQuadrature,
                    Linear, LinearL0Boundary, LinearBoundary,
                    Poly, PolyBoundary)
from pysgpp.extensions.datadriven.uq.operations import getBasis
import numpy as np


def getIntegral(grid, level, index):
    # create new grid
    if grid.getType() in [LinearBoundary, LinearL0Boundary]:
        return np.power(2., -max(1, level))
    elif grid.getType() == Linear:
        # # employ 4/3 rule
        # if gp.isLeaf():
        #     q *= 4. / 3.
        return np.power(2., -level)
    elif grid.getType() == Poly:
        return getBasis(grid).getIntegral(level, index)
    elif grid.getType() == PolyBoundary:
        return getBasis(grid).getIntegral(level, index)
    else:
        raise AttributeError('unsupported grid type %s' % grid.getType())


def getIntegralOfBasisFunction(grid, gp):
    quad = 1.
    for d in xrange(gp.dim()):
        quad *= getIntegral(grid, gp.getLevel(d), gp.getIndex(d))
    return quad


def doQuadrature(grid, alpha):
    try:
        return createOperationQuadrature(grid).doQuadrature(alpha)
    except Exception:
        # import ipdb; ipdb.set_trace()
        s = 0.0
        gs = grid.getStorage()

        # set function values for n_alpha
        for i in xrange(gs.size()):
            gp = gs.get(i)

            q = 1.0
            for d in xrange(gs.dim()):
                q *= getIntegral(grid, gp.getLevel(d), gp.getIndex(d))

            s += alpha[i] * q
        return s
