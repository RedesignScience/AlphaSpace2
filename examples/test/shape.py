"""Testing script for shape functions in OESHAPE"""

from openeye.oeshape import *
from openeye.oechem import *
from openeye.oegrid import *
import numpy as np

"""
Initializing shape query
"""

shape_query = OEShapeQuery()

shape_coords = np.random.randn(10, 3) * 10

for x, y, z in shape_coords:
    shape_gaussian = OEGaussian()
    shape_gaussian.SetCenter(OEFloatArray((x, y, z)))
    shape_query.AddShapeGaussian(shape_gaussian)
    print(x, shape_gaussian.GetX())

shifted_shape_coords = shape_coords * 10
shifted_shape_query = OEShapeQuery()

for x, y, z in shifted_shape_coords:
    shape_gaussian = OEGaussian()

    shape_gaussian.SetCenter(OEFloatArray((x, y, z)))
    shifted_shape_query.AddShapeGaussian(shape_gaussian)

result = OEBestOverlayResults()
overlay = OEOverlay()

overlay.SetupRef(shape_query)
