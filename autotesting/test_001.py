import os
from lb_colloids import PSphere
from lb_colloids import LBImage
import matplotlib.image as mpimg
import numpy as np


ws = os.path.abspath(os.path.dirname(__file__))
data = os.path.join(ws, "..", "data")
examples = os.path.join(ws, "..", "examples", "data")
outpth = os.path.join(ws, "temp")


def test_psphere():
	grain_radius = 20
	porosity = 0.40
	dimension = 200
	sensitivity = 0.01

	psp = PSphere(radius=grain_radius, porosity=porosity,
				  dimension=dimension, sensitivity=sensitivity)

	matrix = psp.get_matrix()
	mpor = PSphere.static_porosity(matrix)

	if (porosity - sensitivity) <= mpor <= (porosity + sensitivity):
		pass
	else:
		raise AssertionError("Porosity out of defined tolerance")

	if not psp.percolates:
		raise AssertionError("Porous media does not percolate")

	porosity = 0.25

	psp = PSphere(radius=grain_radius, porosity=porosity,
				  dimension=dimension, sensitivity=sensitivity)

	mpor = psp.matrix_porosity

	if (porosity - sensitivity) <= mpor <= (porosity + sensitivity):
		pass
	else:
		raise AssertionError("Porosity out of defined tolerance")

	if not psp.percolates:
		raise AssertionError("Porous media does not percolate")

	psp.save_image(os.path.join(outpth, "psphere_test.png"))

	img = mpimg.imread(os.path.join(outpth, "psphere_test.png"))
	img_matrix = np.array(img[:, :, 0], dtype=bool)

	img_por = PSphere.static_porosity(img_matrix)

	if abs(mpor - img_por) > sensitivity:
		raise AssertionError("Porous media not written to disk correctly")


def test_lb_image():
	image_name = "test_image.png"
	img = LBImage.Images(os.path.join(examples, image_name))

	pmin = np.min(img.arr)
	pmax = np.min(img.arr)

	if 0 <= pmin <= 255:
		pass
	else:
		raise AssertionError("Image could not be converted to 8 bit")

	if 0 <= pmax <= 255:
		pass
	else:
		raise AssertionError("Image could not be converted to 8 bit")


def test_lb_boundary_condition():
	image_name = "test_image.png"
	img = LBImage.Images(os.path.join(examples, image_name))

	unique = np.unique(img.arr)

	fluidvx = [unique[0], ]
	solidvx = [unique[1], ]
	nlayers = 5

	bc = LBImage.BoundaryCondition(img.arr, fluidvx, solidvx,
								   nlayers)

	if not np.allclose(bc.grey_values, unique):
		raise AssertionError("Grey values not properly passed to "
							 "BoundaryCondition class")

	if bc.porosity != 0.4044:
		raise AssertionError("Porosity calculation has failed")

	if bc.binarized.shape != (210, 202):
		raise AssertionError("Boundary conditions not added in the correct "
							 "dimension")

	if not all(bc.binarized[:, 0]):
		raise AssertionError("No flow boundary not added properly")

	if not all(bc.binarized[:, -1]):
		raise AssertionError("No flow boundary not added properly")


def test_lb_model_fortran():
	pass


def test_lb_model_python():
	pass

	
def test_colloid_model():
	pass
	

if __name__ == "__main__":
	test_psphere()
	test_lb_image()
	test_lb_boundary_condition()