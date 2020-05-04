from lb_colloids import LBImage
from lb_colloids import LB2DModel
import matplotlib.pyplot as plt
from lb_colloids import cIO
from lb_colloids import ColloidModel
import os
import time


def test_setup():
	if os.path.exists(os.path.join(".", "temp")):
		os.remove(os.path.join(".", "temp"))
	os.mkdir(os.path.join(".", "temp"))

	
def test_imports():
	import lb_colloids
	from lb_colloids import PSphere
	from lb_colloids import LBImage
	from lb_colloids import LB2DModel
	from lb_colloids import Colloids
	import lb_colloids.LB.LB2D


if __name__ == "__main__":
	test_setup()
	test_imports()
