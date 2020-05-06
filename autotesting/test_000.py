import os
import shutil
import platform


def test_setup():
	if os.path.exists(os.path.join(".", "temp")):
		shutil.rmtree(os.path.join(".", "temp"))
	os.mkdir(os.path.join(".", "temp"))

	
def test_imports():
	import lb_colloids
	from lb_colloids import PSphere
	from lb_colloids import LBImage
	from lb_colloids import LBBC
	from lb_colloids import LB2DModel
	from lb_colloids import Colloids
	if platform.system().lower() == "windows":
		import lb_colloids.LB.LB2Dw
	else:
		import lb_colloids.LB.LB2D


if __name__ == "__main__":
	test_setup()
	test_imports()
