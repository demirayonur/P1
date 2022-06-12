import sys
sys.path.insert(1, './src')

import random_instance
import naive_cutting_plane

n_node, n_arc, n_com, my_seed = 50,200,30,1234
instance = random_instance.Instance(n_node, n_arc, n_com, my_seed)

ncp = naive_cutting_plane.CuttingPlane(instance)
ncp.run()
ncp.plot()