import sys
sys.path.insert(1, './src')

import random_instance
import single_level_miqp

n_node, n_arc, n_com, my_seed = 24,80,15,1234
instance = random_instance.Instance(n_node, n_arc, n_com, my_seed)
formulation = single_level_miqp.SingleLevelMIQP(instance)
formulation.build_model()
formulation.write_model()
formulation.solve_model()
formulation.summarize_results()

print(formulation.obj)
print('done')
