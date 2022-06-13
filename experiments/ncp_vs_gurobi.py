from cProfile import label
import sys
sys.path.insert(1, './src')

import random_instance
import single_level_miqp
import naive_cutting_plane
import parameters as params

import matplotlib.pyplot as plt
import numpy as np

# get the instance
n_node, n_arc, n_com, my_seed = 10,40,10,1234
instance = random_instance.Instance(n_node, n_arc, n_com, my_seed)

# gurobi
formulation = single_level_miqp.SingleLevelMIQP(instance)
formulation.build_model()
formulation.solve_model()
gurobi_obj = formulation.obj
gurobi_time = round(formulation.time,2)

# naive cutting plane
ncp = naive_cutting_plane.CuttingPlane(instance)
ncp.run()
ncp_ls_obj = ncp.ls_obj
ncp_time = ncp.time
ncp_row_util = ncp.row_util

improvements_steps, imp_obj = [], []
non_improvement_steps, non_imp_obj = [], []
for i in range(1, len(ncp_ls_obj)):
    if ncp_ls_obj[i] >= ncp_ls_obj[i-1] + params.tolerance:
        improvements_steps.append(i)
        imp_obj.append(ncp_ls_obj[i])
    else:
        non_improvement_steps.append(i)
        non_imp_obj.append(ncp_ls_obj[i])

imp_step_rate = round((len(improvements_steps) / (len(improvements_steps) + len(non_improvement_steps))) * 100, 2)


# plot the results
n_item = len(ncp_ls_obj)
plt.figure(figsize=(14,7))
gurobi_label = 'Gurobi(' + str(gurobi_time) + ' sec)'
plt.plot(list(range(n_item)), np.ones(n_item)*gurobi_obj, c='black', label=gurobi_label)
ncp_label = 'Naive Cutting Plane(' + str(ncp_time) + ' sec)'
plt.plot(list(range(n_item)), ncp_ls_obj, c='orange', label=ncp_label)
plt.fill_between(list(range(n_item)), np.ones(n_item)*gurobi_obj, ncp_ls_obj, color='tab:brown', alpha=0.2, label='GAP')
plt.scatter(improvements_steps, imp_obj, label='improvement steps in NCP', color="red", marker='^',s=100)
plt.scatter(non_improvement_steps, non_imp_obj, label='non-improvement steps in NCP', color="blue", marker='>',s=100)
text1 =  'row utilization: ' + str(ncp_row_util) + '%'
plt.text(n_item - 11, min(ncp_ls_obj)+10, text1, bbox=dict(facecolor='red', alpha=0.5))
text2 = 'improvement step rate: ' + str(imp_step_rate) + '%'
plt.text(n_item - 11, min(ncp_ls_obj)+13, text2, bbox=dict(facecolor='red', alpha=0.5))
plt.xlabel("Iteration Number", fontsize=12)
plt.ylabel("Objective Function Value", fontsize=12)
plt.xlim([0, n_item])
plt.ylim([min(ncp_ls_obj)-5, max(ncp_ls_obj)+5])
title_ = str(n_node) + ' nodes - ' + str(n_arc) + ' arcs - ' + str(n_com) + ' commodities with seed ' + str(my_seed)
plt.title(title_, fontsize=20)
plt.legend()
plt.show()
