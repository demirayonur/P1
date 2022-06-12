import matplotlib.pyplot as plt
import parameters as params
import networkx as nx
import numpy as np
import random 


class Instance(object):

    def __init__(self,n, m, k, my_seed):

        '''
            Generates a problem instance given the 
            number of 
            n: nodes
            m: arcs
            k: commodities
            for a given seed value for the cases where
            we use randomness. 

            Note: all entities are counted starting from zero
                  UNLIKE real-life problem instances. 
        '''
        
        # store the sizes
        self.n, self.m, self.k = n, m, k

        # fix my seeds
        np.random.seed(my_seed)
        random.seed(my_seed)

        # random graph from network X
        self.graph = nx.gnm_random_graph(n, m, seed=my_seed, directed=True)

        # get incidence matrix
        self.incidence_matrix = nx.incidence_matrix(self.graph, oriented=True).todense()

        # random commodities
        demands = np.random.randint(params.min_demand, params.max_demand, k)  # random demands for commodities that will be generated
        self.commodities = {}  # we will store our commoditiies here:  commodity key --> (from_node, to_node, demand)
        counter = 0  # keeps track of how many commodities generated until that time.
        o_d_pairs = []  # we will store OD pairs here.
        while counter < k:  
            from_node = random.choice(list(range(n)))  # from which node
            to_node = random.choice(list(range(n)))    # to which node
            candidate_edge = [from_node, to_node]      # our candidate edge
            if from_node == to_node:   # no cycles are allowed
                continue
            if candidate_edge in o_d_pairs:  # repetation in OD pairs is not allowed
                continue
            if not nx.has_path(self.graph, from_node, to_node):  # there has to be at least one path between origin to destination
                continue
            # if there is no problem from the above issues:
            o_d_pairs.append(candidate_edge)  # accept the candidate edge as an OD pair.
            self.commodities[counter] = (from_node, to_node, demands[counter])  # create the commodity
            counter += 1  # increase the counter.

        # get random parameters on the arcs
        self.alpha = [round(i,2) for i in np.random.uniform(params.min_alpha,params.max_alpha, m)]
        self.beta = [round(i,2) for i in np.random.uniform(params.min_beta,params.max_beta, m)]

        # toll-control
        n_allowed_toll = int(m * params.rate_toll)
        list_of_allowed_tolls = random.sample(list(range(m)), n_allowed_toll)
        self.theta = [params.max_toll if a in list_of_allowed_tolls else 0 for a in range(m)]

