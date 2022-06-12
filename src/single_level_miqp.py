import parameters as params
from gurobipy import GRB
import gurobipy as gp
import pandas as pd
import numpy as np
import time 


class SingleLevelMIQP(object):

    def __init__(self, instance):
        
        self.ins_ = instance
        self.model = gp.Model("SingleLevelMIQP")

        if params.flag == False:
            self.model.setParam('OutputFlag', 0)

        # sets
        self.nodes, self.arcs, self.pairs = range(instance.n), range(instance.m), range(instance.k)

        # continuous decision variables
        self.tau = self.model.addVars(self.arcs, name="tau", lb=np.zeros(instance.m), ub=instance.theta)
        self.f = self.model.addVars(self.arcs, name="f", lb=np.zeros(instance.m))
        self.X = self.model.addVars(self.arcs, self.pairs, name='X', lb=0)
        self.Lambda = self.model.addVars(self.arcs, self.pairs, name='Lambda', lb=0)
        self.Mu = self.model.addVars(self.nodes, self.pairs, name='Mu')  

        # binary decision variables
        self.Xi = self.model.addVars(self.arcs, self.pairs, vtype=GRB.BINARY, name="Xi")

        # fill after solving it
        self.status, self.obj, self.gap, self.time = None, None, None, None

    def add_objective_function(self):

        congestion = gp.quicksum(self.ins_.alpha[a] * (self.f[a] ** 2) + self.ins_.beta[a] * self.f[a] for a in self.arcs)
        self.model.setObjective(congestion, sense=GRB.MINIMIZE)

    def add_f_X_relation(self):
        
        self.model.addConstrs((self.X.sum(a, '*') == self.f[a] for a in self.arcs), "f_X_")
    
    def add_flow_balance_constraints(self):

        counter = 0
        for j in self.ins_.commodities:
            origin, destination, demand = self.ins_.commodities[j]
            for i in self.nodes:
                temp_name = "balance_" + str(counter)
                if i == origin:
                    self.model.addConstr(sum(self.ins_.incidence_matrix[i, a] * self.X[a, j] for a in self.arcs)==-demand, temp_name)
                elif i == destination:
                    self.model.addConstr(sum(self.ins_.incidence_matrix[i, a] * self.X[a, j] for a in self.arcs)==demand, temp_name)
                else:
                    self.model.addConstr(sum(self.ins_.incidence_matrix[i, a] * self.X[a, j] for a in self.arcs)==0, temp_name)
                counter += 1
    
    def add_disjunctions(self):
        
        counter = 0
        for a in self.arcs:
            for j in self.pairs:
                temp_name_x = "disjunction_x_" + str(counter)
                temp_name_lambda = "disjunction_lambda_" + str(counter)
                self.model.addConstr(self.X[a, j] <= 1000000 * self.Xi[a,j], temp_name_x)
                self.model.addConstr(self.Lambda[a,j] <= 1000000 * (1 - self.Xi[a,j]), temp_name_lambda)
                counter +=1

    def add_equilibrium_conditions(self):
        
        counter = 0
        for a in self.arcs:
            for j in self.pairs:
                temp_name = 'equilibrium_' + str(counter)
                free_sum = gp.quicksum(self.Mu[i,j] * self.ins_.incidence_matrix[i,a] for i in self.nodes)
                self.model.addConstr(self.ins_.alpha[a] * self.f[a] + self.ins_.beta[a] + self.tau[a] - self.Lambda[a,j] + free_sum == 0, temp_name)
                counter += 1 
    
    def build_model(self):

        self.add_objective_function()
        self.add_f_X_relation()
        self.add_flow_balance_constraints()
        self.add_disjunctions()
        self.add_equilibrium_conditions()

    def write_model(self):
        
        # do not call this method before building it

        self.model.write("SingleLevelMIQP.lp")
    
    def solve_model(self):

        # do not call this method before building it

        time_start = time.time()
        self.model.optimize()
        time_final = time.time()
        solution_time = time_final - time_start

        self.status, self.obj, self.gap, self.time =  self.model.Status, self.model.ObjVal, self.model.MIPGap, solution_time
    
    def summarize_results(self):
        
        # do not call this method before solving it
        if self.status == 2 or self.status == 9:  # optimal or feasible under time limit
            
            # sheet 0 (solution statistics)
            df_stat = pd.DataFrame()
            df_stat['description'] = ['solver', 'status', 'obj', 'gap', 'time']
            df_stat['value'] = ['gurobi', self.status, self.obj, self.gap, self.time]

             # sheet 1 (arc parameters + variables)
            df_arcs = pd.DataFrame() 
            df_arcs['from'] = [arc[0] for arc in self.ins_.graph.edges]
            df_arcs['to'] = [arc[1] for arc in self.ins_.graph.edges]
            df_arcs['alpha'] = self.ins_.alpha
            df_arcs['beta'] = self.ins_.beta
            df_arcs['is_tollable'] = [1 if i>0 else 0 for i in self.ins_.theta]
            df_arcs['toll'] = [self.tau[a].X for a in self.arcs]
            df_arcs['flow'] = [self.f[a].X for a in self.arcs]
            
            # sheet 2 (X)
            df_X = pd.DataFrame(np.array([self.X[a,j].X for a in self.arcs for j in self.pairs]).reshape((self.ins_.m, self.ins_.k)))

            # sheet 3 (Lambda)
            df_Lambda = pd.DataFrame(np.array([self.Lambda[a,j].X for a in self.arcs for j in self.pairs]).reshape((self.ins_.m, self.ins_.k)))

            # sheet 4 (Xi)
            df_Xi = pd.DataFrame(np.array([1 if self.Xi[a,j].X>=params.tolerance else 0 for a in self.arcs for j in self.pairs]).reshape((self.ins_.m, self.ins_.k)))

            # sheet 5 (Mu)
            df_Mu = pd.DataFrame(np.array([self.Mu[i,j].X for i in self.nodes for j in self.pairs]).reshape((self.ins_.n, self.ins_.k)))

            # Create a Pandas Excel writer using XlsxWriter as the engine.
            file_path = './results_miqp_original.xlsx'
            writer = pd.ExcelWriter(file_path, engine='xlsxwriter')

            # Write each dataframe to a different worksheet. 
            df_stat.to_excel(writer, index=False, sheet_name='solution_stats')
            df_arcs.to_excel(writer, index=False, sheet_name='arcs')
            df_X.to_excel(writer, index=False, sheet_name='X')
            df_Lambda.to_excel(writer, index=False, sheet_name='Lambda')
            df_Xi.to_excel(writer, index=False, sheet_name='Xi')
            df_Mu.to_excel(writer, index=False, sheet_name='Mu')

            # Close the Pandas Excel writer and output the Excel file.
            writer.save()
