import single_level_miqp
import parameters as params
import matplotlib.pyplot as plt


class CuttingPlane(object):

    def __init__(self, instance):

        self.ins_ = instance
        self.iter_no = 0
        self.cuts = []
        self.ls_obj = []
    
    def run(self):

        while True:
            f = single_level_miqp.SingleLevelMIQP(self.ins_, self.iter_no, self.cuts)
            f.build_model()
            f.solve_model()
            self.ls_obj.append(f.obj)
            pair,value=f.get_most_violated_disjunction()
            print(self.iter_no, pair, f.obj)
            if value <= params.tolerance:  # it means that we get the optimal solution
                break
            self.cuts.append(pair)
            self.iter_no += 1
            f.model.reset(0)
        
    def plot(self):

        n = len(self.ls_obj)
        plt.plot(list(range(n)), self.ls_obj)
        plt.show()
