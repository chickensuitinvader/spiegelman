"""
Spatial Variant of Spiegelman's Monster Simulation
Allows templates to move between simulations based on their composition
Author: David Wu
Universty of Auckland
"""

import spiegelman as spg
import run_simulation as rsim
from operator import itemgetter
import pylab

# spatial simulation object
# object contains:
#   - instances: dict of index to SpSim objects
#   - parameters: simulation parameters
#   - history: recorded information about movements
class SptSim(object):
    # initialisation method, sim parameter can be used to specify an object
    # or output file to copy asimulation from
    def __init__(self, sim = None):
        spg.importlib.reload(spg.imports)
        if sim == None:
            self.instances = dict()
            self.parameters = rsim.pSet(spg.imports.parameters)
            self.parameters.set('Cycles', 1)                     # each SpSim runs for one cycle before shuffle process
            self.__set_size__(self.parameters.get('Size'))
            self.history = list()
        elif isinstance(sim, str):
            self.read_from(sim)
        else:
            try:
                self.instances = {x:spg.SpSim(sim.instances[x]) for x in sim.instances}
                self.parameters = sim.parameters
                self.history = sim.history
            except:
                raise(TypeError)
    
    # initialisation helper function, generates instances with indices
    def __set_size__(self, size):
        if not isinstance(size,tuple):
            raise(TypeError)
        ps = list()
        for n in range(size[0]*size[1]): 
            ps.append(rsim.pSet(self.parameters))
            ps[n].set('Instance', n)
            ps[n].set('Changes', ps[n].get('Changes')[n])
            pool = dict(ps[n].get('InitialPool'))
            pool[ps[n].get('Changes')[1]] *= self.parameters.get('ChangeValue')
            ps[n].set('InitialPool', pool)
            self.instances[n] = spg.SpSim(ps[n])
    
    def what(self,f=None):
        if f == None:
            [print(x, ':', self.parameters.list()[x]) for x in self.parameters.list()]
        else:
            print(self.parameters.list(), file = f)
            
    # method to run the simulations
    def run(self, printing = -1):
        start = spg.time.time()
        for i in range(self.parameters.get('Epochs')):
            for sp in self.instances:
                self.instances[sp].run(0)
            self.shuffle(printing)
        if printing:
            print('Elapsed Time:', round(spg.time.time()-start,3), 's')
    
    # method to cause templates from the simulations to move
    # current implementation:
    #       Preference Method:
    #           Move templates by their most common base to the simulation
    #           which ended with the most of that base in the previous epoch
    def shuffle(self,printing = -1):
        if (printing > 0):
            print('Shuffling')
        if self.parameters.get('ShuffleType') == 'Preference':
            ranking = self.getRanking()
            self.history.append(ranking)
            for sp in self.instances:
                for s in self.instances[sp].templates:
                    count = set((x,s.count(x)) \
                                for x in self.parameters.get('InitialPool'))
                    if ((100*spg.random.random()) < self.parameters.get('ShufflePercent')):
                        mover = max(count, key=itemgetter(1))[0]   # most influential
                        target = ranking[mover][0][0]
                        self.instances[target].templates.append(s)
                        self.instances[sp].templates.remove(s)    
        if self.parameters.get('ShuffleType') == 'ComplexPref':
            ranking = self.getRanking()
            self.history.append(ranking)
            for sp in self.instances:
                inst = self.instances[sp].parameters['Instance']
                for s in self.instances[sp].templates:
                    count = {x:s.count(x) \
                                for x in self.parameters.get('InitialPool')}
                    if ((100*spg.random.random()) < self.parameters.get('ShufflePercent')):
                        moveVector = {i:sum([count[b]*(ranking[b][inst][1] - ranking[b][i][1]) \
                                             for b in self.parameters.get('InitialPool')]) \
                                      for i in self.instances}
                        moveTo = max(moveVector, key = moveVector.get)
                        self.instances[moveTo].templates.append(s)
                        self.instances[sp].templates.remove(s)
                        
                                                
            
    def getRanking(self):
        ranking = dict()
        for x in self.parameters.get('InitialPool').keys():
            ranking[x] = [(sp,self.instances[sp].pool.quantities[x]) \
                            for sp in self.instances]
            ranking[x].sort(key = itemgetter(1), reverse = True)
        return ranking
        
    # method to print to console relevant information TODO    
    def toPrint(self):
        for sp in self.instances:
            self.instances[sp].toPrint()
    
    # convenience method
    def plotting(self, key = None, instance = None):
        spg.plt.close()
        if (instance != None):
                self.instances[instance].plotting(key)
        else:
            f,ax = spg.plt.subplots(self.parameters.get('Size')[0], \
                                    self.parameters.get('Size')[1])
            for inst in self.instances:
                axis = ax.ravel()[inst]
                if key == None:
                    print('For Instance', inst)
                if key[0] == '3':
                    im = self.instances[inst].plotting3(key[1:],axis)
                    f.colorbar(mappable = im, ax=axis)
                else:
                    self.instances[inst].plotting(key,axis)
                spg.plt.suptitle(axis.get_title())
                axis.set_title('Instance ' + str(inst))
        
    
    # method to write simulation data to external data files
    # will export to a folder, and each SpSim would be its own file
    def export_to(self, folder = None):
        if folder == None:
            folder =  'Sim_' + spg.time.strftime('%d%m%y_%H%M%S')
        rsim.os.mkdir(folder)
        for n in self.instances:
            file = folder + '/' + str(n) + '.SIMHIST'
            self.instances[n].export_to(file)
        historyFile = open(folder + '/poolhistory.SPTHIST', 'w')
        print(self.history, file = historyFile)
        historyFile.close()
        paramFile = open(folder + '/parameters.SPT', 'w')
        self.what(paramFile)
        paramFile.close()
        return folder
    
    # method to extract data from a folder of external data files
    def read_from(self, folder):
        n = 0
        self.instances = dict()
        for fileName in rsim.os.listdir(folder):
            if fileName[-8:] == '.SIMHIST':
                self.instances[n] = spg.SpSim(folder+ '/'+ fileName)
                n += 1
        if (n > 0):
            try:
                f = open(folder+ '/poolhistory.SPTHIST','r')
                self.history = spg.ast.literal_eval(f.readline())
                f.close()
                f = open(folder+'/parameters.SPT','r')
                self.parameters = rsim.pSet(spg.ast.literal_eval(f.readline()))                    
                f.close()
            except:
                print('Incomplete Information')
        else:
            raise Exception('No Files Found in ' + folder)
        
# end of class SptSim   

# convenience method to run and write data files for a simulation, and then
# return simulation        
def go(inFile = None, printing = -1, outFile = None):
    try:
        s = SptSim(inFile)
        s.run(printing)
        f = s.export_to(outFile)
        autoPlot(s,destination=f)
    except Exception as ex:
        print(ex)
    finally:
        return s

#convenience method
def autoPlot(sim, graphs = ('a','n','p','u%','3h'), destination = '.', ret = False):
    if isinstance(sim,str):
        ret = True
        destination = sim
        s = SptSim(sim)    
    else:
        t = type(sim)
        s = t(sim)
    for g in graphs:
        s.plotting(g)
        figManager = spg.plt.get_current_fig_manager()
        figManager.window.showMaximized()
        pylab.savefig(destination + '/Graph_' + g + '.png')
    try:
        f=open(destination+'/paramters.PREC')
        s.what(f)
        f.close()
    except:
        pass
    if ret:
        return s

