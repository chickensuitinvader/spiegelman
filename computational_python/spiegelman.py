"""
Simulation v4
Implementation of Spiegelman's Monster Experiment
Author: David Wu
University of Auckland
"""

import random, importlib, time, ast
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection
from matplotlib import cm as CM
from operator import itemgetter
import parameters as imports
import cull_function

#print('Current Known Issues: Does not accept negative cull_funtion')
    
# simulation object
# Object contains:
#    parameters - dictionary of simulation parameters
#    replicators - list of Replicator objects
#    templates - list of template strings
#    history - dict of recorded data
#    pool - Pool object, records number of monomers remiaing
class SpSim(object):    
    # object can be intialised using the sim input:
    # Type: sim - SpSim, str
    def __init__(self, sim = None):
        self.history = {'Number':[], 'Average':[], 'Lengths':[], 'Pool':[], 'Uniques':[], 'EarlyQuit':[], 'Progress':[]}
        self.parameters = {}
        self.templates = []
        if sim is None:
            importlib.reload(imports)
            self.parameters = dict(imports.parameters)
            self.templates = [self.makeTemplate() for n in range(self.parameters['InitialTemplates'])]
        elif isinstance(sim, str):
            self.read_from(sim)
        else:
            try:
                self.parameters = dict(sim.parameters)
                self.templates = list(sim.templates)
                self.history = dict(sim.history)
            except Exception as ex:
                if self.templates != []:
                    #print('No History Found... Continuing')
                    pass
                elif self.parameters != {}:
                    #print('No Templates Found... Creating Templates and Continuing' )
                    self.templates = [self.makeTemplate() for n in range(self.parameters['InitialTemplates'])]
                    pass
                else:
                    print('Error in Assimilating Data')
                    print(ex)
                    raise(Exception)
        self.replicators = [Replicator(self.parameters['Pairings'],\
                                       self.parameters['PointMutations'], \
                                       self.parameters['BlockMutations']) \
                                       for n in range(int(self.parameters['Replicators']))]
        self.pool = Pool(self.parameters)
        
    # method for making strings of random length from a language            
    def makeTemplate(self):
        length = int(random.normalvariate(self.parameters['SeedLength'][0],self.parameters['SeedLength'][1]))
        return ''.join(random.choice(list(self.parameters['Pairings'].values())) for n in range(length))

    # method for recording data from current simulation state
    # current recorded data:
    #   - Number of Templates
    #   - List of Lengths of Templates
    #   - Median Length of Templates
    #   - Dictionary of Pool Quantities
    #   - Number of Unique Templates
    def addHistory(self):
        self.history['Number'].append(len(self.templates))
        self.history['Lengths'].append([len(n) for n in self.templates])
        self.history['Lengths'][-1].sort()
        self.history['Average'].append(self.history['Lengths'][-1][int(self.history['Number'][-1]/2)])
        self.history['Uniques'].append(len(set(self.templates)))
        self.history['Pool'].append(self.pool.quantities)
    
    def addProgress(self):
        j = [len(n) for n in self.templates]
        j.sort()
        self.history['Progress'][-1].append(j)

    # method to print information about current simulation state
    # Type: iteration - int
    def toPrint(self, iteration = -1):
        print('========================\nParameters: \n========================')
        for n in self.parameters:
            print('\t', n,' : ',self.parameters[n])
        print('========================\nGenome Length Distribution: \n========================')
        print('Unique Genomes:',self.history['Uniques'][iteration])
        dist = {n: self.history['Lengths'][iteration].count(n) for n in self.history['Lengths'][iteration]}
        for k in dist:
            print(k,': ', dist[k], '|\t ', 'x'*dist[k])    
        print('========================\nPool: \n========================')
        print(self.history['Pool'][iteration])
        if iteration != -1:
            return
        print('========================\nReplicators: \n========================')
        i = 0
        for r in self.replicators:
            print('Replicator ', i)
            i += 1
            r.toPrint()
        print('========================\nTo plot, use plotting()\n========================')

    # method to plot aspects of recorded data
    # Type: specifics - str
    def plotting(self, specifics=None, ax = None):
        if specifics == None:
            print('Interactive Mode\n========================\nGraphs Available:')
            print('a: Median Length over Time')
            print('h: Distribution of Genome Lengths')
            print('n: Number of Genomes Strands over Time')
            print('p: Number of Free Nucleotides Remaining at the End of an Epoch')
            print('u: Uniques Genome Species over Time')
            print('u%: Genome Species as Percentage of Total Number of Genome Strands over Time')
            print('q: Iteration where the Pool Empties (if any)')
            print('l: Distribution of Lengths that pass to the next epoch')
            print('#: 3D Snapshot Views of Data')
            specifics = input('Please enter your selection, or enter nothing to quit: ')
            if specifics == '':
                print('Quitting')
                return
            elif specifics == '#':
                self.plotting3(None, ax)
                return
            elif specifics in 'hp':
                totalIt = len(self.history['Number'])
                print('Epochs Thus Far: ', totalIt)
                sub = input('Please enter the epoch number: ')
                specifics += sub
        
        if ax == None:
            plt.close()
            ax = plt.subplot(1,1,1)
        
        if specifics == 'a':
            ax.plot(self.history['Average'])
            ax.set_title('Median Genome Length in Each Epoch')
            ax.set_xlabel('Epoch')
            ax.set_ylabel('Genome Length')
        elif specifics.startswith('h'):
            try:
                ax.hist(self.history['Lengths'][int(specifics[1:])])
                ax.set_title('Genome Length Distribution in Epoch ' + str(specifics[1:]))
                ax.set_xlabel('Genome Length')
                ax.set_ylabel('Frequency')
            except IndexError:
                print('Unknown Iteration')
                raise(IndexError)
            except ValueError:
                print('Unknown Iteration')
                raise(ValueError)
        elif specifics.startswith('l'):
            try:
                ax.hist(self.history['Progress'][int(specifics[1:])][-1])
                ax.set_title('Passing Genome Length Distribution in Epoch ' + str(specifics[1:]))
                ax.set_xlabel('Genome Length')
                ax.set_ylabel('Frequency')
            except IndexError:
                print('Unknown Iteration')
                raise(IndexError)
            except ValueError:
                print('Unknown Iteration')
                raise(ValueError)
        elif specifics == 'n':
            ax.plot(self.history['Number'])
            ax.set_title('Number of Genomes in Each Epoch')
            ax.set_ylabel('Number of Genomes')
            ax.set_xlabel('Epoch')
        elif specifics == 'p':
            ax.hold('on')
            for k in self.history['Pool'][-1]:
                ax.plot([d[k] for d in self.history['Pool']])
            ax.set_title('Number of Free Nucleotides Left at the End of an Epoch')
            ax.set_ylabel('Number of Free Nucleotides')
            ax.set_xlabel('Epoch')
            ax.legend([k for k in self.history['Pool'][-1]])
        elif specifics == 'u':
            ax.plot(self.history['Uniques'])
            ax.set_title('Number of Unique Genome Species in Each Epoch')
            ax.set_xlabel('Epoch')
            ax.set_ylabel('Number of Genome Species')
        elif specifics == 'u%':
            cyc = len(self.history['Number'])
            ax.plot([self.history['Uniques'][n]/self.history['Number'][n] \
                      for n in range(cyc)])
            ax.set_title('Unique Genome Species Percentage')
            ax.set_xlabel('Epoch')
            ax.set_ylabel('Species (%)')
        elif specifics == 'q':
            x = [i[0] for i in self.history['EarlyQuit']]
            if sum(x) == 0:
                x = range(len(x))
            y = [i[1] for i in self.history['EarlyQuit']]
            ax.plot(x,y)
            ax.set_title('Iteration In Which the Pool Empties')
            ax.set_xlabel('Epoch')
            ax.set_ylabel('Iteration')
        else:
            print('Unknown Graph Type')

                
    def plotting3(self, args = None, ax = None, heat = True, ite = -1):
        if args == None:
            print('Interactive Mode\n========================\nGraphs Available:')
            print('h: Distribution of Lengths')
            print('#: Specific Iteration Plots')
            args = input('Please enter your selection, or enter nothing to quit:  ')
            if args == '':
                print('Quitting')
                return
            elif args == '#':
                self.plotting()
                return
        if args == 'h':
            xLabel = 'Genome Length'
            zLabel = 'Frequency'
            bins = 200
            binseq = np.linspace(min(min(self.history['Lengths'])),max(max(self.history['Lengths'])),bins)
            binseq = np.insert(binseq,0,min(binseq)-1)
            binseq = np.insert(binseq,0,0)
            processed = [np.histogram(x,bins=binseq) for x in self.history['Lengths']]
            dataX = [(binseq[i] + binseq[i+1])/2 if i not in [0,1]  else binseq[i] for i in range(bins+1)]
            dataY = [n[0] for n in processed] 
        elif args == 'p':
            xLabel = 'Genome Length'
            zLabel = 'Frequency'
            bins = 200
            llim = min(self.history['Progress'][ite])[0]
            ulim = max(self.history['Progress'][ite],key=itemgetter(-1))[-1]
            binseq = np.linspace(llim,ulim,bins)
            binseq = np.insert(binseq,0,min(binseq)-1)
            binseq = np.insert(binseq,0,0)
            processed = [np.histogram(x,bins=binseq) for x in self.history['Progress'][ite]]
            dataX = [(binseq[i] + binseq[i+1])/2 if i not in [0,1]  else binseq[i] for i in range(bins+1)]
            dataY = [n[0] for n in processed] 
            dataY.insert(-1,[0]*len(dataX))
            [dataY.append(dataY[-1]) for _ in range(10)]
        else:
            print('Unknown Graph Type')
            return
            
        zs = range(len(dataY))
        if heat:
            try:
                im = self.heatDist(dataX,dataY,zs,ax)
            except Exception as ex:
                im = None
                raise ex
            finally:
                return im
            
        vert = [list(zip(dataX, dY)) for dY in dataY]
            
        if ax == None:
            fig = plt.figure()
            ax = Axes3D(fig)
        
        poly = PolyCollection(vert)
        poly.set_alpha(0.7)
        ax.add_collection3d(poly, zs=zs, zdir='y')
        
        ax.set_xlabel(xLabel)
        ax.set_xlim3d(min(dataX),max(dataX))
        #ax.set_xscale('log')
        ax.set_ylabel('Epoch')
        ax.set_ylim3d(-1,max(zs)+1)
        ax.set_zlabel(zLabel)
        ax.set_zlim3d(0,max([k for n in dataY for k in n]))
        #ax.set_zscale('log')
        
        plt.show()

    #method for heatmap style plot
    # requires parent function to provide no axes, or to input own color bar
    def heatDist(self, dataX, dataY, epochs,ax):
        #process dataY to a Probability
        
        if ax == None:
            plt.close()
            fig, ax = plt.subplots(1,1)
            ax0 = False
        else:
            ax0 = True
        
        im = ax.contourf(dataX, epochs, dataY, 200, cmap = CM.spectral)
        ax.set_title('Distribution of Genome Lengths')
        ax.set_xlabel('Genome Length')
        ax.set_ylabel('Epoch/Iteration')
        ax.set_xlim(0,max(dataX))  
        
        if ax0 == False:
            fig.colorbar(im)
            plt.show()
        else:
            return im
            
    # method to run the simulation with the parameters
    # strategy:
    #   - Start timer and record an iteration of history
    #   - In every iteration:
    #       - Replace the pool
    #       - Perform the Iteration Replications
    #       - Record an Iteration of History
    #       - Cull the Population
    # Type: printing - int
    def run(self, printing=-1, cull = None, progress = False):
        if cull == None:
            importlib.reload(cull_function)
            self.parameters['CullFunction'] = cull_function.decision
            cull = cull_function.cull_function
        start = time.time()
        if (printing > 0):
            print('Simulation Start\n========================')
        for iterations in range(int(self.parameters['Cycles'])):
            if progress:
                self.history['Progress'].append([])
            self.pool.initialise()
            self.doIteration(iterations, progress)
            self.addHistory()
            if (printing == 1):
                print('Iteration: ', iterations)
                self.toPrint()
            self.transfer(cull)
            if progress:
                self.addProgress()
        if (printing == 2):
            self.plotting3('h')
        if (printing == -1):
            print('Elapsed Time: ', round(time.time()-start,3), 's')

    # method to complete one iteration of the simulation
    # possible replicator actions per replication:
    #   - Find a template to join to
    #   - Perform a single base replication
    #   - Release copy it has made, and move to another template
    def doIteration(self, iteration, progress):
        empty = False
        for replications in range(int(self.parameters['MaxReplications'])):
            if (empty or self.pool.isLow()):
                self.history['EarlyQuit'].append((iteration, replications))
                break
            for replicator in self.replicators:
                replicator.timer -= 1
                if (self.pool.isEmpty):
                    empty = True
                    break
                elif replicator.stopped():
                    replicator.move(self.templates)
                    replicator.replicate(self.pool)
                elif not replicator.working():
                    copy = replicator.extract()
                    if self.qualifies(copy):
                        self.templates.append(copy)
                        if progress:
                            self.addProgress()
                    replicator.move(self.templates)
                    replicator.replicate(self.pool)
                    

    # method to determine if a string will be added to the pool of templates
    # strategy: 
    #   - check if the string is longer than the minimum length
    # Type: template - str
    def qualifies(self,template):
        return (len(template) >= self.parameters['MinLength'])

    # method to cull the population between iterations of the simulation
    # strategy:
    #   - if transfer percent is above 50%:
    #       - remove templates from template list
    #   - if transfer percent is below 50%:
    #       - construct temp list
    #       - add templates from template list to temp list
    #       - set template list to be temp list
    def transfer(self, cull):
        scores = [cull(x) for x in self.templates]
        newList = np.random.choice(self.templates, \
                                   int(len(self.templates)*0.01*self.parameters['TransferPercent']), \
                                   p=[x/sum(scores) for x in scores], \
                                   replace=False)
                   
        self.templates = list(newList)
#        if self.parameters['TransferPercent'] >= 50:
#            for n in range(int((1 - self.parameters['TransferPercent']/100) * len(self.templates))):
#                self.templates.remove(random.choice(self.templates))
#        else:
#            temp = list()
#            for n in range(int((self.parameters['TransferPercent']/100) * len(self.templates))):
#                temp.append(random.choice(self.templates))
#            self.templates.clear()
#            self.templates.extend(temp)
        for replicator in self.replicators:
            replicator.release()

    # method to record history onto a text file
    # Type: fileName - str
    def export_to(self,fileName):
        f = open(fileName,'w')
        for cat in self.history:
            print(':=History_Category=:', cat, file = f)
            if cat == 'Progress':
                for n in self.history[cat]:
                    print(n, file = f)
            else:
                print(self.history[cat], file = f)
        print('Templates',[n for n in self.templates], file = f)
        print('Parameters', self.parameters, file = f)
        f.close()
    
    # method to read recorded histroy from file
    # Type: fileName - str
    def read_from(self,fileName):
        f = open(fileName,'r')
        cat = str()
        for line in f:
            if line.startswith(':=History_Category=:'):
                cat = line[21:-1]
                if not self.history.__contains__(cat):
                    self.history[cat] = []
            elif line.startswith('Templates'):
                self.templates = ast.literal_eval(line[10:])
            elif line.startswith('Parameters'):
                self.parameters = ast.literal_eval(line[11:])
            elif not(cat is str()):
                if cat == 'Progress':
                    self.history[cat].append(eval(line))
                else:
                    self.history[cat] = ast.literal_eval(line)
            else:
                print('Format Error')
                raise(TypeError)
        f.close()

#end of class SpSim        
        
# pool object 
# Object contains:
#   initials - parameter corresponding to a full pool
#   quantities - dictionary containg counts of current number of monomers
#   elements - list of possible monomers in pool               
class Pool(object):
    
    # intialise object with parameters
    # Type: parameters - dict
    def __init__(self,parameters):
        self.initials = dict(parameters['InitialPool'])
        self.quantities = dict()
        self.elements = list()
        try:
            self.emptiness = float(parameters['EmptyPool'])
            self.full = sum(self.initials.values())
            self.current = float('inf')
        except:
            print('Legacy Type, No Low Pool Checks, DO NOT RUN')

    # method to reset the quantity of each base in the pool to the intial
    def initialise(self):
        self.quantities = dict(self.initials)
        self.elements = [k for k in self.quantities.keys()]
        self.current = self.full
        self.isEmpty = False
                         
    # method to determine if the pool is low
    # the pool is low if the amount of bases in the pool is below a proportion
    # of the initial number of bases
    def isLow(self):                      
        return (self.current/self.full < self.emptiness)

    # method to reduce the number of bases in the pool
    # Type: mType - str
    def deplete(self,mType):
        for s in mType:
            self.quantities[s] -= 1
            self.current -= 1
        if any(0 >= x for x in self.quantities.values()):
            self.isEmpty = True

    # method of printing current pool information
    def toPrint(self):
        for x in self.elements:
            print(x,' : ', self.quantities[x],'; ',end='\t')
        print()

# end of class Pool        
        
# replicator object
# Object contains:
#   point - parameters corresponding to point mutation types and rates
#   block - parameters corresponding to block mutation types and rates
#   template - template string that is being copied
#   copy - template string being formed
#   position - position of Replicator along template string being copied
class Replicator(object):
    
    # method to intialise object with two parameters
    # Type: point - dict
    #       block - dict    
    def __init__(self, pairings, point, block):
        self.magic = dict(pairings)
        self.bases = tuple(self.magic.values())
        self.point = dict(point)
        self.block = dict(block)
        self.template = str()
        self.copy = str()
        self.timer = 0

    # method of printing current replicator information
    def toPrint(self):
        print('TImer: ', self.timer)
        print('Template: ', self.template)
        print('Copy: ', self.copy)

    # method that returns true if the replicator is currently operating on a template
    def working(self):
        return self.timer > 0

    # method that returns true if the replicator has no template to work with
    def stopped(self):
        return self.template == str()

    # method to randomly choose a template out of a list for the replicator to begin working on
    # Type: templates - list
    def move(self,templates):
        self.template = str(random.choice(templates))
        self.copy = str()
        self.timer = len(self.template)
        
    # method of returning a copy of the string produced by the replicator
    # returns a reversed copy to mimick behaviour of directional copying
    def extract(self):
        return str(self.copy[::-1])

    # method of performing a single replication process
    # strategy:
    #    - get the corresponding pair to the template string base
    #   - point mutate the base
    #   - deplete the base from the pool
    #   - add the base to the end of the copy
    #   - move the position of the replicator via block mutation method
    # Type: pool - Pool
    #       rules - dict
    def replicate(self,pool):
        self.copy = self.template.translate(self.magic)
        self.pointMutation()
        self.blockMutation()
        pool.deplete(self.copy)
        self.setTimer()        
#        base = rules[self.template[self.position]]
#        base = self.pointMutation(base,list(pool.elements))
#        pool.deplete(base)
#        self.copy += (base)
#        self.position = self.blockMutation()

    # method that simulates a point mutation in a single base 
    # Type: base - str
    #       pool - list
    def pointMutation(self):
        points = [(np.random.randint(0,len(self.template)),x) \
                  for x in self.point \
                  for f in range(np.random.binomial(len(self.template), \
                                                    self.point[x]))]
        for c in points:
            if c[1] == 'substitution':
                self.copy = self.copy[:c[0]] + np.random.choice(self.bases) + self.copy[c[0]+1:]
            elif c[1] == 'addition':
                self.copy = self.copy[:c[0]] + np.random.choice(self.bases) + self.copy[c[0]:]
            elif c[1] == 'deletion':
                self.copy = self.copy[:c[0]] + self.copy[c[0]+1:]
                
#        roll = random.random()
#        if roll > max(self.point.values()):
#            return base
#        else:
#            option = random.choice([m for m in self.point if self.point[m] > roll])
#            if option == 'substitution':
#                return random.choice(list(pool))
#            elif option == 'addition':
#                return base + random.choice(list(pool))
#            elif option == 'deletion':
#                return ''
#            else:
#                return base

    # method that simulates a block mutation by moving the position
    def blockMutation(self):
        points = [(np.random.randint(0,len(self.template), size = 2),x) \
                  for x in self.block \
                  for f in range(np.random.binomial(len(self.template), \
                                                    self.block[x]))]
        for c in points:
            if c[1] == 'addition':
                self.copy = self.copy[:c[0][0]] + self.copy[c[0][0]:c[0][1]]*2 + self.copy[c[0][1]:]
            elif c[1] == 'deletion':
                self.copy = self.copy[:c[0][0]] + self.copy[c[0][1]+1:]
#        roll = random.random()
#        if roll > max(self.block.values()):
#            return (self.position + 1)
#        else:
#            option = random.choice([m for m in self.block if self.block[m] > roll])
#            if option == 'deletion':
#                return (self.position + 
#                random.randint(1, max(len(self.template) - self.position-1,1)))
#            elif option == 'addition':
#                return random.randint(0,self.position)
#            else:
#                return (self.position + 1)
    
                
    def setTimer(self):
        self.timer = len(self.copy)
        
    # method that causes the replicator to reset 
    def release(self):
        self.template = str()
        self.copy = str()
        self.timer = 0

# end of class Replicator
        
###################################################


# runs a SpSim and returns the simulation object, exports history data to file
def go(inFile = None, printing = -1, outFile = str(), cull = None, progress = False):
    if outFile is str():
        outFile = 'out'+ time.strftime('%d%m%y_%H%M%S') + '.SIMHIST'
    print('Printing to', outFile)
    if inFile is None:
        print('Running with new simulation')
    else:
        print('Importing data from:', str(inFile))
    try:
        s = SpSim(inFile)
        s.run(printing, cull=cull, progress=progress)
        print('Simulation Finished ...' , end='')
    except Exception as ex:
        print(ex)
    finally:
        try:
            s.export_to(outFile)
            print(' History Exported')
        except:
            print('Export Failed')
        finally:
            return s
