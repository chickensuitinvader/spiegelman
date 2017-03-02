# -*- coding: utf-8 -*-
"""
A program to run multiple Spigelman's Monster simulations in order for 
parametric analysis
"""

import os, time, ast, importlib, pylab
import numpy as np
import spiegelman
import parameters as defaults
from matplotlib.collections import LineCollection
from matplotlib.colors import colorConverter as colourConverter


def multipleRuns(number, inputs = None):
    sim_folder = 'Sim_' + time.strftime('%d%m%y_%H%M%S')
    os.mkdir(sim_folder)
    sims_list = list()
    for N in range(number):
        sims_list.append(spiegelman.go(inputs, 0, 
                                       sim_folder + '/' + str(N) + '.SIMHIST'))
    return sims_list
    
# temp class
class pSet(object):
    
    def __init__(self,p):
        if isinstance(p, pSet):
            self.parameters = dict(p.parameters)
        else:
            self.parameters = dict(p)
        
    def set(self,key,value):
        self.parameters[key] = type(value)(value)
        
    def get(self, key):
        return self.parameters[key]
        
    def list(self):
        return self.parameters
    
# interactive parameter runs
def parameterRuns():
    importlib.reload(defaults)
    sims_list = list()
    
    print('Default Parameters:')
    [print(s,':',defaults.parameters[s]) for s in defaults.parameters]
    
    # extract information
    totalRuns = int(input('How many runs are required?  '))
    print('Parameters Available:')
    print([p for p in defaults.parameters])
    param = input('Which Parameter would you like to change?  ')
    if not defaults.parameters.keys().__contains__(param):
        print('Invalid parameter')
        raise(KeyError)
    elif isinstance(defaults.parameters[param], dict):
        sub_params = [sub for sub in defaults.parameters[param]]
        print('Sub-parameters to tune: ', sub_params)
        prange = dict()
        parameters = dict()
        for n in sub_params:
            prange[n] = input('Please enter the range for the sub-parameter '+ n+
                  ' (seperate with ::)  ').split('::')
            if len(prange[n]) != 2:
                print('Non standard ranges not yet implemented.')
                raise(Exception)
        parameters = [{n : (float(prange[n][0]) + k*(float(prange[n][1]) - \
                                  float(prange[n][0]))/(totalRuns-1)) for \
                                  n in sub_params} for k in range(totalRuns) ]
    elif isinstance(defaults.parameters[param], (float,int)):      
        prange = input('Please enter the range of the parameter (separate with ::)  ')
        inputs = prange.count('::') + 1
        prange = prange.split('::')
        if inputs == totalRuns:
            parameters = list()
            for k in prange:
                parameters.append(ast.literal_eval(k))
        elif inputs == 2:
            pMin = float(prange[0])
            pMax = float(prange[1])
            parameters = [pMin + k*(pMax-pMin)/(totalRuns-1) for k in range(totalRuns)]
        else:
            print('Format has not been implemented. Try MIN::MAX or give one arg/run')
            raise(Exception)
    elif isinstance(defaults.parameters[param], (list,tuple)):
        print('Parameter is a tuple or list, not implemented yet')
        raise(Exception)
    else:
        print('Unknown parameter Type')
        raise(Exception)
            
        
    # make parameter inputs
    psSet = list()
    for n in range(totalRuns):
        df = defaults.parameters
        df[param] = parameters[n]
        psSet.append(pSet(dict(df)))
    print('Simulations Made, now running')

    # run sims
    try:
        sim_folder = 'Sim_' + time.strftime('%d%m%y_%H%M%S')
        os.mkdir(sim_folder)
        for n in range(totalRuns):
            name = sim_folder + '/' + str(n) + '.SIMHIST'
            sims_list.append(spiegelman.go(psSet[n], -1, name))
    except Exception as ex:
        print(ex)
    finally:
        print('Returning Parameter List and Simulation List')
        return [psSet, sims_list]
"""        
#def paramComparisons(paramOutput, comp_pnt = None, position = -1):
#    pList = paramOutput[0]
#    sims_list = paramOutput[1]
#    if comp_pnt == None:
#        print('Categories', list(sims_list[0].history.keys()))
#        comp_pnt = input('What category to compare?  ')
#    if comp_pnt not in sims_list[0].history.keys():
#        print('Key Error, Quitting')
#        return
#    
#    ind_pnt = [k for k in pList[1].parameters
#               if pList[1].parameters[k] != pList[0].parameters[k]]
#    if len(ind_pnt) != 1:
#        print('Parameters that Changed:', ind_pnt)
#        print('All Parameters:', [k for k in pList[1].parameters])
#        ind_pnt = input('Which parameter to index with?  ')
#    else:
#        ind_pnt = ind_pnt[0]
#        
#    print('Comparing:', comp_pnt)
#    indices = list()
#    values = list()
#    for n in sims_list:
#        indices.append(pList[sims_list.index(n)].parameters[ind_pnt])
#        comp_vals = n.history[comp_pnt][position]
#        if isinstance(comp_vals, dict):
#            values.append([comp_vals[k] for k in comp_vals])
#        else:
#            values.append(comp_vals)
#        #print(indices[-1],values[-1])
#    
#    spiegelman.plt.plot(indices,values,'.')
#    spiegelman.plt.xlabel(ind_pnt)
#    spiegelman.plt.ylabel(comp_pnt)
#    spiegelman.plt.title(comp_pnt + ' by ' + ind_pnt + ' in Iteration' + str(position))
#    
#    
#    return [indices,values]
"""    
def plotComparison(paramOutput, heat = False, comp_pnt = None):
    pList = paramOutput[0]
    sList = paramOutput[1]
    multigraph = False
    
    #determine output to compare
    if comp_pnt == None:
        print('Categories', list(sList[0].history.keys()))
        comp_pnt = input('What category to compare?  ')
    if comp_pnt not in sList[0].history.keys():
        print('Key Error, Quitting')
        return
    
    #determine parameter to compare by, and its values
    ind_pnt = [k for k in pList[1].parameters
               if pList[1].parameters[k] != pList[0].parameters[k]]
    if (len(ind_pnt) != 1):
        print('Parameters that Changed:', ind_pnt)
        print('All Parameters:', [k for k in pList[1].parameters])
        ind_pnt = input('Which parameter to index with?  ')
        indices = [p.parameters[ind_pnt] for p in pList]
    elif (len(ind_pnt[0]) != 1):
        print('Parameters that Changed:', ind_pnt[0])
        print('Sub-Parameters:', [k for k in pList[0].parameters[ind_pnt[0]] 
                                  if (pList[0].parameters[ind_pnt[0]][k] != 
                                   pList[1].parameters[ind_pnt[0]][k])])
        sub_ind = input('Which sub-parameter to index with?  ')
        indices = [p.parameters[ind_pnt[0]][sub_ind] for p in pList]
        ind_pnt = ind_pnt[0] + '.' + sub_ind
    else:
        ind_pnt = ind_pnt[0]
        print('Parameter: ', ind_pnt)
        indices = [p.parameters[ind_pnt] for p in pList]
    
    #determine values
    epochs = [range(s.parameters['Cycles']) for s in sList]
    if isinstance(sList[0].history[comp_pnt][0],dict):
        
        values = [[[k[j] for j in k] for k in s.history[comp_pnt]] for s in sList]
        multigraph = True
    else:
        values = [s.history[comp_pnt] for s in sList]
    
    if heat == True:
        heatMap(comp_pnt, ind_pnt, indices, epochs, values)
        return
        
    #determine plot type
    if comp_pnt == 'Lengths':
        #use hist
        print('Histogram not implemented yet')
        return
    elif multigraph:
        print('Multigraph not yet implemented')
        return
    else:
        vert = [list(zip(epochs[values.index(v)],v)) for v in values]
    
    fig = spiegelman.plt.figure()
    ax = spiegelman.Axes3D(fig)           

    colours = [colourConverter.to_rgba(i) for i in 'bgrcmyk']                
    lines = LineCollection(vert, colors=colours)
    lines.set_alpha(0.7)
    lines.set_linewidth(1.5)
    ax.add_collection3d(lines, zs=indices, zdir='y')
    
    ax.set_xlabel('Epoch')
    ax.set_ylabel(ind_pnt)
    ax.set_zlabel(comp_pnt)
    
    ax.set_xlim3d(min([i for j in epochs for i in j]), max([i for j in epochs for i in j]))
    ax.set_ylim3d(min(indices),max(indices))
    ax.set_zlim3d(min([i for j in values for i in j]), max([i for j in values for i in j]))
    spiegelman.plt.show()
    return [indices,vert]
    
def getOutputs(folder):
    start = time.time()
    outputs = list([list(),list()])
    for fileName in os.listdir(folder):
        if fileName[-8:] == '.SIMHIST':
            outputs[1].append(spiegelman.SpSim(folder+'/'+fileName))
            outputs[0].append(pSet(outputs[1][-1].parameters))
    print('Outputs Obtained in', round(time.time() - start, 3), 's')
    return outputs
            
def customRun(specs = None):
    if isinstance(specs,(pSet, str, spiegelman.SpSim)):
        s = spiegelman.go(specs)
    else:
        print('Custom Run, manually changing parameters')
        ps = dict()
        for n in defaults.parameters:
            print('Changing', n,'\nDefault:', defaults.parameters[n])
            got = input('Please enter the value for this parameter:  ')
            ps[n] = ast.literal_eval(got) if got != str() else defaults.parameters[n]
        s = spiegelman.go(pSet(ps))
    return s

def heatMap(comp_pnt, ind_pnt, indices, epochs, values):
    spiegelman.plt.close()
    indices = np.tile(np.array(indices),(len(epochs[0]),1)).transpose()
    epochs = np.array(epochs)
    values = np.array(values)
    spiegelman.plt.contourf(epochs,indices,values)
    spiegelman.plt.colorbar()
    return

def heatDists(outs,folder):
    try:
        os.mkdir(folder+'/Dists')
    except Exception as ex:
        print(ex)
    for i in outs[1]:
        i.plotting3('h', heat=True)
        fm = spiegelman.plt.get_current_fig_manager()
        fm.window.showMaximized()
        pylab.savefig(folder + '/Dists/'+ str(outs[1].index(i)) + '.png')
        
def cull_runs(ns = None, lims = None, plot = False):
    if ns == None:
        ns = int(input('Enter number of runs: '))
    if lims == None:
        llim = ast.literal_eval(input('Enter lower bound: '))
        ulim = ast.literal_eval(input('Enter upper bound: '))
    else:
        llim = min(lims)
        ulim = max(lims)
    ps = np.transpose([np.linspace(llim[i],ulim[i], num=ns) for i in range(len(llim))])
    
    ss = list()
    outFolder = time.strftime('%d%m%y_%H%M')
    os.mkdir(outFolder)
    for p in ps:
        def cullfn(x):
            return spiegelman.cull_function.cull_function(x,None,p)
            
        print('p = ', p)    
        out = outFolder + '/' + time.strftime('%d%m%y_%H%M%S.SIMHIST')
        ss.append(spiegelman.SpSim())
        ss[-1].parameters['CullParameter'] = tuple(p)
        ss[-1].run(cull = cullfn)
        ss[-1].export_to(out)

    if plot:
        cplotting(ss, outFolder)

    return ss

def cplotting(ss, outFolder=None):
    if outFolder == None:
        try: 
            os.mkdir('Temp')
        finally:
            outFolder = 'Temp/'
    os.mkdir(outFolder+'/Pics')
    cull_plot(ss)
    spiegelman.plt.savefig(outFolder+'/Pics/cplot.png')
    for s in ss:
        N = str(ss.index(s))
        for g in ('a','p','u%','n'):
            s.plotting(g)
            figM = spiegelman.plt.get_current_fig_manager()
            figM.window.showMaximized()
            spiegelman.plt.savefig(outFolder+'/Pics/'+N+g+'.png')
        s.plotting3('h')
        figM = spiegelman.plt.get_current_fig_manager()
        figM.window.showMaximized()
        spiegelman.plt.savefig(outFolder+'/Pics/'+N+'h.png')

def cull_plot(ss, cm = spiegelman.CM.RdBu):
    averages = np.transpose([x.history['Average'] for x in ss])
    try:
        params = [np.sqrt(np.sum(q**2))  for x in ss for q in x.parameters['CullParameter']]
        plabel = 'Normed Parameter Value'
    except:
        params = range(np.size(averages,1))
        plabel = 'Parameter Number'
    epochs = range(np.size(averages,0))
    #return(averages,params,epochs)
    spiegelman.plt.close()
    extents = [params[0],params[-1],epochs[0],epochs[-1]]
    f,ax = spiegelman.plt.subplots(1,1,figsize=(6,6))
    im = ax.imshow(np.log10(averages), extent=extents, interpolation='None', cmap = cm, aspect='auto')
    f.colorbar(im)
    ax.set_title('Logged Average Lengths')
    ax.set_xlabel(plabel)
    ax.set_ylabel('Epoch')