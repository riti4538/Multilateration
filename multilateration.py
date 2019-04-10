#Python code implementing the ICH (Information Content Heuristic) approximation algorithm for metric dimension (see [1])
#This code also applies to general multilateration problems (see [2])
#Code implementing a constructive algorithm for a bound on the metric dimension of Hamming graphs is included (see [2]) along with functions for testing resolvability
#[1] Hauptmann, M., Schmied, R., and Viehmann, C. Approximation complexity of metric dimension problem, Journal of Discrete Algorithms 14(2012), 214-222.
#[2] Tillquist, R. C., and Lladser, M. E. Low-dimensional representation of genomic sequences. Journal of Mathematical Biology (Mar 2019).

import numpy as np
import networkx as nx
from scipy.stats import entropy
from scipy.misc import comb
from itertools import product, combinations
import multiprocessing as mp
import pickle

########################
### READ/WRITE STATE ###
########################
#Save a snapshot of the algorithm state as a list of tuples (tags, chosen elements) using Python pickle
#input: tags - a dictionary of tags for each element based on the chosen columns
#       chosen - a list of the columns chosen so far
#       saveFile - name of the file to write to
#       overwrite - if true overwrite the contents of the file, otherwise append, default is False
def saveState(tags, chosen, saveFile, overwrite=False):
  L = [] if overwrite else readState(saveFile)
  L.append((tags, chosen))
  with open(saveFile, 'wb') as o:
    pickle.dump(L, o, protocol=pickle.HIGHEST_PROTOCOL)

#Read a list of saved states of the algorithm 
#input: inFile - the file to read
#return: a list of (tag, chosen) pairs
def readState(inFile):
  L = []
  with open(inFile, 'rb') as f:
    L = pickle.load(f)
  return L

#####################
### ICH Algorithm ###
#####################
#Apply the ICH algorithm to approximate metric dimension or multilateration
#input: M - a list of lists, dictionary of dictionaries, or networkx graph object on which to perform the ICH algorithm
#           if a list of lists, the ICH algorithm will be applied directly
#           if a dictionary of dictionaries, the colNames argument will be used as a complete list of columns. if this is empty, an arbitrary key list from an elemt of M will be used
#           if a networkx graph object, the distance matrix will be used
#       colNames - optional list of columns. if M is a dictionary of dictionaries and colNames is not empty it will be used as the set of columns to consider, defaults to empty
#       dictDefault - optional default value to use if M is a dictionary of dictionaries, defaults to -1
#       useFullDistance - optional boolean value. if true and M is a graph, convert M to a list of lists before continuing. note that distances in the graph are assumed to be positive.
#                         if false, individual columns of a distance matrix are generated on demand in the colEntropy function
#       name - optional prefix to give to a file in which to save current state of the algorithm, defaults to empty
#       stateFile - optional name of a file to read current state from, default is empty
#       randOrder - optional boolean value. if true randomize the order in which columns are checked, defaults to true
#       procs - optional number of processes to run, defaults to 1
#return: a list of chosen columns representing a resolving set
def ich(M, colNames=[], dictDefault=-1, useFullDistance=True, name='', stateFile='', randOrder=True, procs=1):
  progressFile = name+'_ich_progress'
  if name: saveState({}, [], progressFile, overwrite=True)
  if isinstance(M, nx.classes.graph.Graph) and useFullDistance:
    nodes = sorted(M.nodes())
    M = nx.floyd_warshall(M)
    M = [[int(M[u][v]) if (v in M and M[u][v]!=np.inf) else -1 for v in nodes] for u in nodes]
  distr = {}
  tags = {}
  chosen = []
  if stateFile:
    (tags, chosen) = readState(stateFile)[-1]
    for t in tags:
      if tags[t] not in distr: distr[tags[t]] = []
      distr[tags[t]].append(t)
  elif isinstance(M, list): tags = {i:'' for i in xrange(len(M))}
  elif isinstance(M, dict): tags = {i:'' for i in M}
  elif isinstance(M, nx.classes.graph.Graph): tags = {i:'' for i in M.nodes()}
  n = len(tags)
  check = colNames if (isinstance(M, dict) and colNames) else (sorted(M[M.keys()[0]].keys()) if isinstance(M, dict) else sorted(tags.keys()))
  for col in chosen: check.remove(col)
  check = list(np.random.permutation(check)) if randOrder else sorted(check)
  while len(distr) < n and len(chosen) < n:
    (distr, tags, _, chosen) = pickColumn(M, tags, check, chosen=chosen, dictDefault=dictDefault, procs=procs)
    check.remove(chosen[-1])
    if name: saveState(tags, chosen, progressFile, overwrite=True)
  if len(distr) < n and len(chosen) == n: return 'NO SOLUTION EXISTS'
  return chosen

#Determine the unchosen column maximizing change in entropy
#input: M - the object on which to perform multilateration
#       tags - tages given the columns already chosen
#       check - a list of columns left to check
#       chosen - optional list of already chosen columns, defaults to empty
#       dictDefault - optional default value to use if M is a dictionary of dictionaries, defaults to -1
#       procs - optional argument specifying the number of processes to use, defaults to 1
#return: the tag distribution with the new column, the new tags, the entropy with the new column, an updated list of chosen columns
def pickColumn(M, tags, check, chosen=[], dictDefault=-1, procs=1):
  (eMax, cMax, distr) = (-1, -1, {})
  if procs > 1:
    pool = mp.Pool(processes=procs)
    results = pool.map_async(colEntropy, [(col, M, tags, dictDefault) for col in check])
    results = results.get()
    pool.close()
    pool.join()
    (eMax, cMax, distr) = max(results)
  else:
    for col in check:
      (e, col, condDistr) = colEntropy((col, M, tags, dictDefault))
      if e > eMax: (eMax, cMax, distr) = (e, col, condDistr)
  chosen.append(cMax)
  if isinstance(M, list):
    for t in tags: tags[t] += ';'+str(M[t][cMax])
  elif isinstance(M, dict):
    for t in tags: tags[t] += ';'+str(M[t].get(cMax, dictDefault))
  elif isinstance(M, nx.classes.graph.Graph):
    Dcol = nx.single_source_shortest_path_length(M, cMax)
    for t in tags: tags[t] += ';'+str(Dcol.get(t, -1))
  return (distr, tags, eMax, chosen)

#Determine the joint entropy resulting from adding a given column
#Input is given as a single tuple so map_async may be used from pickColumn
#input: col - the column to add
#       M - the object on which to perform multilateration
#       tags - tags given the columns already chosen
#       dictDefault - default value to use if M is a dictionary of dictionaries
#return: the joint entropy, the column, and the distribution of tags
def colEntropy((col, M, tags, dictDefault)):
  getSymbol = lambda r,c: ''
  Dcol = {}
  if isinstance(M, list): getSymbol = lambda r,c: M[r][c]
  elif isinstance(M, dict): getSymbol = lambda r,c: M[r].get(c, dictDefault)
  elif isinstance(M, nx.classes.graph.Graph):
    Dcol = nx.single_source_shortest_path_length(M, col)
    getSymbol = lambda r,c: Dcol.get(r, -1)
  jointDistr = {}
  for elem in tags:
    t = tags[elem]+';'+str(getSymbol(elem, col))
    if t not in jointDistr: jointDistr[t] = 0
    jointDistr[t] += 1
  e = entropy(jointDistr.values(), base=2)
  return (e, col, jointDistr)

######################
### HAMMING GRAPHS ###
######################
#Computes the hamming distance or number of mismatches between two strings
#If one string is longer the other, only its prefix is used
#input: a, b - two sequences to compare
#return: the hamming distance between a and b
def hammingDist(a, b):
  return sum(1 for (x,y) in zip(a,b) if x!=y)

#Given a resolving set of a Hamming graph H(k, a), determine a resolving set for H(k+1, a) (see [2])
#input: resSet - a resolving set for H(k, a)
#       alphabet - the alphabet from which to draw characters for the new resolving set
#       rand - optional boolean, if true randomize resSet and alphabet order, default is false
#return: a resolving set for H(k+1, a)
def hammingConstruction(resSet, alphabet, rand=False):
  alphabet = [[a] for a in alphabet]
  if len(resSet)==0: return alphabet[:-1]
  if rand:
    resSet = map(list, np.random.permutation(resSet))
    alphabet = map(list, np.random.permutation(alphabet))
  newResSet = [r+alphabet[2*i] if 2*i<len(alphabet) else r+alphabet[0] for i,r in enumerate(resSet)]
  num = len(alphabet) / 2
  for i in xrange(num):
    v = resSet[i]+alphabet[2*i+1]
    newResSet.append(v)
  return newResSet

#Find all resolving sets of a Hamming graph via a brute force search for a particular size
#This may be extremely slow even for small values of k and alphabet
#input: k - the length of strings in the hamming graph
#       alphabet - the alphabet to use in the hamming graph
#       size - the size of sets to check
#       verbose - optional bool. if true, print percent of sets checked every 10000 sets, default is false
#return: all resolving sets of the given size for the specified hamming graph
def hammingAllResolving(k, alphabet, size, verbose=False):
  if isinstance(alphabet, str): alphabet = list(alphabet)
  resSets = []
  kmers = product(alphabet, repeat=k)
  numCombos = comb(int(np.power(len(alphabet), k)), size)
  for i,R in enumerate(combinations(kmers, size)):
    if verbose and i%10000==0: print('Brute force progress: ', i / numCombos)
    if checkResolvingHamming(R, k, alphabet, verbose=verbose): resSets.append(R)
  return resSets

#Check that a given set of strings is resolving for a specified Hamming graph
#This may be extremely slow even for small values of k and alphabet
#input: R - a set of strings to check as resolving
#       k - length of strings
#       alphabet - characters that strings are composed of
#       verbose - optional bool. if true, print percent of strings checked every 10000 sets, default is false
#return: true if R is resolving on H(k, |alphabet|) and false otherwise
def checkResolvingHamming(R, k, alphabet, verbose=False):
  if isinstance(alphabet, str): alphabet = list(alphabet)
  tags = {}
  tot = float(np.power(len(alphabet), k))
  for i,seq in enumerate(product(alphabet, repeat=k)):
    if verbose and i%10000==0: print('Check resolving Hamming progress: ', i/tot)
    tag = ';'.join(map(str, [hammingDist(list(seq), r) for r in R]))
    if tag in tags: return False
    tags[tag] = 1
  return True

###########################
### CHECK RESOLVABILITY ###
###########################
#Given a set of columns and an object on which to check resolvability, check that the set is resolving
#input: R - a set of columns
#       M - an object on which to check the resolvability of R
#       colNames - optional list of columns. if M is a dictionary of dictionaries and colNames is not empty it will be used as the set of columns to consider, defaults to empty
#       dictDefault - optional default value to use if M is a dictionary of dictionaries, defaults to -1
#return: true if R is resolving and false otherwise
def checkResolving(R, M, colNames=[], dictDefault=-1):
  tags = {}
  elements = []
  if isinstance(M, list): elements = xrange(len(M))
  elif isinstance(M, dict): elements = M.keys()
  elif isinstance(M, nx.classes.graph.Graph): elements = M.nodes()
  for elem in elements:
    tag = []
    if isinstance(M, list): tag = [M[elem][r] for r in R]
    elif isinstance(M, dict): tag = [M[elem].get(r, dictDefault) for r in R]
    elif isinstance(M, nx.classes.graph.Graph): #tag = [nx.shortest_path_length(M, source=elem, target=r) for r in R]
      for r in R:
        try:
          tag.append(nx.shortest_path_length(M, source=elem, target=r))
        except:
          tag.append(-1)
    tag = ';'.join(map(str, tag))
    if tag in tags: return False
    tags[tag] = 1
  return True






