from multilateration import *

#determine resolving sets for H(3, 4) and H(12, 4), DNA 3mers and 12mers, using the ICH algorithm for 3mers and the constructive algorithm for 12mers
def findResSets():
  kmers = sorted([''.join(x) for x in product('ACGT', repeat=3)])
  D = [[hammingDist(a, b) for a in kmers] for b in kmers]
  R = ich(D)
  R = [kmers[r] for r in R]
  print('3mer resolving set via ICH', len(R), R)
  #['AAA', 'ACC', 'CAG', 'GGC', 'CGT', 'AAC']

  T = [r for r in R]
  for i in xrange(20):
    R = [t for t in T]
    while len(R[0]) < 7:
      R = hammingConstruction(R, 'ACGT', rand=True)
      if not checkResolvingHamming(R, len(R[0]), 'ACGT'):
        print('Construction failed', len(R[0]), R)
        exit(0)
      else: print('Success', i, len(R), len(R[0]), R)

  while len(R[0]) < 12:
    R = hammingConstruction(R, 'ACGT', rand=True)
    print('construction', len(R[0]), len(R), R)
    if not checkResolvingHamming(R, len(R[0]), 'ACGT'):
      print('Construction failed building 12-mer', len(R[0]), R)
      exit(0)
    else: print('Success', len(R[0]), len(R), R)
  print('12mer resolving set via construction', len(R), R)
  #['AAAAAAAAAAAA', 'ACCGGGGGGGGG', 'CAGAAAAAAAAA', 'GGCAAAAAAAAA', 'CGTAAAAAAAAA', 'AACAAAAAAAAA', 'AAACAAAAAAAA', 'ACCTAAAAAAAA', 'AAAACAAAAAAA', 'ACCGTAAAAAAA', 'AAAAACAAAAAA', 'ACCGGTAAAAAA', 'AAAAAACAAAAA', 'ACCGGGTAAAAA', 'AAAAAAACAAAA', 'ACCGGGGTAAAA', 'AAAAAAAACAAA', 'ACCGGGGGTAAA', 'AAAAAAAAACAA', 'ACCGGGGGGTAA', 'AAAAAAAAAACA', 'ACCGGGGGGGTA', 'AAAAAAAAAAAC', 'ACCGGGGGGGGT']

  orig3mer = ['TGG', 'CGG', 'GTG', 'GCG', 'GGT', 'GGC']
  orig12mer = ['GGGGGGGGGTGG', 'GGGGGGGGGCGG', 'GGGGGGGGGGTG', 'GGGGGGGGGGCG', 'GGGGGGGGGGGT', 'GGGGGGGGGGGC', 'GGGGGGGGTGGC', 'GGGGGGGGCGGC', 'GGGGGGGTCGGC', 'GGGGGGGCCGGC', 'GGGGGGTCCGGC', 'GGGGGGCCCGGC', 'GGGGGTCCCGGC', 'GGGGGCCCCGGC', 'GGGGTCCCCGGC', 'GGGGCCCCCGGC', 'GGGTCCCCCGGC', 'GGGCCCCCCGGC', 'GGTCCCCCCGGC', 'GGCCCCCCCGGC', 'GTCCCCCCCGGC', 'GCCCCCCCCGGC', 'TCCCCCCCCGGC', 'CCCCCCCCCGGC']
  #print('original 3mer res set used, TGG, CGG, GTG, GCG, GGT, GGC', checkResolvingHamming(orig3mer, 3, 'ACGT')) #True
  #print('original 12mer res set used '+','.join(orig12mer), checkResolvingHamming(orig12mer, 12, 'ACGT')) #True


#############
### TESTS ###
#############
if __name__=='__main__':
  if True: #find res set, testing hamming construction etc
    print('find res sets')
    findResSets()

  if True: #test finding all resolving sets for a given hamming graph
    print('TEST ALL HAMMING')
    alphabet = '01'
    for k in xrange(1, 4):
      for size in [s for s in range(1, 4) if s<=k]:
        L = hammingAllResolving(k, alphabet, size, verbose=False)
        print(k, size, len(L))
        #for l in L: print('   ', l)
    alphabet = ['first', 'second']
    for k in xrange(1, 4):
      for size in [s for s in range(1, 4) if s<=k]:
        L = hammingAllResolving(k, alphabet, size, verbose=False)
        print(k, size, len(L))
        #for l in L: print('   ', l)

  if True: #test hamming construction
    print('TEST HAMMING CONSTRUCTION')
    alphabet = '01'
    R = []
    for k in xrange(1, 6):
      R = hammingConstruction(R, alphabet, rand=True)
      if not checkResolvingHamming(R, k, alphabet):
        print('construction failed', R, k, alphabet)
        exit(0)
      else: print('Success', R, k, alphabet)

    alphabet = '012'
    R = []
    for k in xrange(1, 5):
      R = hammingConstruction(R, alphabet, rand=True)
      if not checkResolvingHamming(R, k, alphabet):
        print('construction failed', R, k, alphabet)
        exit(0)
      else: print('Success', R, k, alphabet)

    alphabet = ['first', 'second', 'third']
    R = []
    for k in xrange(1, 5):
      R = hammingConstruction(R, alphabet, rand=True)
      if not checkResolvingHamming(R, k, alphabet):
        print('construction failed', R, k, alphabet)
        exit(0)
      else: print('Success', R, k, alphabet)

  n = 200
  repeats = 50
  if True: #test ich
    print('TEST ICH')
    #graphs, matrix, dict
    for i in xrange(repeats):
      print('ER random graphs', i)
      p = np.random.choice(np.arange(0.1, 1., 0.1))
      G = nx.gnp_random_graph(n, p)
      nodes = sorted(G.nodes())
      print('   make structures')
      M = nx.floyd_warshall(G)
      M = [[int(M[u][v]) if (v in M[u] and M[u][v]!=np.inf) else -1 for v in nodes] for u in nodes]
      D = {u:{v:M[u][v] for v in xrange(len(M[u])) if M[u][v]!=-1} for u in nodes}
      useDist = np.random.choice([True, False], p=[0.1, 0.9])
      print('   ich graph')
      GR = ich(G, randOrder=False, useFullDistance=useDist)
      if not checkResolving(GR, G):
        print('graph not resolving', useDist)
        exit(0)
      print('   ich matrix')
      MR = ich(M, randOrder=False)
      if not checkResolving(MR, M):
        print('matrix not resolving')
        exit(0)
      print('   ich dict')
      DR = ich(D, colNames=D.keys(), randOrder=False)
      if not checkResolving(DR, D):
        print('dict not resolving')
        exit(0)
      if GR!=MR or MR!=DR:
        print('Fail', useDist, GR, MR, DR)
        exit(0)

  if True:
    for i in xrange(repeats):
      print('Pref attach', i)
      m = np.random.choice(range(1, 10, 1))
      G = nx.barabasi_albert_graph(n, m)
      nodes = sorted(G.nodes())
      print('   make structures')
      M = nx.floyd_warshall(G)
      M = [[int(M[u][v]) if (v in M[u] and M[u][v]!=np.inf) else -1 for v in nodes] for u in nodes]
      D = {u:{v:M[u][v] for v in xrange(len(M[u])) if M[u][v]!=-1} for u in nodes} 
      useDist = np.random.choice([True, False], p=[0.1, 0.9])
      print('   ich graph')
      GR = ich(G, randOrder=False, useFullDistance=useDist)
      if not checkResolving(GR, G):
        print('graph not resolving', useDist)
        exit(0)
      print('   ich matrix')
      MR = ich(M, randOrder=False)
      if not checkResolving(MR, M):
        print('matrix not resolving')
        exit(0)
      print('   ich dict')
      DR = ich(D, colNames=D.keys(), randOrder=False)
      if not checkResolving(DR, D):
        print('dict not resolving')
        exit(0)
      if GR!=MR or MR!=DR:
        print('Fail', useDist, GR, MR, DR)
        exit(0)

  if True:
    for i in xrange(repeats):
      print('Watts Strogatz', i)
      k = np.random.choice(range(2, 5, 1))
      p = np.random.choice(np.arange(0.1, 1., 0.1))
      G = nx.watts_strogatz_graph(n, k, p)
      nodes = sorted(G.nodes())
      print('   make structures', k, p)
      M = nx.floyd_warshall(G)
      M = [[int(M[u][v]) if (v in M[u] and M[u][v]!=np.inf) else -1 for v in nodes] for u in nodes]
      D = {u:{v:M[u][v] for v in xrange(len(M[u])) if M[u][v]!=-1} for u in nodes}
      useDist = np.random.choice([True, False], p=[0.1, 0.9])
      print('   ich graph')
      GR = ich(G, randOrder=False, useFullDistance=useDist)
      if not checkResolving(GR, G):
        print('graph not resolving', useDist)
        exit(0)
      print('   ich matrix')
      MR = ich(M, randOrder=False)
      if not checkResolving(MR, M):
        print('matrix not resolving')
        exit(0)
      print('   ich dict')
      DR = ich(D, colNames=D.keys(), randOrder=False)
      if not checkResolving(DR, D):
        print('dict not resolving')
        exit(0)
      if GR!=MR or MR!=DR:
        print('Fail', useDist, GR, MR, DR)
        exit(0)
    
  if True: #matrix, dict
    for i in xrange(repeats):
      print('matrix and dict', i)
      x = np.random.choice([4,10,20])
      M = [[np.random.choice(range(x)) for _ in xrange(n)] for _ in xrange(n)]
      D = {u:{v:M[u][v] for v in xrange(len(M[u])) if M[u][v]!=-1} for u in xrange(len(M))}
      print('   ich matrix')
      MR = ich(M, randOrder=False)
      if not checkResolving(MR, M):
        print('matrix not resolving')
        exit(0)
      print('   ich dict')
      DR = ich(D, colNames=D.keys(), randOrder=False)
      if not checkResolving(DR, D):
        print('dict not resolving')
        exit(0)
      if MR!=DR:
        print('Fail', MR, DR)
        exit(0)

  if True: #test procs
    procs = 5
    for i in xrange(repeats):
      print('test procs', procs)
      p = np.random.choice(np.arange(0.1, 1., 0.1))
      G = nx.gnp_random_graph(n, p)
      nodes = sorted(G.nodes())
      print('   make structures')
      M = nx.floyd_warshall(G)
      M = [[M[u][v] if (v in M[u] and M[u][v]!=np.inf) else -1 for v in nodes] for u in nodes]
      D = {u:{v:M[u][v] for v in xrange(len(M[u])) if M[u][v]!=-1} for u in nodes}
      useDist = np.random.choice([True, False], p=[0.1, 0.9])
      print('   ich graph')
      GR = ich(G, randOrder=False, useFullDistance=useDist, procs=procs)
      if not checkResolving(GR, G):
        print('graph not resolving', useDist)
        exit(0)
      print('   ich matrix')
      MR = ich(M, randOrder=False, procs=procs)
      if not checkResolving(MR, M):
        print('matrix not resolving')
        exit(0)
      print('   ich dict')
      DR = ich(D, colNames=D.keys(), randOrder=False, procs=procs)
      if not checkResolving(DR, D):
        print('dict not resolving')
        exit(0)
      print('testing proc on er', i, len(GR), len(MR), len(DR))
    
  if True: #test save... fully resolving (make sure none added)
    for i in xrange(5):
      print('ER random graphs, save full', i)
      p = np.random.choice(np.arange(0.1, 1., 0.1))
      G = nx.gnp_random_graph(n, p)
      nodes = sorted(G.nodes())
      print('   make structures')
      M = nx.floyd_warshall(G)
      M = [[int(M[u][v]) if (v in M[u] and M[u][v]!=np.inf) else -1 for v in nodes] for u in nodes]
      D = {u:{v:M[u][v] for v in xrange(len(M[u])) if M[u][v]!=-1} for u in nodes}
      useDist = np.random.choice([True, False], p=[0.1, 0.9])
      print('   ich graphs')
      GR = ich(G, randOrder=False, useFullDistance=useDist, name='gr')
      GR2 = ich(G, randOrder=False, useFullDistance=useDist, stateFile='gr_ich_progress')
      if GR!=GR2 or not checkResolving(GR, G):
        print('graph not resolving', useDist)
        exit(0)
      print('   ich matrices')
      MR = ich(M, randOrder=False, name='mr')
      MR2 = ich(M, randOrder=False, stateFile='mr_ich_progress')
      if MR!=MR2 or not checkResolving(MR, M):
        print('matrix not resolving')
        exit(0)
      print('   ich dicts')
      DR = ich(D, colNames=D.keys(), randOrder=False, name='dr')
      DR2 = ich(D, colNames=D.keys(), randOrder=False, stateFile='dr_ich_progress')
      if DR!=DR2 or not checkResolving(DR, D):
        print('dict not resolving')
        exit(0)
      if GR!=MR or MR!=DR:
        print('Fail', useDist, GR, MR, DR)
        print('....', useDist, GR2, MR2, DR2)
        exit(0)
    
  if True: #test save... partially resolving random small subset of nodes
    for i in xrange(5):
      print('Pref attach, save partial', i)
      m = np.random.choice(range(1, 10, 1))
      G = nx.barabasi_albert_graph(n, m)
      nodes = sorted(G.nodes())
      print('   make structures')
      M = nx.floyd_warshall(G)
      M = [[int(M[u][v]) if (v in M[u] and M[u][v]!=np.inf) else -1 for v in nodes] for u in nodes]
      D = {u:{v:M[u][v] for v in xrange(len(M[u])) if M[u][v]!=-1} for u in nodes} 
      chosen = list(np.random.choice(xrange(len(M)), size=4))
      tags = {u:';'.join(map(str, [M[u][r] for r in chosen])) for u in xrange(len(M))}
      saveState(tags, chosen, 'temp_file', overwrite=True) 
      useDist = np.random.choice([True, False], p=[0.1, 0.9])
      print('   ich graph')
      GR = ich(G, randOrder=False, useFullDistance=useDist, stateFile='temp_file')
      if not checkResolving(GR, G):
        print('graph not resolving', useDist)
        exit(0)
      print('   ich matrix')
      MR = ich(M, randOrder=False, stateFile='temp_file')
      if not checkResolving(MR, M):
        print('matrix not resolving')
        exit(0)
      print('   ich dict')
      DR = ich(D, colNames=D.keys(), randOrder=False, stateFile='temp_file')
      if not checkResolving(DR, D):
        print('dict not resolving')
        exit(0)
      if GR!=MR or MR!=DR:
        print('Fail', useDist, GR, MR, DR)
        exit(0)

  print('DONE')



