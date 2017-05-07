"""Detection of Interspersed Signatures in SEQuence Data

general usage:
disseqd -d datafile -m modelfile [options]

basic options:
-h / --help             prints this message
-d / --data datafile    file with input sequence
-m / --model modelfile  file with model parameters

output options:
-v / --verbose          more detailed output
-o / --ouput prefix     outputs to files prefix.model.txt and prefix.decode.txt

fitting options:
-f / --fit              fit model rather than decode data
-r / --rounds nrounds   fits model for specified number of rounds
-l / --lock s/t/e       locks any combination of start/transitions/emissions

override options: (overrides options in modelfile)
-k / --kmer kmer                                order of model plus 1
-n / --nstates (nstates / state1,state2[,...])  number or names of states
-a / --alphabet alphabet                        alphabet used for model
-s / --start start1,start2[,...]                start probabilities
-t / --transitions length1,length2[,...]        expected lengths of states
-e / --emissions datafile1,datafile2[,...]      files with state observations

"""

import getopt, sys
import numpy
import itertools
import re
from os.path import expanduser

config = None

def main():
    """Reads model and disseqds data."""
    global config

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hvo:fr:l:d:m:k:n:a:s:t:e:', ['help', 'verbose', 'output=', 'fit', 'rounds=', 'lock=', 'data=', 'model=', 'kmer=', 'nstates=', 'alphabet=', 'start=', 'transitions=', 'emissions='])
    except getopt.GetoptError as err:
        print('unable to disseqd: {0}'.format(err)) # will print something like 'option -a not recognized'
        print(__doc__)
        sys.exit(2)

    for o, a in opts:
        if o in ('-h', '--help'):
            print(__doc__)
            sys.exit()
        elif o in ('-v', '--verbose'):
            verbose = True
        elif o in ('-o', '--output'):
            output = a
        elif o in ('-f', '--fit'):
            fit = True
        elif o in ('-r', '--rounds'):
            rounds = a
        elif o in ('-l', '--lock'):
            lock = a
        elif o in ('-d', '--data'):
            datafile = a
        elif o in ('-m', '--model'):
            modelfile = a
        elif o in ('-k', '--kmer'):
            kmer = a
        elif o in ('-n', '--nstates'):
            nstates = a
        elif o in ('-a', '--alphabet'):
            alphabet = a
        elif o in ('-s', '--start'):
            start = a
        elif o in ('-t', '--transitions'):
            transitions = a
        elif o in ('-e', '--emissions'):
            emissions = a
        else: # does this ever happen?
            assert False, 'unhandled option'

    config = dict(locals())
    validate_input()
    validate_data()

    if(config.get('fit') == True):
        for r in range(config.get('rounds')):
            print('Fitting round {0}'.format(r))
            forward()
            backward()
            gamma()
            xi()
            update()
            decode()
            print(v,ll)
    else:
        decode()
    print(v,ll)
    output_model()

def open_file(f,mode='r'):
    """Returns a filehandle to the file, opened in specified mode."""
    try:
        fh = open(expanduser(f),mode)
    except OSError as err:
        print('{0} {1}'.format(f,err.strerror))
        sys.exit(2)
    return fh

def output_model():
    """Writes the model to file."""
    global config
    global pi,A,B,literal_emission
    print('Writing model to file.')
    config['start'] = pi
    config['transitions'] = A
    config['emissions'] = B

    f = str(config.get('output'))+'.model.txt'
    fh = open_file(f,'w')
    fh.write('>kmer\n{0}\n'.format(config.get('kmer')))
    fh.write('>nstates\n{0}\n'.format('\n'.join(str(literal) for numeric,literal in sorted(literal_state.items()))))
    fh.write('>alphabet\n{0}\n'.format(config.get('alphabet')))
    fh.write('>start\n{0}\n'.format('\n'.join([ str(_) for _ in config.get('start') ])))
    fh.write('>transitions\n')
    [ fh.write('{0}\n'.format('\t'.join(str(prob) for prob in row))) for row in config.get('transitions') ]
    fh.write('>emissions\n')
    [ fh.write('{0}\t{1}\n'.format(label,'\t'.join(str(prob) for prob in row))) for label,row in zip(literal_emission,config.get('emissions').transpose()) ]
    fh.write('\n')

    fh.close()

    f = str(config.get('output'))+'.decode.txt'
    fh = open_file(f,'w')
    fh.write('>loglikelihood:'+str(ll))
    [ fh.write('\t{0}:{1}'.format(literal,numeric)) for numeric,literal in sorted(literal_state.items()) ]
    fh.write('\n'+v+'\n')
    fh.close()

def validate_input():
    """Validates the input."""
    global config

    if(config.get('output') == None):
        print('No output specified, using "{0}".'.format(config.setdefault('output','out')))

    if(config.get('datafile')):
        print('Using data in {0}'.format(config.get('datafile')))
        if(config.get('modelfile')):
            print('Using model {0}'.format(config.get('modelfile')))
            read_model()
        validate_model()
    else:
        print('missing --data option')
        sys.exit(2)

    config.setdefault('lock','-')
    if(re.search('s',config.get('lock')) and re.search('t',config.get('lock')) and re.search('e',config.get('lock'))):
        print('No fitting done when locking all parameters, consider removing at least one or the fitting option.')
        sys.exit(2)
    if(config.get('fit')):
        try:
            default_rounds = 20
            config['rounds'] = int(config.setdefault('rounds',default_rounds))
        except ValueError:
            print('Erroneous rounds option ({0}), defaulting to {1}'.format(config.get('rounds'),default_rounds))
            config['rounds'] = default_rounds
        print('Fitting to {0} with {1} rounds. Output written with prefix "{2}"'.format(config.get('datafile'),config.get('rounds'),config.get('output')))
        [ print('{0} probabilities'.format('Locking '+x if re.search(x, config.get('lock')) else 'Updating '+x)) for x in ['s','t','e'] ]
    else:
        print('Disseqding {0}. Output written with prefix "{1}"'.format(config.get('datafile'),config.get('output')))

def validate_data():
    """Validates data before disseqtion."""
    global config
    obs = ''

    print('Reading data from {0}'.format(config.get('datafile')))
    fh = open_file(config.get('datafile'),'r')
    for line in fh:
        line = line.rstrip('\n')
        if(line[0:1] == '#'): # comments
            continue
        if(line[0:1] == '>'): # headers
            header = line[1:]
            continue
        obs += line
    fh.close()

    config.setdefault('obs',obs)
    config.setdefault('nobs',len(obs)-config.get('kmer')+1)
    print(config.get('obs'))

def read_model():
    """Reads model from file."""
    global config
    # global pi,A,B
    global literal_transition,literal_emission,numeric_transition,numeric_emission
    header = None
    pi = None
    A = None
    B = None
    literal_transition = None
    literal_emission = None
    numeric_transition = None
    numeric_emission = None

    print('Reading model from {0}'.format(config.get('modelfile')))
    fh = open_file(config.get('modelfile'),'r')
    for line in fh:
        line = line.rstrip('\n')
        if(line[0:1] == '#'): # comments
            continue
        if(line[0:1] == '>'): # headers
            header = line[1:]
            if(header == 'nstates'): nstates = list()
            if(header == 'start'): pi = list()
            if(header == 'transitions'): A = list(); literal_transition = list()
            if(header == 'emissions'): B = dict(); literal_emission = list()
            continue
        if(header == 'kmer' and not config.get('kmer')):
            config.setdefault('kmer',line)
            print('read kmer: {0}'.format(config.get('kmer')))
            continue
        if(header == 'nstates' and not config.get('nstates')):
            nstates += [line]
            print('read nstates: {0}'.format(line))
            continue
        if(header == 'alphabet' and not config.get('alphabet')):
            config.setdefault('alphabet',line)
            print('read alphabet: {0}'.format(line))
            continue
        if(header == 'start' and not config.get('start')):
            tmp = line.split('\t')
            pi += [tmp[0]]
            print('read start: {0}'.format(tmp[0]))
            continue
        if(header == 'transitions' and not config.get('transitions')):
            tmp = line.split('\t')
            A += [tmp]
            print('read transitions: {0}'.format(tmp))
            continue
        if(header == 'emissions' and not config.get('emissions')):
            tmp = line.split('\t')
            literal_emission += [tmp[0]]
            B[tmp[0]] = tmp[1:]
            print('read emissions: {0}: {1}'.format(tmp[0],tmp[1:]))
            continue
    fh.close()

    config.setdefault('nstates',','.join(nstates))
    config.setdefault('start',pi)
    config.setdefault('transitions',A)
    config.setdefault('emissions',B)

def validate_model():
    """Validates model parameters."""
    global config
    global pi,A,B
    global literal_state,literal_transition,literal_emission,numeric_state,numeric_transition,numeric_emission

    if(not config.get('kmer')):
        print('unable to disseqd: no option "kmer" given')
        sys.exit(2)
    else:
        try:
            config['kmer'] = int(config.get('kmer'))
            if(config.get('kmer')<1):
                raise ValueError
        except ValueError:
            print('unable to disseqd: option "kmer" must be a positive integer')
            sys.exit(2)

    if(not config.get('nstates')):
        print('unable to disseqd: no option "nstates" given')
        sys.exit(2)
    elif(isinstance(config.get('nstates'),str) and re.search(r',',config.get('nstates'))):
        config['nstates'] = config.get('nstates').split(',')
    if(isinstance(config.get('nstates'),list)):
        literal_state = dict(zip(numpy.arange(len(config.get('nstates'))),config.get('nstates')))
        numeric_state = dict(zip(config.get('nstates'),numpy.arange(len(config.get('nstates')))))
        config['nstates'] = len(config.get('nstates'))
    else:
        try:
            config['nstates'] = int(config.get('nstates'))
            if(config.get('nstates')<2):
                raise ValueError
            numeric_state = dict(zip(numpy.arange(config.get('nstates')),numpy.arange(config.get('nstates'))))
            literal_state = numeric_state
        except ValueError:
            print('unable to disseqd: option "nstates" must be a positive integer larger than 1')
            sys.exit(2)
    rng = numpy.arange(config.get('nstates'))

    if(not config.get('alphabet')):
        print('unable to disseqd: no option "alphabet" given')
        sys.exit(2)
    else:
        generate_literal_emissions(config.get('alphabet'),config.get('kmer'))
        numeric_emission = dict(zip(literal_emission,range(len(literal_emission))))

    pi = config.get('start')
    default_pi = list(itertools.repeat(1/config.get('nstates'),config.get('nstates')))
    if(not pi):
        print('no option "start" given, defaulting to uniform probabilities {0}'.format(1/config.get('nstates')))
        pi = default_pi
    elif(isinstance(pi,str) and re.search(r',',pi)):
        pi = pi.split(',')
    if(len(pi) == config.get('nstates')):
        try:
            pi = [ float(x) for x in pi ]
        except ValueError:
            print('option "start" must be {0} probabilities'.format(config.get('nstates')))
    else:
        print('no valid option "start" given, defaulting to uniform probabilities {0}'.format(1/config.get('nstates')))
        pi = default_pi
    pi = numpy.array(pi)
    config['start'] = pi

    A = config.get('transitions')
    if(not A):
        print('unable to disseqd: no option "transitions" given')
        sys.exit(2)
    elif(isinstance(A,str) and re.search(r',',A)):
        A = list([_] for _ in A.split(','))
    try:
        A = numpy.array(A,dtype=numpy.float64)
    except ValueError:
        print('unable to disseqd: transitions ({0}) must be numbers'.format(A))
        sys.exit(2)
    if(A.shape==(config.get('nstates'),1)):
        diag = [ 1-1/_ for _ in A[:,0] ]
        off_diag = [ (1-_)/(config.get('nstates')-1) for _ in diag ]
        A = numpy.tile(off_diag,(config.get('nstates'),1)).transpose()
        A[rng,rng] = diag
    if(not A.shape==(config.get('nstates'),config.get('nstates'))):
        print('unable to disseqd: mismatching number of states ({0}) and transitions ({1})'.format(config.get('nstates'),A))
        sys.exit(2)
    literal_transition = literal_state
    numeric_transition = numeric_state
    config['transitions'] = A

    B = config.get('emissions')
    if(not B):
        print('unable to disseqd: no option "emissions" given')
        sys.exit(2)
    elif(isinstance(B,str) and re.search(r',',B)):
        B = dict.fromkeys(list(B.split(',')),[])
    if(len(B.keys()) > 0):
        if(len(list(B.values())[0]) == config.get('nstates')):
            discard = [ x for x in B.keys() if not x in literal_emission ]
            print('discarding entries: ',discard)
            include = [ x for x in literal_emission if not x in B.keys() ]
            print('including entries: ',include)
            for x in discard:
                del B[x]
            for x in include:
                B[x] = list(itertools.repeat(sys.float_info.epsilon,config.get('nstates')))
            try:
                B = [ [ float(y) for y in B[x] ] for x in sorted(B) ]
            except ValueError:
                print('unable to disseqd: option "emission" must be {0} probabilities per emission'.format(config.get('nstates')))
                sys.exit(2)
            B = numpy.array(B)
        elif(len(B.keys()) == config.get('nstates')):
            files = list(itertools.chain(list(B.keys())))
            B = numpy.transpose(numpy.array([ read_emissions(f,literal_transition[t],config.get('kmer')) for f,t in zip(files,rng) ]))
        else:
            print('unable to disseqd: the different states need unique filenames')
            sys.exit(2)
    else:
        print('unable to disseqd: no vaild option "emissions" given')
        sys.exit(2)
    B = numpy.transpose(B)
    config['emissions'] = B

    pi = normalize(pi)
    A = normalize(A)
    B = normalize(B,desc='emissions')
    print('pi: {0}'.format(pi),pi.dtype)
    print('A: {0}'.format(A),A.dtype)
    print('B: {0}'.format(B),B.dtype)
    print('literal_transition: {0}'.format(literal_transition))
    print('literal_emission: {0}'.format(literal_emission))

def generate_literal_emissions(alphabet='ACGT',k=2):
    """Generates a list of possible emissions."""
    global literal_emission
    alphabet = ''.join(sorted(set(list(alphabet))))
    literal_emission = list(''.join(_) for _ in itertools.product(alphabet, repeat=k))

def read_emissions(f='',s=0, k=2):
    """Reads emissions from file."""
    global config
    global literal_emission
    emissions = dict.fromkeys(literal_emission,0)

    print('calculating {0}-mer emission probabilities for state {1} from file {2}'.format(k,s,f))

    fh = open_file(expanduser(f),'r')
    for line in fh:
        line = line.rstrip('\n')
        if(re.match(line,'>')):
            continue
        else:
            for kmer in range(len(line)-k+1):
                try:
                    emissions[line[kmer:kmer+k]]+=1
                except KeyError:
                    print('Skipping {0}; contains symbols not in alphabet "{1}"'.format(line[kmer:kmer+k],config.get('alphabet')))
    fh.close()
    return(list(v for k,v in sorted(emissions.items())))

def normalize(prob,desc=False):
    """Normalizes probabilities. Asserts minimum probability of emissions."""
    global config
    if(desc=='emissions'):
    # reshape needed for emission probabilities of >0-order HMMs (kmer>1)
        prob = prob.reshape((len(config.get('alphabet'))**(config.get('kmer')-1)*config.get('nstates'),len(config.get('alphabet'))))
    prob = (prob.T / prob.sum(axis=prob.ndim-1)).T
    # asserts no start or transition of only zeros
    prob[numpy.isnan(prob)] = 1/len(config.get('alphabet'))
    if(desc=='emissions'):
        # asserts a minimum probability to emissions
        minprob = 1/len(config.get('alphabet'))/1000
        marginal = prob.min(axis=prob.ndim-1) < minprob
        addprob = minprob/(1-minprob*len(config.get('alphabet')))
        prob[marginal] += addprob
        prob = (prob.T / prob.sum(axis=prob.ndim-1)).T
        prob = prob.reshape((config.get('nstates'),len(literal_emission)))
    return prob

def forward():
    """Performs forward-algorithm of update procedure."""
    global config
    global pi,A,B,c,f

    f = numpy.zeros(shape=(config.get('nobs')+1,config.get('nstates')),dtype=numpy.float64) # row vector
    c = numpy.ones(shape=(config.get('nobs')+1),dtype=numpy.float64)
    f[0,:] = pi
    # print(0,'-',c[0],f[0,:],sum(f[0,:]))
    symbol = config.get('obs')[0:config.get('kmer')]
    f[1,:] = f[0,:].dot(numpy.diag(B[:,numeric_emission[symbol]])) # no transition before first observation
    c[1] = sum(f[1,:])
    f[1,:] = f[1,:]/c[1]
    # print(t,symbol,c[t],f[t,:],sum(f[t,:]))
    for t in range(2,config.get('nobs')+1): #from 2 if first observation treated differently, but then sum of g[0,:] NOT equal to one
        symbol = config.get('obs')[t-1:t+config.get('kmer')-1]
        f[t,:] = f[t-1,:].dot(A).dot(numpy.diag(B[:,numeric_emission[symbol]]))
        c[t] = sum(f[t,:])
        f[t,:] = f[t,:]/c[t]
    #     print(t,symbol,c[t],f[t,:],sum(f[t,:]))
    # print(f)

def backward():
    """Performs bakcward-algorithm of update procedure."""
    global config
    global A,B,c,b

    b = numpy.ones(shape=(config.get('nstates'),config.get('nobs')+1),dtype=numpy.float64) # column vector
    # print(config.get('nobs'),'-',c[0],b[:,config.get('nobs')],sum(b[:,config.get('nobs')]))
    for t in range(config.get('nobs'),0,-1):
        symbol = config.get('obs')[t-1:t+config.get('kmer')-1]
        b[:,t-1] = A.dot(numpy.diag(B[:,numeric_emission[symbol]])).dot(b[:,t])
        b[:,t-1] = b[:,t-1]/c[t]
    #     print(t-1,symbol,c[t],b[:,t-1],sum(b[:,t-1]))
    # print(b.T)

def gamma():
    """Performs posterior-algorithm of update procedure."""
    global config
    global f,b,g

    g = numpy.zeros(shape=(config.get('nobs')+1,config.get('nstates')),dtype=numpy.float64) # row vector

    for t in range(0,config.get('nobs')+1):
        symbol = config.get('obs')[t-1:t+config.get('kmer')-1]
        g[t,:] = f[t,:]*b[:,t]
        # g[t,:] = g[t,:]/sum(g[t,:]) # if scaling f and b separately
        # print(t,symbol,c[t],g[t,:],sum(g[t,:]))
    # g = normalize(g) # if scaling f and b separately
    # print(g)

def xi():
    """Performs joint posterior-algorithm of update procedure."""
    global config
    global A,B,c,f,b,x

    x = numpy.zeros(shape=(config.get('nobs'),config.get('nstates'),config.get('nstates')),dtype=numpy.float64) # vector of matrices

    for t in range(0,config.get('nobs')):
        symbol = config.get('obs')[t:t+config.get('kmer')]
        x[t,:,:] = numpy.tile(f[t,:],(config.get('nstates'),1)).transpose()*numpy.tile(b[:,t+1],(config.get('nstates'),1))*A*B[:,numeric_emission[symbol]]/c[t+1]
        # x00 = f[t,0]*b[0,t+1]*A[0,0]*B[0,numeric_emission[symbol]]/c[t+1]
        # x01 = f[t,0]*b[1,t+1]*A[0,1]*B[1,numeric_emission[symbol]]/c[t+1]
        # x10 = f[t,1]*b[0,t+1]*A[1,0]*B[0,numeric_emission[symbol]]/c[t+1]
        # x11 = f[t,1]*b[1,t+1]*A[1,1]*B[1,numeric_emission[symbol]]/c[t+1]
        # xall = numpy.matrix([[x00,x01],[x10,x11]])
        # print(t,symbol,'\n',x[t,:,:],x[t,:,:].sum(),'\n',xall,xall.sum())
        # print(t,symbol,'\n',x[t,:,:],x[t,:,:].sum())
    # print(x)

def update():
    """Updates all probabilities unless other specified."""
    global config

    if(not re.search('s',config.get('lock'))): update_start()
    if(not re.search('t',config.get('lock'))): update_transmissions()
    if(not re.search('e',config.get('lock'))): update_emissions()
    # print('pi',pi,pi.sum(axis=0))
    # print('A',A,A.sum(axis=1))
    # print('B',B.transpose(),B.transpose().sum(axis=0))

def update_start():
    """Updates start probabilities."""
    global config
    global pi,g

    # pi = g[1,:]
    # pi needs normalization if no transition before first state
    pi = normalize(g[0,:])
    # print(x.sum(axis=0),'\n',numpy.tile(g.sum(axis=0),(2,1)).transpose())

def update_transmissions():
    """Updates transition probabilities."""
    global config
    global A,x,g

    # A does not need normalization beacause of scaling
    A = (x.sum(axis=0)-x[0,:,:])/numpy.tile((g.sum(axis=0)-g[0,:]-g[-1,:]),(config.get('nstates'),1)).transpose()

def update_emissions():
    """Updates emission probabilities."""
    global config
    global B,g

    B = numpy.zeros_like(B)
    for t in range(1,config.get('nobs')+1):
        symbol = config.get('obs')[t-1:t+config.get('kmer')-1]
        B[:,numeric_emission[symbol]] += g[t,:]
    # B needs normalization
    B = normalize((B.transpose()/g.sum(axis=0)).transpose(),desc='emissions')

def decode():
    """Decodes data with model parameters."""
    global config
    global pi,A,B,v,ll

    pil = numpy.nan_to_num(numpy.log(pi))
    Al = numpy.nan_to_num(numpy.log(A))
    Bl = numpy.nan_to_num(numpy.log(B))

    v = ['']*config.get('nstates')
    ll = numpy.zeros(shape=(config.get('nstates')),dtype=numpy.float64)

    symbol = config.get('obs')[0:config.get('kmer')]
    tmpl = pil+Bl[:,numeric_emission[symbol]]
    ll = tmpl
    v = [ _ + str(numpy.argmax(ll,axis=0)) for _ in v ]
    for t in range(1,config.get('nobs')):
        symbol = config.get('obs')[t:t+config.get('kmer')]
        tmpl = numpy.tile(ll,(config.get('nstates'),1)).transpose()+Al+Bl[:,numeric_emission[symbol]]
        ll = numpy.max(tmpl,axis=0)
        v = [ v[a] + str(b) for a,b in zip(numpy.argmax(tmpl,axis=0),range(config.get('nstates'))) ]

    v = v[numpy.argmax(ll)]
    ll = numpy.max(ll)

if __name__ == "__main__":
    main()
