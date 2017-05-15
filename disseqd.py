import argparse
import sys
import numpy
import itertools
import re
from os.path import expanduser

def main():
    global config

    parser = argparse.ArgumentParser(description='Detection of Interspersed Signatures in SEQuence Data')

    basic = parser.add_argument_group('basic arguments')
    basic.add_argument('-d','--data' ,help='file with input fasta sequence')
    basic.add_argument('-m','--model', help='file with model parameters')

    output = parser.add_argument_group('output arguments')
    output.add_argument('-v', '--verbose', action='store_true', help='increases verbosity level')
    output.add_argument('-o', '--output', metavar='PREFIX', help='outputs to files PREFIX.model.txt and PREFIX.decode.txt')

    train = parser.add_argument_group('training arguments')
    train.add_argument('-f', '--fit', metavar='STE', help='fits any combination of start (S), transmission (T), and emission (E) probabilities')
    train.add_argument('-r', '--rounds', type=int, help='fits model for ROUNDS number of rounds')

    override = parser.add_argument_group('override arguments')
    override.add_argument('-k', '--kmer', type=int,  help='order of model plus 1')
    override.add_argument('-n', '--nstates', nargs='+',  help='number or names of states')
    override.add_argument('-a', '--alphabet',  help='alphabet used for model')
    override.add_argument('-s', '--start', nargs='+', metavar='PROB', help='NSTATES start probabilities')
    override.add_argument('-t', '--transitions', nargs='+', metavar='LENGTH',  help='NSTATES expected lenghts of states')
    override.add_argument('-e', '--emissions', nargs='+', metavar='FILE',  help='NSTATES files with state observations')

    args = parser.parse_args()
    config = vars(args)
    format_config()
    config.get('verbose') and print(config.items())

    validate_input()
    validate_data()

    if(config.get('fit')):
        for r in range(config.get('rounds')):
            config.get('verbose') and print('Fitting round {0}'.format(r))
            forward()
            backward()
            gamma()
            xi()
            update()
            decode()
            config.get('verbose') and print(v,ll)
    else:
        decode()
    print(v,ll)
    output_model()

def format_config():
    global emission_entries

    if(isinstance(config.get('transitions'),list)): config['transitions'] = [ [x] for x in config.get('transitions') ]
    if(isinstance(config.get('emissions'),list)): emission_entries = config.get('emissions')

def read_model():
    """Reads model from file."""
    global config
    global emission_entries

    nstates = None
    pi = None
    A = None
    B = None
    emission = None

    config.get('verbose') and print('Reading model from {0}'.format(config.get('model')))
    fh = open_file(config.get('model'),'r')
    for line in fh:
        line = line.rstrip('\n')
        if(line[0:1] == '#'): # comments
            continue
        if(line[0:1] == '>'): # headers
            header = line[1:]
            if(header == 'nstates'): nstates = list()
            if(header == 'start'): pi = list()
            if(header == 'transitions'): A = list();
            if(header == 'emissions'): B = list(); emission = list()
            continue
        if(header == 'kmer' and config.get('kmer') is None):
            config['kmer'] = line
            config.get('verbose') and print('read kmer: {0}'.format(line))
            continue
        if(header == 'nstates' and config.get('nstates') is None):
            nstates += [line]
            config.get('verbose') and print('read nstates: {0}'.format(line))
            continue
        if(header == 'alphabet' and config.get('alphabet') is None):
            config['alphabet'] = line
            config.get('verbose') and print('read alphabet: {0}'.format(line))
            continue
        if(header == 'start' and config.get('start') is None):
            tmp = line.split('\t')
            pi += [tmp[0]]
            config.get('verbose') and print('read start: {0}'.format(tmp[0]))
            continue
        if(header == 'transitions' and config.get('transitions') is None):
            tmp = line.split('\t')
            A += [tmp]
            config.get('verbose') and print('read transitions: {0}'.format(tmp))
            continue
        if(header == 'emissions' and config.get('emissions') is None):
            tmp = line.split('\t')
            emission += [tmp[0]]
            B += [tmp[1:]]
            config.get('verbose') and print('read emissions: {0}: {1}'.format(tmp[0],tmp[1:]))
            continue
    fh.close()

    if(config.get('nstates') is None): config['nstates'] = nstates
    if(config.get('start') is None): config['start'] = pi
    if(config.get('transitions') is None): config['transitions'] = A
    if(config.get('emissions') is None):
        config['emissions'] = B
        emission_entries = emission

def validate_input():
    """Validates the input."""
    global config

    if(config.get('output') == None):
        config['output'] = 'out'
        print('No option "output" specified, using "{0}".'.format(config.get('output')))

    if(config.get('data')):
        config.get('verbose') and print('Using data in {0}'.format(config.get('data')))
        if(config.get('model')):
            config.get('verbose') and print('Using model in {0}'.format(config.get('model')))
            read_model()
        validate_model()
    else:
        print('missing "data" option')
        sys.exit(2)

    if(config.get('fit')):
        if(re.search('(?i)[ste]',config.get('fit'))):
            try:
                default_rounds = 20
                if(config.get('rounds') is None): config['rounds'] = default_rounds
                config['rounds'] = int(config.get('rounds'))
            except ValueError:
                print('Erroneous "rounds" option ({0}), defaulting to {1}'.format(config.get('rounds'),default_rounds))
                config['rounds'] = default_rounds
                config.get('verbose') and print('Fitting to {0} with {1} rounds. Output written with prefix "{2}"'.format(config.get('data'),config.get('rounds'),config.get('output')))
            [ config.get('verbose') and print('{0} probabilities.'.format('Updating '+x if re.search('(?i)'+x, config.get('fit')) else 'Locking '+x)) for x in ['S','T','E'] ]
        else:
            print('No fitting done when specifying none of S/T/E parameters.')
            sys.exit(2)
    else:
        config.get('verbose') and print('Disseqding {0}. Output written with prefix "{1}".'.format(config.get('data'),config.get('output')))

def validate_model():
    """Validates model parameters."""

    validate_kmer()
    validate_nstates()
    validate_alphabet()
    validate_start()
    validate_transitions()
    validate_emissions()

def validate_kmer():
    """Validates kmer parameter."""

    # kmer is None, integer or misspecified
    if(config.get('kmer') is None): # None
        print('Unable to disseqd: no option "kmer" given.')
        sys.exit(2)
    else: # integer
        try:
            config['kmer'] = int(config.get('kmer'))
            if(config.get('kmer')<1): raise ValueError
        except ValueError: # misspecified
            print('Unable to disseqd: option "kmer" must be a positive integer, not {0}.'.format(config.get('kmer')))
            sys.exit(2)

    config.get('verbose') and print('kmer: {0}'.format(config.get('kmer')))

def validate_nstates():
    """Validates nstates parameter and generates state name lookups."""
    global literal_state,numeric_state

    # nstates is None, number, list with number, list with names
    if(config.get('nstates') is None): # None
        print('Unable to disseqd: no option "nstates" given.')
        sys.exit(2)
    if(not isinstance(config.get('nstates'),list)): # number
        config['nstates'] = [ config.get('nstates') ]
    if(len(config.get('nstates'))==1): # list with number
        try:
            config['nstates'] = int(config.get('nstates')[0])
            if(config.get('nstates')<2): raise ValueError
            numeric_state = dict(zip(numpy.arange(config.get('nstates')),numpy.arange(config.get('nstates'))))
            literal_state = numeric_state
        except ValueError:
            print('Unable to disseqd: option "nstates" must be a positive integer larger than 1.')
            sys.exit(2)
    else: # list with names
        numeric_state = dict(zip(config.get('nstates'),numpy.arange(len(config.get('nstates')))))
        literal_state = dict(zip(numpy.arange(len(config.get('nstates'))),config.get('nstates')))
        config['nstates'] = len(config.get('nstates'))

    config.get('verbose') and print('nstates: {0}'.format(config.get('nstates')))

def validate_alphabet():
    """Validates alphabet string and generates list of emissions from it."""
    global numeric_emission

    # alphabet is None or string
    if(config.get('alphabet') is None): # None
        print('Unable to disseqd: no option "alphabet" given.')
        sys.exit(2)
    else:
        generate_literal_emissions(config.get('alphabet'),config.get('kmer'))
        numeric_emission = dict(zip(literal_emission,range(len(literal_emission))))

    config.get('verbose') and print('alphabet: {0}'.format(config.get('alphabet')))

def generate_literal_emissions(alphabet='ACGT',k=2):
    """Generates a list of possible emissions."""
    global literal_emission

    alphabet = ''.join(sorted(set(list(alphabet))))
    literal_emission = list(''.join(x) for x in itertools.product(alphabet, repeat=k))

def validate_start():
    """Validates start parameters."""

    # start is None or list of probabilities or misspecified
    pi = config.get('start')
    default_pi = list(itertools.repeat(1/config.get('nstates'),config.get('nstates')))
    if(pi is None): # None
        config.get('verbose') and print('No option "start" given, defaulting to uniform probabilities {0}.'.format(1/config.get('nstates')))
        pi = default_pi
    if(len(pi) == 1 or len(pi) == config.get('nstates')): # list of probabilities
        try:
            pi = [ float(x) for x in pi ]
        except ValueError:
            print('option "start" must be probabilities, not {1}.'.format(pi))
    else: # misspecified
        config.get('verbose') and print('Option "start" must be probabilities, not {0}, defaulting to uniform probabilities {1}.'.format(pi,1/config.get('nstates')))
        pi = default_pi
    pi *= config.get('nstates')
    pi = numpy.array(pi[0:config.get('nstates')])
    # pi = numpy.array(pi)
    config['start'] = normalize(pi)

    config.get('verbose') and print('start: {0}'.format(config.get('start')))

def validate_transitions():
    """Validates transition probabilities and generates transition lookups."""

    # transitions is None, list of lengths, list of list of probabilities or misspecified
    A = config.get('transitions')
    if(A is None or not isinstance(A,list)): # None
        print('Unable to disseqd: no option "transitions" given.')
        sys.exit(2)
    try: # list
        A = numpy.array(A,dtype=numpy.float64)
    except ValueError:
        print('Unable to disseqd: transitions ({0}) must be numbers'.format(A))
        sys.exit(2)
    # print(A.shape,A.size,A)
    if(A.shape==(config.get('nstates'),1)): # lengths
        diag = [ 1-1/i for i in A[:,0] ]
        off_diag = [ (1-i)/(config.get('nstates')-1) for i in diag ]
        A = numpy.tile(off_diag,(config.get('nstates'),1)).transpose()
        rng = numpy.arange(config.get('nstates'))
        A[rng,rng] = diag
    if(not A.shape==(config.get('nstates'),config.get('nstates'))): # misspecified
        print('Unable to disseqd: mismatching number of states ({0}) and transitions ({1})'.format(config.get('nstates'),A))
        sys.exit(2)
    config['transitions'] = A

    config.get('verbose') and print('transitions: {0}'.format(config.get('transitions')))

def validate_emissions():
    """Validates emission probabilities."""
    global emission_entries

    #emissions is None, list of files, list of list of probabilities or misspecified
    B = config.get('emissions')
    if(B is None or not isinstance(B,list)): # None
        print('Unable to disseqd: no option "emissions" given.')
        sys.exit(2)
    if(emission_entries): # defined/exists
        if(len(B[0]) == config.get('nstates')): # list of list of probabilities
            tmp = dict(zip(emission_entries,B))
            discard = [ x for x in emission_entries if not x in literal_emission ]
            config.get('verbose') and print('discarding entries: ',discard)
            include = [ x for x in literal_emission if not x in emission_entries ]
            config.get('verbose') and print('including entries: ',include)
            for x in discard:
                del tmp[x]
            for x in include:
                tmp[x] = list(itertools.repeat(sys.float_info.epsilon,config.get('nstates')))
            try:
                tmp = [ [ float(y) for y in tmp[x] ] for x in sorted(tmp) ]
            except ValueError:
                print('Unable to disseqd: option "emissions" must be {0} probabilities per emission.'.format(config.get('nstates')))
                sys.exit(2)
            B = numpy.array(tmp)
        elif(len(emission_entries) == config.get('nstates')): # list of file names
            rng = numpy.arange(config.get('nstates'))
            B = numpy.transpose(numpy.array([ read_emissions(f,literal_state[t],config.get('kmer')) for f,t in zip(emission_entries,rng) ]))
        else:
            print('Unable to disseqd: please supply {0} unique filenames, not {1}.'.format(config.get('nstates'),emission_entries))
            sys.exit(2)
    else: # undefined
        print('Unable to disseqd: no vaild option "emissions" given')
        sys.exit(2)
    B = numpy.transpose(B)
    B = normalize(B,'emissions')
    config['emissions'] = B

    config.get('verbose') and print('emissions: {0}'.format(config.get('emissions')))

def validate_data():
    """Validates data before disseqtion."""
    global config

    config.get('verbose') and print('Reading data from {0}'.format(config.get('data')))

    obs = ''
    fh = open_file(config.get('data'),'r')
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
    # print(config.get('obs'))

def forward():
    """Performs forward-algorithm of update procedure."""
    global config
    global c,f

    pi = config.get('start')
    A = config.get('transitions')
    B = config.get('emissions')

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
    global c,b

    A = config.get('transitions')
    B = config.get('emissions')

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
    global c,f,b,x

    A = config.get('transitions')
    B = config.get('emissions')

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

    if(re.search('(?i)S',config.get('fit'))): update_start()
    if(re.search('(?i)T',config.get('fit'))): update_transmissions()
    if(re.search('(?i)E',config.get('fit'))): update_emissions()
    # print('pi',pi,pi.sum(axis=0))
    # print('A',A,A.sum(axis=1))
    # print('B',B.transpose(),B.transpose().sum(axis=0))

def update_start():
    """Updates start probabilities."""
    global config
    global g

    # pi = g[1,:]
    # pi needs normalization if no transition before first state
    config['start'] = normalize(g[0,:])
    # print(g.sum(axis=0),'\n',numpy.tile(g.sum(axis=0),(2,1)).transpose())

def update_transmissions():
    """Updates transition probabilities."""
    global config
    global x,g

    # A does not need normalization beacause of scaling
    config['transitions'] = (x.sum(axis=0)-x[0,:,:])/numpy.tile((g.sum(axis=0)-g[0,:]-g[-1,:]),(config.get('nstates'),1)).transpose()

def update_emissions():
    """Updates emission probabilities."""
    global config
    global g

    B = numpy.zeros_like(config.get('emissions'))
    for t in range(1,config.get('nobs')+1):
        symbol = config.get('obs')[t-1:t+config.get('kmer')-1]
        B[:,numeric_emission[symbol]] += g[t,:]
    # B needs normalization
    B = normalize((B.transpose()/g.sum(axis=0)).transpose(),desc='emissions')
    config['emissions'] = B

def decode():
    """Decodes data with model parameters."""
    global config
    global v,ll

    pil = numpy.nan_to_num(numpy.log(config.get('start')))
    Al = numpy.nan_to_num(numpy.log(config.get('transitions')))
    Bl = numpy.nan_to_num(numpy.log(config.get('emissions')))
    v = ['']*config.get('nstates')
    ll = numpy.zeros(shape=(config.get('nstates')),dtype=numpy.float64)

    symbol = config.get('obs')[0:config.get('kmer')]
    tmpl = pil+Bl[:,numeric_emission[symbol]]
    ll = tmpl
    v = [ x + str(numpy.argmax(ll,axis=0)) for x in v ]
    for t in range(1,config.get('nobs')):
        symbol = config.get('obs')[t:t+config.get('kmer')]
        tmpl = numpy.tile(ll,(config.get('nstates'),1)).transpose()+Al+Bl[:,numeric_emission[symbol]]
        ll = numpy.max(tmpl,axis=0)
        v = [ v[a] + str(b) for a,b in zip(numpy.argmax(tmpl,axis=0),range(config.get('nstates'))) ]

    v = v[numpy.argmax(ll)]
    ll = numpy.max(ll)

def normalize(prob,desc=''):
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

def read_emissions(f='',s=0, k=2):
    """Reads emissions from file."""
    global config
    global literal_emission

    config.get('verbose') and print('calculating {0}-mer emission probabilities for state {1} from file {2}'.format(k,s,f))

    emissions = dict.fromkeys(literal_emission,0)
    fh = open_file(expanduser(f),'r')
    for line in fh:
        line = line.rstrip('\n')
        if(re.match('>',line)):
            continue
        else:
            for kmer in range(len(line)-k+1):
                try:
                    emissions[line[kmer:kmer+k]]+=1
                except KeyError:
                    config.get('verbose') and print('Skipping {0}; contains symbols not in alphabet "{1}"'.format(line[kmer:kmer+k],config.get('alphabet')))
    fh.close()
    return(list(v for k,v in sorted(emissions.items())))

def output_model():
    """Writes the model to file."""
    global config
    global literal_emission

    f = str(config.get('output'))+'.model.txt'
    config.get('verbose') and print('Writing model to file {0}.'.format(f))
    fh = open_file(f,'w')
    fh.write('>kmer\n{0}\n'.format(config.get('kmer')))
    fh.write('>nstates\n{0}\n'.format('\n'.join(str(literal) for numeric,literal in sorted(literal_state.items()))))
    fh.write('>alphabet\n{0}\n'.format(config.get('alphabet')))
    fh.write('>start\n{0}\n'.format('\n'.join([ str(i) for i in config.get('start') ])))
    fh.write('>transitions\n')
    [ fh.write('{0}\n'.format('\t'.join(str(prob) for prob in row))) for row in config.get('transitions') ]
    fh.write('>emissions\n')
    [ fh.write('{0}\t{1}\n'.format(label,'\t'.join(str(prob) for prob in row))) for label,row in zip(literal_emission,config.get('emissions').transpose()) ]
    fh.write('\n')
    fh.close()

    f = str(config.get('output'))+'.decode.txt'
    config.get('verbose') and print('Writing decoded data to file {0}.'.format(f))
    fh = open_file(f,'w')
    fh.write('>loglikelihood:'+str(ll))
    [ fh.write('\t{0}:{1}'.format(literal,numeric)) for numeric,literal in sorted(literal_state.items()) ]
    fh.write('\n'+v+'\n')
    fh.close()

def open_file(f,mode='r'):
    """Returns a filehandle to the file, opened in specified mode."""
    try:
        fh = open(expanduser(f),mode)
    except OSError as err:
        print('{0} {1}'.format(f,err.strerror))
        sys.exit(2)
    return fh

if __name__ == "__main__":
    main()
