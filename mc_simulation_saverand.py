import numpy,math,time, optparse, sys, pickle
import random as orig_random

class randwrap(object):
    def __init__(self, path):
        self.f = open(path, 'w')

    def _log(self, t, a, b, x):
        print >>self.f, t, a, b, x

    def choice(self, seq):
        l = len(seq)
        x = orig_random.randint(0, l - 1)
        self._log('I', 0, l - 1, x)
        return seq[x]

    def randint(self, a, b):
        x = orig_random.randint(a, b)
        self._log('I', a, b, x)
        return x

    def uniform(self, a, b):
        x = orig_random.uniform(a, b)
        self._log('D', a, b, x)
        return x



random = None
#accepted move vectors
#MOVES = numpy.array([[0,0,1], [0,1,0], [1,0,0], [0,0,-1], [0,-1,0], [-1,0,0], [1,0,1], [0,1,1], [1,1,0], [-1,0,-1], [0,-1,-1], [-1,-1,0], [1,0,-1], [0,1,-1], [1,-1,0], [-1,0,1], [0,-1,1], [-1,1,0], [-1,1,-1], [-1,-1,-1], [1,1,1], [1,-1,1], [1,-1,-1], [-1,-1,1], [-1,1,1], [1,1,-1]])
MOVES = numpy.array([[0,0,1], [0,1,0], [1,0,0], [0,0,-1], [0,-1,0], [-1,0,0], [1,0,1], [0,1,1], [1,1,0], [-1,0,-1], [0,-1,-1], [-1,-1,0], [1,0,-1], [0,1,-1], [1,-1,0], [-1,0,1], [0,-1,1], [-1,1,0]])
#accepted matching positions of binding sites
BMOVES = numpy.array([[1,0,0], [-1,0,0], [0,1,0], [0,-1,0], [0,0,1], [0,0,-1]])
#BMOVES = numpy.array([[1,0,0], [-1,0,0], [0,1,0], [0,-1,0], [0,0,1], [0,0,-1], [1,1,0], [-1,-1,0], [-1,1,0], [1,-1,0], [1,0,1], [-1,0,1], [1,0,-1], [-1, 0, -1], [0,1,-1], [0,-1,-1], [0,-1,1], [0,1,1]])

# radius of the nucleus
R = 10
# 2 x radius + a fringe, because lamin barrier has to be hermetic
BOUND = 2 * R + 2
ran_seed = 1
orig_random.seed(ran_seed)
print ran_seed

EMPTY = 0
BINDER = 1
LAMIN = 2
BSITE_R = 3
BSITE_L = 4
REGDNA = 5

DIST = 3

# max distance of good neighbors
GOOD_NEIGH = 3

def pars_inp():
    ## read/validate command-line arguments
    optparser = optparse.OptionParser(usage = "%prog regular_bsites.txt lamin_bsites.txt rand_out.txt [<options>]")

    optparser.add_option('-p', type = "string",
        dest = "Out_str",
        default = '',
        help = "An output filename with MC trajectory")
    optparser.add_option('-i', type = "string",
        dest = "In_str",
        default = '',
        help = "An input state (pickled file)")
    optparser.add_option('-s', type = "int",
        dest = "Steps",
        default = 100000,
        help = "Number of steps for simulation (default 100000)")
    optparser.add_option('-l', type = "int",
        dest = "Ch_lenght",
        default = 512,
        help = "Lenght of the chain (default 512)")
    optparser.add_option('-b', type = "int",
        dest = "Binders",
        default = 256,
        help = "Number of binders (default 256)")
			
    (opts, args) = optparser.parse_args()

    if len(args) < 3 and opts.In_str =='':
        print optparser.format_help() #prints help if no arguments
        sys.exit(1)
    return opts	

def init_dist_matrix(max_d = GOOD_NEIGH + 1):
    dist_matrix = numpy.zeros((max_d, max_d, max_d), dtype = numpy.float32)
    for i in range(max_d):
        for j in range(max_d):
            for k in range(max_d):
                dist_matrix[i][j][k] = math.sqrt(i**2 + j**2 + k**2)
    return dist_matrix
DIST_MATRIX = init_dist_matrix()

def initialize_import(f):
    list_ob = pickle.load(open(f))
    ch = list_ob[0]
    b = list_ob[1]
    a = list_ob[2]
    state = list_ob[3]
    return ch, b, a, state

def dist_from_mi(x, y, z, mi):
        return math.sqrt((x - mi)**2 + (y - mi)**2 + (z - mi)**2)

def getStateWithLamins(bound, f):
    state = numpy.zeros((bound, bound, bound), dtype = numpy.int)
    MIDDLE = bound / 2
    lam_name = f.split('.pdb')[0] + '_lamin.pdb'
    save_lam = open(lam_name, "w")
    save_lam.write("HEADER LAMINA")
    at_nr = 1

    for x in range(BOUND):
        for y in range(BOUND):
            for z in range(BOUND):
                border = abs(dist_from_mi(x, y, z, MIDDLE) - MIDDLE + 1)
                if border <= 2:
                    state[x, y, z] = LAMIN
                    if border == 1:
                        #print border
                        line = "\nATOM  " + str(at_nr).rjust(5) + " " + "P".center(4) + " " + "LAM" + "  " + str(at_nr).rjust(4) + "    " + str(round(x * DIST, 3)).rjust(8) + str(round(y * DIST, 3)).rjust(8) + str(round(z * DIST, 3)).rjust(8) + "  0.00 00.00"
                        at_nr += 1
                        save_lam.write(line)
                    
    save_lam.close()
    return state


def dist(r1,r2):
    return DIST_MATRIX[tuple(abs(r1 - r2))]
    

def cross(cha, po1, po2): # if are chain neighbours - return True
    # pier = numpy.where((cha == po1).all(axis=1))
    # dru = numpy.where((cha == po2).all(axis=1))
    pier = numpy.nonzero((cha == po1).all(axis=1))
    dru = numpy.nonzero((cha == po2).all(axis=1))
    if len(pier) == 1 and len(dru) == 1:
        if pier[0] == dru[0]+1 or pier[0] == dru[0]-1:
            return True
        else:
            return False
    else: print "Two atoms of the chain have the same coordinates. Please check ", po1, po2, pier[0], dru[0]
    
def intersect(new_p, next, sta, ch): #check if new position can cause intersection of the polimer
    if  dist(next, new_p) == 1: 
        return False
    elif dist(next, new_p) > 1:
        differ = new_p - next
        if differ[0] == 0:
            pos1 = numpy.array(([next[0], next[1], new_p[2]]))
            pos2 = numpy.array(([next[0], new_p[1], next[2]]))
            if sta[tuple(pos1)] not in [EMPTY, LAMIN, BINDER] and sta[tuple(pos2)] not in [EMPTY, LAMIN, BINDER]:
                return cross(ch, pos1, pos2)
            else: return False         
        elif differ[1] == 0:
            pos1 = numpy.array(([next[0], next[1], new_p[2]]))
            pos2 = numpy.array(([new_p[0], next[1], next[2]]))
            if sta[tuple(pos1)] not in [EMPTY, LAMIN, BINDER] and sta[tuple(pos2)] not in [EMPTY, LAMIN, BINDER]:
                return cross(ch, pos1, pos2)
            else: return False
        elif differ[2] == 0:
            pos1 = numpy.array(([next[0], new_p[1], next[2]]))
            pos2 = numpy.array(([new_p[0], next[1], next[2]]))
            if sta[tuple(pos1)] not in [EMPTY, LAMIN, BINDER] and sta[tuple(pos2)] not in [EMPTY, LAMIN, BINDER]:
                return cross(ch, pos1, pos2)
            else: return False
        else: print new_p, next, "The distance between these two positions are not = sqrt(2)"

def initialize_random(n, m, fa, bound = BOUND):
    chain = numpy.zeros((n, 3), dtype = numpy.int)
    binders = numpy.zeros((m, 3), dtype = numpy.int)
    state = getStateWithLamins(bound, fa)
    attached_to_lamins = []

    chain[0] = [bound / 2] * 3

    def get_site_type_list(fpath, length):
        positions = [0] * length
        for l in open(fpath):
            positions[int(l) -1] = 1
        return positions

    regular_bsites = get_site_type_list(sys.argv[1], n)
    lamin_bsites   = get_site_type_list(sys.argv[2], n)

    def get_site_type(i, regular_bsites, lamin_bsites): # BSITE_R interacts with binders whereas BSITE_L interacts both with lamins and binders
        if regular_bsites[i] == 1 and lamin_bsites[i] == 1:
            print "The lamin site are the same as regular! Please change it and rerun the program"
            sys.exit(1)
        elif regular_bsites[i] == 1:
            return BSITE_R
        elif lamin_bsites[i] == 1:
            return BSITE_L
        else:
            return REGDNA

    cur = chain[0]
    state[tuple(cur)] = get_site_type(0, regular_bsites, lamin_bsites)
        
    for i in range(1, n):
        mov = random.choice(MOVES)
        tries = 0
        while tries < 100 and (not (no_collisions(tuple(cur + mov), state)) or intersect(cur, cur+mov, state, chain)):
            mov = random.choice(MOVES)
            tries += 1    
            
        assert tries != 100, "unable to find initialization"
        chain[i] = cur + mov
        state[tuple(chain[i])] = get_site_type(i, regular_bsites, lamin_bsites)

        if state[tuple(chain[i])] == BSITE_L and count_bonds(chain[i], [LAMIN], state) > 0:
            attached_to_lamins.append(tuple(chain[i]))

        cur = chain[i]
        
    mid = bound/2
    for i in range(m):
        x = random.randint(0, bound)
        y = random.randint(0, bound)
        z = random.randint(0, bound)
        tries = 0
        distance = dist_from_mi(x, y, z, mid)
        while distance > mid-3 or (not (no_collisions((x, y, z), state)) and tries < 100):
            x = random.randint(0, bound)
            y = random.randint(0, bound)
            z = random.randint(0, bound)
            tries += 1
            distance = dist_from_mi(x, y, z, mid)
        binders[i] = [x, y, z]
        state[tuple(binders[i])] = BINDER
    

    return chain, binders, attached_to_lamins, state


def good_neighbors(x, i, chain):
    #check neighbors (chain)
    if i > 0:
        d1 = dist(chain[i - 1], x)
        if d1 > 0 and d1 < GOOD_NEIGH:
            pass
        else:
            return False
    if i < len(chain) - 1:
        d2 = dist(chain[i + 1], x)
        if d2 > 0 and d2 < GOOD_NEIGH:
            pass
        else:
            return False
    return True
    
def no_collisions(x, state):
    if state[tuple(x)] != EMPTY:
        return False
    else:
        return True

def bonds(chain, stat):

    bonds = 0
    for j in range(chain.shape[0]):
        molecule_pos = chain[j]
        molecule = stat[tuple(molecule_pos)]
        #print j
        if molecule == REGDNA:
            continue
 
        elif molecule == BSITE_R:
            binding = [BINDER]
        elif molecule == BSITE_L:
            binding = [LAMIN, BINDER]
        else: print "INNY STAN!!!", j, molecule
        #print j, molecule
        
        one_ch_at = 0
        for bmove in BMOVES: 
            new = molecule_pos + bmove
            enc = stat[tuple(new)]
            if enc in binding:
                bonds += 1

    return bonds

def modify(chain, binders, state, bound = BOUND):
    #move binders
    if random.randint(0, 1):
        i = random.randint(0, len(binders) - 1)
        move = random.choice(MOVES)
        new = move + binders[i]

        if no_collisions(tuple(new), state):
            return False, i, move
    #move residues
    else:
        i = random.randint(0, len(chain) - 1)
        move = random.choice(MOVES)
        new = move + chain[i]
        if good_neighbors(new, i, chain) and no_collisions(tuple(new), state):  # test if there is no collisions (the same place by different atoms) and no intersect of bonds
            if i != len(chain) - 1:
                if dist(chain[numpy.absolute(i-1)], new) <= numpy.sqrt(2) and dist(chain[numpy.absolute(i+1)], new) <= numpy.sqrt(2) and not intersect(new, chain[numpy.absolute(i-1)], state, chain) and not intersect(new, chain[numpy.absolute(i+1)], state, chain):
                    return True, i, move
                else: pass
                         
            elif dist(chain[numpy.absolute(i-1)], new) <= numpy.sqrt(2) and not intersect(new, chain[numpy.absolute(i-1)], state, chain):
                return True, i, move
                
            else:
                pass
        else:
            pass 
    return None

DIST = 3
def write_as_pdb(chain, binders, attached_to_lamins, state, f, nb, metr_step, name = "chromosome and binders"):

    l = chain.shape[0]
    n = binders.shape[0]
    f.write("HEADER %i %d step %i\nTITLE %s" % (nb, l + n, metr_step, name))
    at_nr = 0

    def pdb_line(at_name, at_nr, desc, pos):
        return "\nATOM  " + str(at_nr).rjust(5) + " " + at_name.center(4) + " " + desc + "  " + str(at_nr).rjust(4) + "    " + str(round(pos[0] * DIST, 3)).rjust(8) + str(round(pos[1] * DIST, 3)).rjust(8) + str(round(pos[2] * DIST, 3)).rjust(8) + "  0.00 00.00"

    for i in range(l):
        at_n = 'C'
        cur_chain = chain[i]
        if state[tuple(cur_chain)] == REGDNA:
            r = "UNB"
        elif state[tuple(cur_chain)] == BSITE_R:
            r = "BOU"
        else: 
            r = "LAM" #BSITE_L
        at_nr += 1
        f.write(pdb_line(at_n, at_nr, r, chain[i]))

    chain_at = at_nr

    for i in range(n):
        at_nr += 1
        r = "BIN"
        at_n = 'O'
        f.write(pdb_line(at_n, at_nr, r, binders[i]))


    for i in range(1, chain_at - 1):
        line = "\nCONECT" + str(i).rjust(5) +  str(i + 1).rjust(5)
        f.write(line)
    f.write("\nEND   \n")

def count_bonds(pos, accepted, state):
    bonds = 0
    for bmove in BMOVES:
        if state[tuple(pos + bmove)] in accepted:
             bonds += 1
    #print bonds
    return bonds

def radius_gyr(chai):
    length = chai.shape[0]
    r = 0.0
    for at in range(length):
        for ato in range(at):
            r = r + math.sqrt(numpy.sum((chai[at] - chai[ato])**2))
    r_gyr = 2*r/(length*(length-1))
    return r_gyr


DELTA=2
# GYRATION = True
GYRATION = False
CHECK_E = False
def metropolis(chain, binders, attached_to_lamins, state, out_fname, name = "chromosome", n = 100):

    out_file = open(out_fname, "w")
    st_nr = 0
    E = bonds(chain, state)
    print "Starting energy:", E
    write_as_pdb(chain, binders, attached_to_lamins, state, out_file, st_nr, 0, name + ";bonds=" + str(E))

    for step in range(n):

        resp = None
        while not resp:
            resp = modify(chain, binders, state)

        ch = numpy.array(chain, copy = True)
        b = numpy.array(binders, copy = True)
        if CHECK_E:
            st = numpy.array(state,copy=True)
        else: pass
        if resp:
            in_chain, i, move = resp
            if in_chain:
                #print i, "change"
                old = numpy.copy(ch[i])
                ch[i] = ch[i] + move
                if CHECK_E:
                    st[tuple(ch[i])] = st[tuple(old)]
                    st[tuple(old)] = EMPTY
                else: pass
                if state[tuple(old)] == BSITE_R:
                    Enew = E + count_bonds(ch[i], [BINDER], state) - count_bonds(old, [BINDER], state)
                    if CHECK_E:
                        Eslow = bonds(ch, st)
                        if Enew != Eslow:
                            print "R", 'Enew', Enew, 'Eslow', Eslow
                    else: pass
                elif state[tuple(old)] == BSITE_L:
                    Enew = E + count_bonds(ch[i], [LAMIN, BINDER], state) - count_bonds(old, [LAMIN, BINDER], state)
                    if CHECK_E:
                        Eslow = bonds(ch, st)
                        if Enew != Eslow: 
                            print "L", 'Enew', Enew, 'Eslow', Eslow
                    else: pass
                    if tuple(ch[i]) not in attached_to_lamins and count_bonds(ch[i], [LAMIN], state) > 0:
                        attached_to_lamins.append(tuple(ch[i]))
                    elif tuple(ch[i]) in attached_to_lamins and count_bonds(ch[i], [LAMIN], state) == 0:
                        attached_to_lamins.remove(ch[i])

                    else:
                        pass 


                else: # REGDNA
                    Enew = E
            else:
                old = numpy.copy(b[i])
                b[i] = b[i] + move
                if CHECK_E:
                    st[tuple(b[i])] = st[tuple(old)]
                    st[tuple(old)] = EMPTY
                else: pass
                Enew = E + count_bonds(b[i], [BSITE_R, BSITE_L], state) - count_bonds(old, [BSITE_R, BSITE_L], state) 
                if CHECK_E:
                    Eslow = bonds(ch, st)
                    if Enew != Eslow:
                        print "B", 'Enew', Enew, 'Eslow', Eslow
                else: pass
        else:
            Enew = E

        if Enew >= E or random.uniform(0.0, 1.0) < math.exp((Enew - E) * DELTA): #accept

            chain = ch
            binders = b
            #write_as_pdb(chain, binders, attached_to_lamins, state, out_file, st_nr, name + ";bonds=" + str(E))
            if resp:
                state[tuple(old + move)] = state[tuple(old)]
                state[tuple(old)] = EMPTY

            if E != Enew:

                E = Enew
                st_nr += 1
                if GYRATION:
                     print "iter", step, "step", st_nr, "energy:", E, "R_gyr ", radius_gyr(chain)
                else:
                    print "iter", step, "step", st_nr, "energy:", E
                
                write_as_pdb(chain, binders, attached_to_lamins, state, out_file, st_nr, step, name + ";bonds=" + str(E))
                #print "WRITE!!!"

    # dump the last state to the pickle
    l_obj = [chain, binders, attached_to_lamins, state]
    pickle_fname = out_fname.split('.pdb')[0] + ".pick"
    pickle_file = open(pickle_fname, 'w')
    pickle.dump(l_obj, pickle_file)

    out_file.close()


def output_name(ou, m, n):
    if ou == '':
        f_n = "MC_traj_%ibin_%ichain.pdb" % (m, n)
    elif "." in ou:
        f_n = ou.split('.')[0] + ".pdb"
    else:
        f_n = ou + ".pdb"
    return f_n

opts = pars_inp()
random = randwrap(sys.argv[3])

if opts.In_str == '':
    rand_init = True
else: 
    rand_init = False


t1 = time.time()
if rand_init: 
    N = opts.Ch_lenght # length of the chain
    M = opts.Binders # nr of binders
    fn = output_name(opts.Out_str, M, N)
    c, b, a, state = initialize_random(N, M, fn)
else: 
    c, b, a, state = initialize_import(opts.In_str)
    N = c.shape[0]
    M = b.shape[0]
    fn = output_name(opts.Out_str, M, N)

print "The lenght of the chain is ", N, ", the number of binders is ", M, ", the nucleus radius is ", R, ", the number of steps is ", opts.Steps, "random_seed ", ran_seed

t2 = time.time()
print "initialization: ", t2 - t1
BOUND = numpy.max(c)

a = []
metropolis(c, b, a, state, fn, n = opts.Steps)
t1 = t2
t2 = time.time()
print "metropolis: ", t2 - t1
