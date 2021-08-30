import get_data
import numpy
import os


def get_distance(c1, c2):
    """
    Calculates the distance between two coordinate points
    :param c1: array of x,y,z
    :param c2: array of x,y,z
    :return: distance between two ponts
    """
    return numpy.linalg.norm(c1 - c2)


def get_centroid(p):
    # print (len(p),p)
    x = [i[0] for i in p]
    y = [i[1] for i in p]
    z = [i[2] for i in p]
    c = [sum(x) / len(p), sum(y) / len(p), sum(z) / len(p)]
    return numpy.array(c)


def find_angle( p, c, d):
    pc,p = p
    v1 = p[1] - pc
    v2 = p[2] - pc
    nv = numpy.cross(v1, v2)
    #nv2 = numpy.cross(v2, v1)
    cv = c - pc
    #cv2 = c - cn
    nnv = nv / numpy.linalg.norm(nv)
    ncv = cv / numpy.linalg.norm(cv)
    #ncv2 = cv2 / numpy.linalg.norm(cv2)
    dp = abs(numpy.dot(nnv, ncv))
    #dp2 = abs(numpy.dot(nnv, ncv2))
    ang = numpy.arccos(dp)
    #ang2 = numpy.arccos(dp2)
    ang_deg = (180 / numpy.pi) * ang
    #ang_deg2 = (180 / numpy.pi) * ang2
    # print(ang_deg)
    s_ang = solid_angle(ang_deg, d)
    return ang_deg, s_ang


def solid_angle(a_deg, r):
    s = 1.4 #C-C bond length
    A=((3.0*numpy.sqrt(3))/2.0)*s*s #area of the hexagonal plane
    a = (numpy.pi / 180) * a_deg
    # sa2=2*numpy.pi*(1.0-1.0/(numpy.sqrt(1+(A*numpy.cos(a)/(numpy.pi*r1**r1)))))
    sa = 2 * numpy.pi * (1.0 - 1.0 / (numpy.sqrt(1 + (A*numpy.cos(a) / (numpy.pi * r * r)))))
    # print (a_deg)
    sa_deg = (180.0 / numpy.pi) * sa
    # sa_deg2 = (180.0 / numpy.pi) * sa2
    # print (a_deg,sa_deg,sa_deg2)
    return sa_deg

def get_aromatic_info(pdb_data):
    aromatic_atoms = {
        'PHE': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'HD1', 'HD2', 'HE1', 'HE2', 'HZ'],
        'TYR': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'HD1', 'HD2', 'HE1', 'HE2', 'HH'],
        'TRP': ['CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2', 'HE3', 'HZ2', 'HZ3', 'HH2', 'HE1'],
        'HIS': ['CG', 'ND1', 'CD2', 'CE1', 'NE2', 'HD1', 'HD2', 'HE1', 'HE2', 'xx', 'yy']  # if needed un comment
    }
    ring_atoms = {
        'PHE': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
        'TYR': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
        'TRP': ['CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
        'HIS': ['CG', 'ND1', 'CD2', 'CE1', 'NE2', ]  # if needed un comment
    }
    aromtic_residues = sorted(
        list(set([(int(i[0]), i[1], i[2]) for i in pdb_data[1].keys() if i[2] in aromatic_atoms.keys()])))
    ar_info = {}
    for m in pdb_data.keys():
        ar_info[m] = {}
        for ar_res in aromtic_residues:
            p = []
            for atm in ring_atoms[ar_res[2]]:
                p.append(pdb_data[m][(str(ar_res[0]), ar_res[1], ar_res[2], atm)])
            c = get_centroid(p)
            ar_info[m][(str(ar_res[0]), ar_res[1], ar_res[2])] = (c, p)
    return ar_info


def calculate_interaction(pdb,bmrb):
    pdb_data_actual = get_data.get_pdb_data(pdb)
    bmrb_data_actual = get_data.get_bmrb_data(bmrb)
    pdb_data_auth = get_data.get_pdb_data(pdb)
    bmrb_data_auth = get_data.get_bmrb_data(bmrb)
    pdb_actual_keys = [i for i in pdb_data_actual[1].keys() if i[3] == 'H']
    pdb_auth_keys = [i for i in pdb_data_auth[1].keys() if i[3] == 'H']
    bmrb_actual_keys = [i for i in bmrb_data_actual[0].keys() if i[3]=='H']
    bmrb_auth_keys = [i for i in bmrb_data_auth[0].keys() if i[3] == 'H']
    actual_actual = [i for i in bmrb_actual_keys if i in pdb_actual_keys]
    actual_auth = [i for i in bmrb_actual_keys if i in pdb_auth_keys]
    auth_actual = [i for i in bmrb_auth_keys if i in pdb_actual_keys]
    auth_auth= [i for i in bmrb_auth_keys if i in pdb_auth_keys]
    if len(bmrb_actual_keys)>0:
        seq_len = float(len(bmrb_actual_keys))
        actual_actual_match = float(len(actual_actual))/seq_len
        actual_auth_match = float(len(actual_auth))/seq_len
        auth_actual_match = float(len(auth_actual))/seq_len
        auth_auth_match = float(len(auth_auth))/seq_len
        if actual_actual_match > 0.7:
            pdb_data = pdb_data_actual
            bmrb_data = bmrb_data_actual
            tag_match = 'ORIG_ORIG'
        elif actual_auth_match > 0.7:
            pdb_data = pdb_data_auth
            bmrb_data = bmrb_data_actual
            tag_match = 'ORIG_AUTH'
        elif auth_actual_match > 0.7:
            pdb_data = pdb_data_actual
            bmrb_data = bmrb_data_auth
            tag_match = 'AUTH_ORIG'
        elif auth_auth_match > 0.7:
            pdb_data = pdb_data_auth
            bmrb_data = bmrb_data_auth
            tag_mathc = 'AUTH_AUTH'
        else:
            raise KeyError
        amide_chemical_shift, aromatic_chemical_shift, entity_size, assembly_size = bmrb_data
        ar_info = get_aromatic_info(pdb_data)
        dist,angle,solid_angle = analyze_enzemble(ar_info,pdb_data)
        if not os.path.isdir('./data'):
            os.system('mkdir ./data')
        if not os.path.isdir('./data/output'):
            os.system('mkdir ./data/output')
        fo_name = 'data/output/{}-{}-{}.csv'.format(pdb,bmrb,tag_match)
        fo=open(fo_name,'w')
        header = '#pdb_id,bmrb_id,entity_size,assembly_size,' \
                 'amide_seq_id,amide_chain_id,amide_res,amide_cs,amide_z,' \
                 'aro_seq_id,aro_chain_id,aro_res,' \
                 'mean_d,meidan_d,std_d,min_d,max_d,' \
                 'mena_ang,median_ang,std_and,min_ang,max_ang,' \
                 'mean_sang,median_sang,std_sang,min_sang,max_sang\n'
        fo.write(header)
        for atm in amide_chemical_shift.keys():
            # print (atm,amide_chemical_shift[atm],dist[atm][0],angle[atm][dist[atm][0][0]],
            #        solid_angle[atm][dist[atm][0][0]] )
            fo.write('{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n'.format(pdb,bmrb,
                                                                             entity_size,assembly_size,
                                                                             atm[0],atm[1],atm[2],
                   amide_chemical_shift[atm][0],amide_chemical_shift[atm][1],
                   dist[atm][0][0][0],dist[atm][0][0][1],dist[atm][0][0][2],
                   ','.join([str(i) for i in dist[atm][0][1]]),
                   ','.join([str(i) for i in angle[atm][dist[atm][0][0]]]),
                   ','.join([str(i) for i in solid_angle[atm][dist[atm][0][0]]])))
        fo.close()

def analyze_enzemble(ar_info,pdb_data):
    amide_res = {}
    for m in pdb_data.keys():
        amide_res[m]={}
        for atm in pdb_data[m].keys():
            if atm[3]=='H':
                amide_res[m][atm]=pdb_data[m][atm]
    dist={}
    angle={}
    solid_angle={}
    for atm in amide_res[1].keys():
        dist[atm]={}
        angle[atm]={}
        solid_angle[atm]={}
        for ar_res in ar_info[1].keys():
            d=[]
            a=[]
            sa=[]
            for m in amide_res.keys():
                d1=get_distance(ar_info[m][ar_res][0],amide_res[m][atm])
                d.append(d1)
                ang = find_angle(ar_info[m][ar_res],amide_res[m][atm],d1)
                a.append(ang[0])
                sa.append(ang[1])
            dist[atm][ar_res]=[round(numpy.mean(d),3),round(numpy.median(d),3),round(numpy.std(d),3),round(min(d),3),round(max(d),3)]
            angle[atm][ar_res]=[round(numpy.mean(a),3),round(numpy.median(a),3),round(numpy.std(a),3),round(min(a),3),round(max(a),3)]
            solid_angle[atm][ar_res]=[round(numpy.mean(sa),3),round(numpy.median(sa),3),round(numpy.std(sa),3),round(min(sa),3),round(max(sa),3)]
    s_dist={}
    for aa in dist.keys():
        sorted_dist = sorted(dist[aa].items(), key=lambda kv: kv[1][0])
        s_dist[aa]=sorted_dist
    return s_dist,angle,solid_angle

if __name__=="__main__":
    plot_3d('data/output/1WYO-11086-ORIG_ORIG.csv')
    #calculate_interaction('1WYO','11086')