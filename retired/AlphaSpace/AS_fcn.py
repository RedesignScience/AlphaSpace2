#AlphaSpace v1.0 - python functions

from retired.AlphaSpace import *
from retired.AlphaSpace import pdbinfo


def runAlphaSpace(pdb_lines, param_file = "AS_param.txt",option_file = "AS_option.txt", default = False, isWrite = True):
    from retired.AlphaSpace import readPDB
    # Default parameters
    param_default = {'core_cutoff': 100.0, 'min_num_alph': 1.0, 'hit_dist': 1.6, \
    'max_r': 5.4, 'min_r': 3.2, 'clust_dist': 4.7, 'aux_cutoff': 30.0, 'tight_communites_cutoff': 8.5, \
    'beta_clust_cutoff': 1.5, 'beta_class_cutoff': 4.2, 'beta_high_cutoff': 20, \
    'beta_mid_cutoff': 10, 'lig_clust_dist': 4.0} 

    option_default = {'screen_by_score': False, 'min_score': 100.0, 'use_naccess': True,\
    'get_face': True, 'expand_around_cntct': False, \
    'screen_by_perc_rank': False, 'max_desolv_perc': 0.95,\
    'do_reverse': False, 'screen_out_subsurf': False, \
    'screen_by_face': True, 'use_exact_ACSA': False,\
    'output_dir': '.', 'screen_face_perc': 0.5, \
    'output_to_screen': False, 'min_perc_rank': 0.9, \
    'use_TER': True, 'screen_by_lig_cntct': False,\
    'use_NP_wt': True, 'screen_by_res': False, \
    'pocket_communities': True, 'tight_option': True, \
    'contact_score': True, 'beta_cluster': True, 'get_beta_vol': False,\
    'screen_by_seam': False, 'exapand_around_seam': False}



    #get the location of the AlphaSpace.py program
    program_dir = os.path.dirname(os.path.realpath(__file__))

    #read in parameters (and options) from file into a dict and put the keys into an ordered list (for ordered output)
    
    if default:
        param_dict = param_default
        param_list = list(param_dict.keys())
        option_dict = option_default
        option_list = list(option_dict.keys())

    else:
        param_dict, param_list = get_params(param_file, program_dir)
        option_dict, option_list = get_params(option_file, program_dir)


    #create output directory if does not exist
    if not os.path.isdir(option_dict["output_dir"]):
        os.system("mkdir " + option_dict["output_dir"])


    #needed in working directory for Naccess calculations with alpha-atoms
    os.system("cp " + program_dir + "/vdw.radii .")
    
    #copy option file to working directory
    if not os.path.isfile('AS_option.txt'):
        os.system("cp " + program_dir + "/AS_option.txt .")

    #copy option file to working directory
    if not os.path.isfile('AS_param.txt'):
        os.system("cp " + program_dir + "/AS_param.txt .")

    #delete any existing log file
    if os.path.isfile(option_dict["output_dir"] + '/log.txt'):
        os.system("rm " + option_dict["output_dir"] + "/log.txt")


    #log program and version
    log("\n\nAlphaSpace v1.0\n\n",option_dict,to_screen=True)

    #output parameter list to log file
    log("--Parameter settings--\n",option_dict,to_screen=option_dict["output_to_screen"])
    for key in param_list:
        log(key + " = " + str(param_dict[key]), option_dict, to_screen=option_dict["output_to_screen"])
    log("\n",option_dict,to_screen=option_dict["output_to_screen"])

    #output option list to log file and screen
    log("--Option settings--\n",option_dict,to_screen=option_dict["output_to_screen"])
    for key in option_list:
        log(key + " = " + str(option_dict[key]), option_dict, to_screen=option_dict["output_to_screen"])
    log("\n\n",option_dict,to_screen=option_dict["output_to_screen"])

    if option_dict['use_naccess']:
        if not is_tool("naccess"):
            print("Notice: Naccess cannot be found!\n")
            print("Try:")
            print("1: follow the tutorial to make Naccess executable from command line by typing 'naccess'")
            print("2: run AlphaSpace without the Naccess program by setting the option: 'use_naccess = False'\n")
            sys.exit()

    #load the PDB file
    pdb = readPDB(pdb_lines,isTER = option_dict['use_TER'],isReverse = option_dict['do_reverse'], \
      useNaccess = option_dict['use_naccess'],getFace = option_dict['get_face'])

    
    #perform the Voronoi tessellation
    cutvert,cutsimp,cutrad = genVoro(pdb, min_r = param_dict['min_r'], max_r = param_dict['max_r'])

    #check for screen by resid
    #if screening by residue, grab the resid list from the input file
    if option_dict['screen_by_res']:
        res_list_all = [pdbinfo.resi(line) for line in open(option_dict['res_file'], 'r') if pdbinfo.isAtom(line)]
        resid_list = list(set(res_list_all))
    else:
        resid_list = []

    #log the start of mapping
    log("\nmapping the interface...\n",option_dict,to_screen=True,to_file=False)

    #perform fragmnet-centric topographical mapping (FCTM)
    ss = genPocket_wVoro(pdb, cutvert, cutsimp, cutrad, is_use_naccess = option_dict['use_naccess'], \
        is_get_face = option_dict['get_face'], clust_dist = param_dict['clust_dist'], \
        min_num_sph = param_dict['min_num_alph'], screen_by_face_atoms = option_dict['screen_by_face'], \
        screen_face_perc = option_dict['screen_face_perc'], hit_dist = param_dict['hit_dist'], \
        screen_by_score=option_dict['screen_by_score'],min_score=option_dict['min_score'], \
        screen_out_subsurf = option_dict['screen_out_subsurf'], subsurf_perc = option_dict['max_desolv_perc'], \
        is_use_exact_ACSA = option_dict['use_exact_ACSA'],screen_by_lig_contact = option_dict['screen_by_lig_cntct'], \
        expand_around_cntct = option_dict['expand_around_cntct'],is_use_NP_wt = option_dict['use_NP_wt'],clust_method = 'average',\
        screen_by_resid = option_dict['screen_by_res'], resid_list = resid_list, screen_by_seam = option_dict['screen_by_seam'],\
        expand_around_seam = option_dict['expand_around_seam'])

    #detect communities and calculate features
    if option_dict['pocket_communities']:
        pocket_communities = genCommunityPocket(pdb, ss, tight_option=option_dict["tight_option"], \
            tight_cutoff=param_dict["tight_communities_cutoff"])

    if option_dict['contact_score']:
        lig_resids = get_resids(option_dict['lig_resid_file'],ss.lig_pdb)
        cntct_score_arr = ss.get_cntct(pdb_lines,lig_resids,option_dict,param_dict)

        #4/4/16: adding this to get "new" cntct output with total space by residue etc
        cntct_dict = ss.get_cntct_NEW(option_dict,param_dict)

        e_cntct_dict = ss.get_e_cntct_NEW(pdb_lines,option_dict,param_dict)

    if isWrite:
        print("\nWriting...\n")

        #copy some files
        #needed in working directory for "draw" script
        if not os.path.isfile('colors_chimera.txt'):
            os.system("cp " + program_dir + "/colors_chimera.txt .")
        #needed in working directory for writing score table
        if not os.path.isfile('colors_table.txt'):
            os.system("cp " + program_dir + "/colors_table.txt .")
        #needed in output directory for chimera visualization
        if not os.path.isfile(option_dict["output_dir"] + '/colors_chimera.txt'):
            os.system("cp " + program_dir + "/colors_chimera.txt " + option_dict["output_dir"])
        #needed in output directory for chimera visualization
        if not os.path.isfile(option_dict["output_dir"] + '/AS_chimera.py'):
            os.system("cp " + program_dir + "/AS_chimera.py .")
        if not os.path.isfile(option_dict["output_dir"] + '/AS_chimera_beta.py'):
            os.system("cp " + program_dir + "/AS_chimera_beta.py .")
        if not os.path.isfile(option_dict["output_dir"] + '/AS_chimera_resid.py'):
            os.system("cp " + program_dir + "/AS_chimera_resid.py .")        


        #write out protein, binder, and interface PDB files
        ss.write_pdb(option_dict['output_dir'],param_dict['lig_clust_dist'])

        #write out pocket PDB files for map visualization
        ss.write_pockets(param_dict['core_cutoff'],param_dict['aux_cutoff'],option_dict['output_dir'])

        #generate file used to visualize alpha-spaces
        ss.write_tetras(option_dict['output_dir'])

        #write out score table (includes: pocket rank, color, score, space, percent occupied, and percent nonpolar)
        ss.write_score_table(option_dict['output_dir'],write_header=True)

        #remove some unneeded files from working dir and output dir    
        #os.system("rm vdw.radii colors_table.txt")
        os.system("rm vdw.radii")

        #write out pocket communities
        if option_dict['pocket_communities']:
            writeCommunityPocket(pocket_communities)

        #write out contact scores by binder resids
        if option_dict["contact_score"]:
            ss.write_cntct(cntct_score_arr,dir_name=option_dict['output_dir'])

            #4/4/16: adding this to get "new" cntct output with total space by residue etc
            ss.write_cntct_NEW(cntct_dict,dir_name=option_dict['output_dir'])
            ss.write_cntct_NEW_2(cntct_dict,e_cntct_dict,dir_name=option_dict['output_dir'])

            ss.write_resid_pockets(cntct_dict,dir_name=option_dict['output_dir'])
            ss.write_resid_pockets_2(cntct_dict,dir_name=option_dict['output_dir'])


        if option_dict["beta_cluster"]:
            ss.fill_beta_frag(param_dict['beta_clust_cutoff'])
            ss.write_beta_frag(param_dict['beta_class_cutoff'],option_dict['output_dir'],\
                param_dict['hit_dist'], param_dict['core_cutoff'],param_dict['aux_cutoff'],\
                param_dict['beta_high_cutoff'],param_dict['beta_mid_cutoff'],\
                option_dict['get_beta_vol'])

    return pdb, ss




def combine_prot_lig(pdb_file_list):
    """
    combine all pdb files from the list into one, with TER between first and all the others 
    """

    new_pdb = [line for line in open(pdb_file_list[0]) if not line.startswith("TER")]

    new_pdb.append("TER")

    for pdb in pdb_file_list[1:]:
        tmp_pdb = [line for line in open(pdb)]
        new_pdb = new_pdb + tmp_pdb

    return new_pdb




def get_resids(lig_resid_file,lig_pdb):
    """
    return a list of the resids to use in calculating contact scores
    """
    if not os.path.isfile(lig_resid_file):
        resid_list_all = [line[17:26].strip() for line in lig_pdb]
        print('\nNo specific "resID" file found > will calculate CNTCT SCORE for all binder residues...\n')
        #print("Can't find the binder resid file!! So not calculating any contact scores!\n\n")
    else:
        resid_list_all = [line[17:26].strip() for line in open(lig_resid_file)]
        print('\n"resID" list file found for CNTCT SCORE output...\n')
    lig_resids = []
    for r in resid_list_all:
        if r not in lig_resids:
            lig_resids.append(r)
    # #4-16-16: merging the resid and resname early now for the new contact score
    # lig_resids_merge = []
    # for r in lig_resids:
    #     merge_resid = r.split()[0] + r.split()[-1]
    #     lig_resids_merge.append(merge_resid)

    return lig_resids


# def get_resids(lig_resid_file,lig_pdb):
#     """
#     return a list of the resids to use in calculating contact scores
#     """
#     if not os.path.isfile(lig_resid_file):
#         resid_list_all = [line[17:26].strip() for line in lig_pdb]
#         print('\nNo specific "resID" file found > will calculate CNTCT SCORE for all binder residues...\n')
#         #print("Can't find the binder resid file!! So not calculating any contact scores!\n\n")
#     else:
#         resid_list_all = [line[17:26].strip() for line in open(lig_resid_file)]
#         print('\n"resID" list file found for CNTCT SCORE output...\n')
#     lig_resids = []
#     for r in resid_list_all:
#         if r not in lig_resids:
#             lig_resids.append(r)
#     #4-16-16: merging the resid and resname early now for the new contact score
#     lig_resids_merge = []
#     for r in lig_resids:
#         merge_resid = r.split()[0] + r.split()[-1]
#         lig_resids_merge.append(merge_resid)

#     return lig_resids_merge


def log(astr, option_dict, to_screen=True, to_file=True):
    '''
    Function: Output AlphaSpace statements, either to the screen or to a file
    Arguments: the string to output, the options dictionary, bool for screen option
    '''
    #print it to screen
    if to_screen:
        print (astr)
    # write it to the output file
    if to_file:
        try:
            f = open(option_dict["output_dir"] + '/log.txt', 'a')
            f.write(astr + "\n")
            f.close()
        except: pass


def get_params(param_file,program_dir):
    """
    Function: read in parameters (or options), evaluate/assign variable type, 
              return a parameter dictionary and an ordered list of the dict's keys
    Arguments: name of parameter file, path to AlphaSpace program directory
    """
    param_dict = {}
    
    #use paramter and option file from current working directory if exist there 
    #otherwise, use the parameter and option file in the program directory
    if os.path.isfile(param_file):
        param_path = param_file
    elif os.path.isfile(program_dir + '/' + param_file):
        param_path = program_dir + '/' + param_file 
    else:
        sys.exit('ERROR - exiting AlphaSpace: Cannot find required input file: "' + param_file + '"')

    # if os.path.isfile(option_file):
    #     param_path = option_file
    # elif os.path.isfile(program_dir + '/' + option_file):
    #     param_path = program_dir + '/' + option_file 
    # else:
    #     sys.exit('ERROR - exiting AlphaSpace: Cannot find option file "' + option_file + '"')

    param_list = []
    for line in (lines.rstrip().replace(' ','').split("=") for lines in open(param_path) \
      if (not lines.startswith('#') and '=' in lines)):

        #check if the parameter is a boolian value
        if any(x in line[1] for x in ["True","true"]):
            param_dict[line[0]] = True
            param_list.append(line[0])
        elif any(x in line[1] for x in ["False","false"]):
            param_dict[line[0]] = False
            param_list.append(line[0])
        #check if the parameter is numerical (first character of parameter is a digit)
        elif line[1][0].isdigit():
            param_dict[line[0]] = float(line[1])
            param_list.append(line[0])
        #otherwise just keep the parameter as a string
        else:
            param_dict[line[0]] = line[1]
            param_list.append(line[0])
    return param_dict, param_list


def genVoro(pdb, min_r, max_r):
    """
    pass a pdb object, get face_idx, run voronoi, and return values to pass to genPocket
    this is to be used for scanning clustering parameters, like distance to get modularity
    """

    coord = [pdbinfo.coord(l) for l in pdb.prot]
    coord = np.array(coord)


    from scipy.spatial import Voronoi, Delaunay
    # use qhull to do the calculation
    vor = Voronoi(coord)
    tess = Delaunay(coord)
   
    #calculate radius for each alpha sphere
    rad = calcRadius(coord, tess.simplices, vor.vertices)

    #filter out spheres below a radius min or above a radius max
    ind = cutIndex(rad,min_r,max_r)

    #filter out spheres without all 4 atoms on the surface (ie: asa > 0)
    #surf_idx = pdb.get_surf_idx()
    tess_simplices = tess.simplices
    #ind = filtIndex(ind,surf_idx,tess_simplices)

    #get voronoi information for final filtered vertices_coords
    cutvert = get_fr_arr(vor.vertices, ind)
    cutsimp = get_fr_arr(tess.simplices, ind)
    cutrad = get_fr_arr(rad, ind)

    return [cutvert,cutsimp,cutrad]



def genPocket_wVoro(pdb, cutvert, cutsimp, cutrad, is_use_naccess = False, is_get_face = False, \
              clust_dist = 4.7, min_num_sph = 1, screen_by_face_atoms = True,\
              screen_face_perc = 0.5, screen_by_lig_contact = False, hit_dist = 1.6, \
              expand_around_cntct = False, is_use_NP_wt = False,\
              is_use_exact_ACSA = False, \
              screen_by_score = False, min_score = 100, \
              screen_by_perc_rank = False, min_perc_rank = 0.9,\
              screen_out_subsurf = False, subsurf_perc = 0.95, clust_method = 'average',\
              screen_by_resid = False, resid_list = [], screen_by_seam = False, expand_around_seam = False):

    """
    Generate Pocket by
    (1) Voronoi tessellation
    (2) Refine by cutoff radius of alpha sphere, surfatom by NACCESS (ASA > 0)
    (3) Subpocket: linkage clustering
    (4) Superpocket: clustering the centoid of subpocekt

    Return: sub_pocket list, super_pocket_list

    this now to generate pockets with preestablished voronoi definitions

    """
    from retired.AlphaSpace import Pocket,Snapshot
    #clust_method = 'average'
    #clust_method = 'complete'

    #only do the clustering if there are vertices_coords to cluster... otherwise return empty Snapshot
    if len(cutvert) > 1:

        #do the clustering  
        zmat = hier.linkage(cutvert, method=clust_method)
        cluster = hier.fcluster(zmat,clust_dist,criterion='distance')
    
    elif len(cutvert) == 1:
        cluster = [1]

    else: 
        return Snapshot()


    p_verts = [ [] for i in range(max(cluster)) ]
    p_atms = [ [] for i in range(max(cluster)) ]
    
    fill_p_verts(p_verts,cluster)
    fill_p_atms(p_atms,p_verts,cutsimp)

    sub_pockets = []
    for i,p in enumerate(p_verts):
        sub_pockets.append(Pocket(p,p_atms[i]))

    #filter pockets by size
    if min_num_sph > 1:
        sub_pockets = filter_by_size(sub_pockets, min_num_sph)

    #screen pockets based on the percentage of each pocket-atom list 
    #that is within the list of interface atoms
    #Edit 11/16/15: took out requirement of "is_use_naccess" because I've added face calculation
    if screen_by_face_atoms and is_get_face:
        sub_pockets = screen_p_by_face(sub_pockets,pdb.face_idx,screen_face_perc)
    elif screen_by_face_atoms:
        print ("Notice: option 'get_face = False' so NOT assigning interface atoms...\n\
        To properly screen pockets using interface atoms, set option: 'get_face = True' ")

    ss = Snapshot(pdb,cutvert,cutsimp,cutrad,sub_pockets)


    if screen_by_lig_contact:
        cntct_p = []
        for p in ss.sub_pockets:
            if p.is_contact(hit_dist):
                cntct_p.append(p)
        #option to expand lig contact pocket list to include pockets with any atoms shared with 
        #lig contact pocket
        if expand_around_cntct:
            cntct_p_exp = []
            for p in cntct_p:
                if p not in cntct_p_exp:
                    cntct_p_exp.append(p)
                for p2 in ss.sub_pockets:
                    if p2 not in cntct_p_exp:
                        if any(i in p.atoms for i in p2.atoms):
                            cntct_p_exp.append(p2)
            ss.sub_pockets = cntct_p_exp
        else:
            ss.sub_pockets = cntct_p            

    if screen_by_resid:
        new_pockets = []
        for p in ss.sub_pockets:
            found = 0
            for a_1 in p.atoms:
                for resid in resid_list:
                    if resid == pdbinfo.resi(ss.prot_pdb[a_1]):
                        found = 1
                        new_pockets.append(p)
                        break
                if found:
                    break
        ss.sub_pockets = new_pockets

    #screen for pockets that contain atoms from more than 1 protein chain
    #note: probably want to use screen_by_subsurf with this to avoid getting buried interfacial pockets
    if screen_by_seam:
        seam_pockets = []
        for p in ss.sub_pockets:
            chain_IDs = []
            for a in p.atoms:
                curr_chain = pdbinfo.chid(ss.prot_pdb[a])
                if curr_chain not in chain_IDs:
                    chain_IDs.append(curr_chain)
            if len(chain_IDs) > 1:
                seam_pockets.append(p)
        if expand_around_seam:
            seam_p_expand = []
            for p in seam_pockets:
                if p not in seam_p_expand:
                    seam_p_expand.append(p)
                for p2 in ss.sub_pockets:
                    if p2 not in seam_p_expand:
                        if any(i in p.atoms for i in p2.atoms):
                            seam_p_expand.append(p2)
            ss.sub_pockets = seam_p_expand 
        else:
            ss.sub_pockets = seam_pockets


    #assign score features to pockets

    #if we are using non-polar weighting in the pocket score
    if is_use_NP_wt and is_use_naccess:
        #assign asa by atom for the pockets, either pocket by pocket or, 
        #if "not exact" calculate once with all pockets
        if not is_use_exact_ACSA:
            ss.set_p_asa_prot()
        for p in ss.sub_pockets:
            if is_use_exact_ACSA:
                p.asa_prot = p.get_exact_asa_prot()
            else:
                p.asa_prot = p.ss.p_asa_prot
            p.set_space_score_occ(hit_dist)

        #sort pockets by score
        list.sort(ss.sub_pockets, key=lambda x: x.score, reverse=True)

    #otherwise use just alpha-space as the pocket score
    else:
        for p in ss.sub_pockets:
            
            #assign 1.0 for atomic ACSAs for all protein atoms; can still perform Match but it's binary in/out of pocket
            p.asa_prot = [1.0] * len(p.ss.prot_pdb)
            
            #this option to bypass non-polar weighting altogether
            #p.set_space_as_score_occ(hit_dist)
            
            #this option to use non-polar weighting even though not calculating atomistic ASAs (all assigned as 1)
            p.set_space_score_occ(hit_dist)

        list.sort(ss.sub_pockets, key=lambda x: x.score, reverse=True)

    #if requested, screen pockets by a minimum score
    if screen_by_score:
        screen_pockets = []
        for p in ss.sub_pockets:
            if p.score >= min_score:
                screen_pockets.append(p)
        ss.sub_pockets = screen_pockets

    #if requested, screen remaining pockets by there percentage rank 
    if screen_by_perc_rank:
        screen_pockets = []
        for i,p in enumerate(ss.sub_pockets):
            if (1.0 - ((i + 1.0)/float((len(ss.sub_pockets))))) < min_perc_rank:
                break
            else:
                screen_pockets.append(p)
        ss.sub_pockets = screen_pockets

    #if requested, screen out pockets if not enough exposure to the solvent (i.e. beneath the surface)
    if screen_out_subsurf and is_use_naccess:
        screen_pockets = []
        perc_sa_cov = ss.get_clusters_percSA()
        for i,p in enumerate(ss.sub_pockets):
            if perc_sa_cov[i] < subsurf_perc:
                screen_pockets.append(p)
        ss.sub_pockets = screen_pockets

    return ss



def fill_p_verts(p_verts,cluster):
    """
    Function: fill an empty list of lists (length equal to the number of clusters) with vertice indices according to the clustering (becomes centroids in super clustering)
    Parameters
    p_verts: empty list of lists
    cluster: cluster indices 
    """
    
    for i,clust_idx in enumerate(cluster):
        p_idx = clust_idx - 1
        p_verts[p_idx].append(i)
    return


def fill_p_atms(p_atms,p_verts,simp):
    """
    Function: fill an empty list of lists (length equal to the number of clusters) with non-redundant atom lists according to the clustering
    Parameters
    p_atms: empty list of lists
    p_verts: the prefilled pocket vertice list of lists
    simp: the filtered list of simplex_idx_list
    """
    
    for i,p_vert in enumerate(p_verts):
        for j in p_vert:
            for k in simp[j]:
                if k not in p_atms[i]:
                    p_atms[i].append(k)
        p_atms[i].sort()
    return


def filter_by_size(pockets,min_sph):
    """
    Funtion: only keep pockets defined by a minimum number of alpha spheres
    """
    filt_p = []
    for p in pockets:
        if len(p.vertices) >= min_sph:
            filt_p.append(p)
    return filt_p


def calc_asa(a_list):
    """
    Function
    write the atom list to a file and call Naccess to calculate the ASA

    Parameters
    a_list: atom list in pdb format 
    
    Return
    list of ASA for each atom
    """
    
    name = "tmpp"
    f = open(name + ".pdb",'w')
    for l in a_list:
        #change HETATM to ATOM for protein-binder complexes
        if pdbinfo.isHETATM(l):
            f.write("ATOM  " + l[6:] + '\n')
        else:
            f.write("%s\n" % l) 
    f.close()
    
    os.system("naccess " + name + ".pdb")
    #os.system("naccess -p 1.0 " + name + ".pdb")
    output = name + ".asa"
    with open(output,'r') as f:
        pdb_asa = f.read().splitlines() 
    asa_list = [float(l[54:62]) for l in pdb_asa] 
    
    suffix = ['.pdb','.rsa','.asa','.log']
    for s in suffix:
        os.system("rm " + name + s)
    
    return asa_list

def calc_asa2(a_list):
    """
    Function
    write the atom list to a file and call Naccess to calculate the ASA
    *** comparing with calc_asa, this will also return radius list ***
    *** from naccess method ***

    Parameters
    a_list: atom list in pdb format 
    
    Return
    list of ASA for each atom
    """    
    name = "tmpp"
    f = open(name + ".pdb",'w')
    for l in a_list:
        #change HETATM to ATOM for protein-binder complexes
        if pdbinfo.isHETATM(l):
            f.write("ATOM  " + l[6:] + '\n')
        else:
            f.write("%s\n" % l) 
    f.close()
    
    os.system("naccess " + name + ".pdb")
    #os.system("naccess -p 1.0 " + name + ".pdb")
    output = name + ".asa"
    with open(output,'r') as f:
        pdb_asa = f.read().splitlines() 
    asa_list = [float(l[54:62]) for l in pdb_asa] 
    asa_radius = [float(l[64:69]) for l in pdb_asa]
    suffix = ['.pdb','.rsa','.asa','.log']
    for s in suffix:
        os.system("rm " + name + s)
    
    return asa_list, asa_radius


def calc_simp_vol(pts):
    """
    calculate the volume of the simplex defined by 4 points in an array
    """
    v = [pts[i] - pts[0] for i in range(1,4,1)]
    
    return np.round(abs(np.dot(v[0], np.cross(v[1], v[2]))) / 6.0, 4)




def calcRadius(coord, tess_simplices, vor_vertices):
    """
    Parameters
    Coord:  input atoms of coord
    tess_simplices : four points form alpha sphere
    vor_vertices: the coord of center of four points form alpha sphere
    
    Return
    Radius of each alpha sphere in the order of input tess_simplices
    """
    
    radius = []
    lenalpha = len(tess_simplices)
    for i in range(lenalpha):
        radius.append(calcDist(coord[tess_simplices[i][0]], vor_vertices[i]))
    return np.array(radius)



def calcDist(v1, v2):
    """
    Function to calculation the distance between two points
    """
    
    return sum((v1 - v2)**2)**0.5


def cutIndex(radius,min_r,max_r):
    """
    Function to get the index of radius range from min_r to max_r angstrom
    """
    
    ind = []
    rlen = len(radius)
    for i in range(rlen):
        if radius[i] >= min_r and radius[i] <= max_r:
           #if radius[i] <= 6.0:
            ind.append(i)
    return ind


def get_fr_arr(arr, ind):
    """
    get the selected array items from the index list
    """
    cut_arr = []
    arr = list(arr)
    for i in ind:
        cut_arr.append(arr[i])
    cut_arr = np.array(cut_arr) 
    return cut_arr


def screen_p_by_face(pockets,face,percent):
    """
    Function: screen the pockets by percentage of lining atoms within the interface. (note: "face" is an atom index list)
    """
    tmp_p = []
    for p in pockets:
        face_hits = 0
        tot = 0
        for a_idx in p.atoms:
            tot += 1
            if a_idx in face:
                face_hits += 1
        curr_perc = float(face_hits) / float(tot)
        if curr_perc >= percent:
            tmp_p.append(p)
    return tmp_p


def estVolume(tess_c, tess_r, niter=1000):
    """
    Calculate the volume by direct monte carlo method
    """
    ntess = len(tess_r)
    # get the position of the box for volume estimation
    # by get the max[x,y,z] and min[x,y,z]
    # calculate the volume of box
    xyzmax = []
    xyzmin = []
    for i in range(3):
        xyzmax.append(max(tess_c[:,i] + tess_r))
        xyzmin.append(min(tess_c[:,i] - tess_r))
    vbox = (xyzmax[0] - xyzmin[0])*(xyzmax[1] - xyzmin[1])*(xyzmax[2] - xyzmin[2])

    # initial hit to be zero
    # square the radius
    hit = 0.
    tess_r2 = tess_r ** 2
    
    for i in range(niter):
        randpt = np.array([random.uniform(xyzmin[i], xyzmax[i]) for i in range(3)])
        for j in range(ntess):
            dist2 = sum((tess_c[j] - randpt)**2)
            if dist2 <= tess_r2[j]:
                hit = hit +  1
                break
                
    return vbox * hit / niter



#for Pocket Communities

def writeCommunityPocket(sortedcompkt):
    """
    Write Pocket communities output table
   
    """
    f = open("table_community.dat", "w")
    #f.write("#Rank\t\t Score\t\t CoreID\t\t AuxID\t\t MinorID\n")
    f.write("#Rank\t\t Score\t\t OccScore\t\t UnoccScore\t\t CoreID\t\t AuxID\t\t MinorID\n")
    for i,p in enumerate(sortedcompkt):
        tmp, k = p
        #f.write(str(i+1) + "\t\t " + str(int(k[3])) + "\t\t ")
        f.write(str(i+1) + "\t\t " + str(int(k[3])) + "\t\t " + str(int(k[4])) + "\t\t " + str(int(k[5])) + "\t\t ")

        for m in range(3):
            if len(k[m]) > 0:
                tmp = [str(j+1) for j in k[m]]
                f.write(",".join(tmp) + "\t\t ")
            else:
                f.write("--\t\t ")
        f.write("\n")

    f.close()
    return None

def CoreCluster(pkt_core_list, ftdict, nclust, clust, CC_cut = 8.5):
    """Core Cluster by average linkage using distance cutoff
    
    :param pkt_core_list: the core pocket id list
    :param ftdict: dict of pocket information
    :param nclust: number of original cluster
    :param clust: cluste index
    :param CC_cut: the cutoff distance 
    """
    
    pktdict = dict(list(zip(pkt_core_list, clust)))
    for i in range(nclust):
        tmp_core_list = [key for key, value in pktdict.items() if value == i]
        if len(tmp_core_list) > 1:
            tmp_core_center = [ftdict[j][2] for j in tmp_core_list]
            tmp_zmat = hier.linkage(tmp_core_center, method="average")
            tmp_cluster = hier.fcluster(tmp_zmat, CC_cut, criterion='distance')
            ## tmp_cluster is 1 indexed 
            tmp_ncluster = max(tmp_cluster)
            if tmp_ncluster > 1:
                for m, n in zip(tmp_core_list, tmp_cluster):
                    if n > 1:
                        pktdict[m] = nclust + n - 2 
                nclust = nclust + tmp_ncluster - 1    
                          
            #print(tmp_core_center)
            #print(tmp_cluster)
            #print(tmp_core_list)
    #print(pktdict)
    newcluster = [pktdict[i] for i in pkt_core_list]
 
    return nclust, newcluster



def genCommunityPocket(pdb, ss, corecut = 100, auxcut = 30, tight_option = False, tight_cutoff = 8.5):
    """
    Generate Pocket Communities for one PDB   *this new version has "tight" option
    
    :param pdb: the pdb instance
    :param ss: pocket class
    :param corecut: cutoff score for core pocket
    :param auxcut: cutoff score for  auxiliary pocket
    """
    
    natom = len(pdb.prot)
    ssp = ss.sub_pockets
    pkt_all_arr = fill_pocket(natom, ssp)
    
    pkt_list = list(range(len(ssp)))
    pkt_core_list = []
    ftdict = {}
    
    for j,p in enumerate(ssp):
        ftdict[j] = [pkt_all_arr[j], p.atm_cntr(), \
                p.vert_cntr(), p.space, p.score, p.occ_cntct_score, p.unocc_cntct_score]
        if p.score >= corecut:
            pkt_core_list.append(j)
 
    compkt = {}      
    if len(pkt_core_list) > 1:
        coremat = []
        for i in pkt_core_list[:-1]:
            for j in pkt_core_list[(i+1):]:
                #print i,j
                #print(sum((ftdict[i][2] -ftdict[i][1]) * (ftdict[j][2] - ftdict[j][1])))
                
                #for pocket to "overlap" they must, first, share at least one atom, and, second,
                #if they are pointing away from each other (pocket atom centroids are closer to each other
                #than alpha-cluster centroids), then the angle between the directional pocket vectors must
                #be between 0 and 90 degrees
                if sum(pkt_all_arr[i] * pkt_all_arr[j]) > 0 and \
                    (sum((ftdict[i][2] -ftdict[i][1]) * (ftdict[j][2] - ftdict[j][1])) > 0 or \
                    (np.linalg.norm(ftdict[i][1] - ftdict[j][1]) > np.linalg.norm(ftdict[i][2] - ftdict[j][2]))):
                        coremat.append(1)
                else:
                    coremat.append(0) 
                    
        coremat = scidist.squareform(np.array(coremat))
        nclust, clust =  csg.connected_components(coremat, directed = False)
        
        
        if tight_option:
            # Do average linkage clustering
            nclust, clust = CoreCluster(pkt_core_list, ftdict, nclust, clust, CC_cut = tight_cutoff)
        
        
        for i in range(nclust):
            compkt[i] = [[],[],[], 0, 0, 0]

        for i,j in zip(clust, pkt_core_list):
            compkt[i][0].append(j)
            
    elif len(pkt_core_list) == 1:
        nclust = 1
        compkt[0] = [[pkt_core_list[0]],[],[], 0, 0, 0]
        
    elif len(pkt_core_list) == 0:
        nclust = 0
        print ("No Pocket with score > %d!" %(corecut))


    for i in range(nclust):
        for j in compkt[i][0]:
            for k in pkt_list:
                if k not in compkt[i][0]:
                    
                    #print i,j
                    #print(sum((ftdict[i][2] -ftdict[i][1]) * (ftdict[j][2] - ftdict[j][1])))
                    if tight_option:
                        if sum(pkt_all_arr[j] * pkt_all_arr[k]) > 0 and \
                        (sum((ftdict[j][2] -ftdict[j][1]) * (ftdict[k][2] - ftdict[k][1])) > 0 or \
                        (np.linalg.norm(ftdict[j][1] - ftdict[k][1]) > np.linalg.norm(ftdict[j][2] - ftdict[k][2]))) and \
                        np.linalg.norm(ftdict[j][2] - ftdict[k][2]) < tight_cutoff :
                            if ftdict[k][4] >= auxcut:
                                if k not in compkt[i][1]:
                                    compkt[i][1].append(k)
                            elif k not in compkt[i][2]:
                                compkt[i][2].append(k)

                    else:
                        if sum(pkt_all_arr[j] * pkt_all_arr[k]) > 0 and \
                        (sum((ftdict[j][2] -ftdict[j][1]) * (ftdict[k][2] - ftdict[k][1])) > 0 or \
                        (np.linalg.norm(ftdict[j][1] - ftdict[k][1]) > np.linalg.norm(ftdict[j][2] - ftdict[k][2]))):
                            if ftdict[k][4] >= auxcut:
                                if k not in compkt[i][1]:
                                    compkt[i][1].append(k)
                            elif k not in compkt[i][2]:
                                compkt[i][2].append(k)






                
    #print compkt
    for i in range(nclust):
        for j in compkt[i][0]:
            compkt[i][3] += ftdict[j][4]
            compkt[i][4] += ftdict[j][5]
            compkt[i][5] += ftdict[j][6]
        for j in compkt[i][1]:
            compkt[i][3] += ftdict[j][4]
            compkt[i][4] += ftdict[j][5]
            compkt[i][5] += ftdict[j][6]
        #don't count minor pocket scores into community score (because still want to preference more consolidated communities)
        #for j in compkt[i][2]:
        #    compkt[i][3] += ftdict[j][4]       
               
    sortedcompkt = sorted(list(compkt.items()), key = lambda x:x[1][3], reverse=True)
    #writeCommunityPocket(sortedcompkt)
    
    return sortedcompkt



###this to be removed if no longer needed
def genCommunityPocket_ORIG(pdb, ss, corecut = 100, auxcut = 30):
    """
    Generate Pocket Communities for one PDB
    
    :param pdb: the pdb instance
    :param ss: pocket class
    :param corecut: cutoff score for core pocket
    :param auxcut: cutoff score for  auxiliary pocket
    """
    
    natom = len(pdb.prot)
    ssp = ss.sub_pockets
    pkt_all_arr = fill_pocket(natom, ssp)
    
    pkt_list = list(range(len(ssp)))
    pkt_core_list = []
    ftdict = {}
    
    for j,p in enumerate(ssp):
        ftdict[j] = [pkt_all_arr[j], p.atm_cntr(), \
                p.vert_cntr(), p.space, p.score, p.occ_cntct_score, p.unocc_cntct_score]
        if p.score >= corecut:
            pkt_core_list.append(j)
 
    compkt = {}      
    if len(pkt_core_list) > 1:
        coremat = []
        for i in pkt_core_list[:-1]:
            for j in pkt_core_list[(i+1):]:
                #print i,j
                #print(sum((ftdict[i][2] -ftdict[i][1]) * (ftdict[j][2] - ftdict[j][1])))
                
                #for pocket to "overlap" they must, first, share at least one atom, and, second,
                #if they are pointing away from each other (pocket atom centroids are closer to each other
                #than alpha-cluster centroids), then the angle between the directional pocket vectors must
                #be between 0 and 90 degrees
                if sum(pkt_all_arr[i] * pkt_all_arr[j]) > 0 and \
                    (sum((ftdict[i][2] -ftdict[i][1]) * (ftdict[j][2] - ftdict[j][1])) > 0 or\
                    (np.linalg.norm(ftdict[i][1] - ftdict[j][1]) > np.linalg.norm(ftdict[i][2] - ftdict[j][2]))):
                        coremat.append(1)
                else:
                    coremat.append(0) 
                    
        coremat = scidist.squareform(np.array(coremat))
        nclust, clust =  csg.connected_components(coremat, directed = False)
        
        for i in range(nclust):
            compkt[i] = [[],[],[], 0, 0, 0]

        for i,j in zip(clust, pkt_core_list):
            compkt[i][0].append(j)
            
    elif len(pkt_core_list) == 1:
        nclust = 1
        compkt[0] = [[pkt_core_list[0]],[],[], 0, 0, 0]
        
    elif len(pkt_core_list) == 0:
        nclust = 0
        print ("No Pocket with score > %d!" %(corecut))
          
    for i in range(nclust):
        for j in compkt[i][0]:
            for k in pkt_list:
                if k not in compkt[i][0]:
                    
                    #print i,j
                    #print(sum((ftdict[i][2] -ftdict[i][1]) * (ftdict[j][2] - ftdict[j][1])))
                    
                    if sum(pkt_all_arr[j] * pkt_all_arr[k]) > 0 and \
                    (sum((ftdict[j][2] -ftdict[j][1]) * (ftdict[k][2] - ftdict[k][1])) > 0 or\
                    (np.linalg.norm(ftdict[j][1] - ftdict[k][1]) > np.linalg.norm(ftdict[j][2] - ftdict[k][2]))):
                        if ftdict[k][4] >= auxcut:
                            if k not in compkt[i][1]:
                                compkt[i][1].append(k)
                        elif k not in compkt[i][2]:
                            compkt[i][2].append(k)
                
    #print compkt
    for i in range(nclust):
        for j in compkt[i][0]:
            compkt[i][3] += ftdict[j][4]
            compkt[i][4] += ftdict[j][5]
            compkt[i][5] += ftdict[j][6]
        for j in compkt[i][1]:
            compkt[i][3] += ftdict[j][4]
            compkt[i][4] += ftdict[j][5]
            compkt[i][5] += ftdict[j][6]
        #don't count minor pocket scores into community score (because still want to preference more consolidated communities)
        #for j in compkt[i][2]:
        #    compkt[i][3] += ftdict[j][4]       
               
    sortedcompkt = sorted(list(compkt.items()), key = lambda x:x[1][3], reverse=True)
    #writeCommunityPocket(sortedcompkt)
    
    return sortedcompkt



##for Pocket matching

def checkPDB(pdblist, debug=True):
    """
    Function to check PDB for the same number of atoms
    
    :param pdblist: list of pdb file

    .. :: note
    """


    if os.path.isdir("Pocket_matching"):
        os.system("rm -r Pocket_matching")

    os.system("mkdir Pocket_matching")   

    newpdblist = []

    for i,fn in enumerate(pdblist):
        newd = "Pocket_matching/AS" + str(i+1) + "/"
        os.system("mkdir " + newd)
        os.system("cp " + fn  + " " + newd)
        
        outpdb = "orig" + str(i+1) + ".pdb"
        cmd =  "sed 's/HIE/HIS/g; s/HID/HIS/g; s/HIP/HIS/g' " + \
            fn + " > " + newd + outpdb
        os.system(cmd)
    
    cwd = os.getcwd()
    
    numsnap = len(pdblist)    
    
    pdblist = [ cwd + "/Pocket_matching/AS" + str(i+1) + "/orig" + str(i+1) + ".pdb" \
        for i in range(numsnap)]
    outpdblist = [ cwd + "/Pocket_matching/AS" + str(i+1) + "/mod" + str(i+1) + ".pdb" \
        for i in range(numsnap)]
    pdbdict = []
    
    ## Intersection of heavy atom list in all input PDBs
    ## Based on AtomName-ResNum-ResName
    for i,pdbnm in enumerate(pdblist):
        pdb_lines = [line for line in open(pdbnm)]
        pdb = readPDB(pdb_lines, isTER=True)
        if  0 == i:
            carryatomlist = []
            tmpdict = {}
            for lines in pdb.prot:
                idf = lines[12:16] + lines[17:21] + lines[22:26]
                carryatomlist.append(idf)
                tmpdict[idf] = lines
            carryatomlist = set(carryatomlist)
            pdbdict.append(tmpdict)
        else:
            tmpatomlist = []
            tmpdict = {}
            for lines in pdb.prot:
                idf = lines[12:16] + lines[17:21] + lines[22:26]
                tmpatomlist.append(idf)
                tmpdict[idf] = lines
            carryatomlist = carryatomlist.intersection(set(tmpatomlist))
            pdbdict.append(tmpdict)
            
    carryatomlist = list(carryatomlist)
   
    ## Output ordered heavy atom based on first PDB
    ## pdb named as tmpfix#.pdb with format
    ##Prot
    ##TER
    ##Lig
   
    for i,pdbnm in enumerate(pdblist):
        pdb_lines = [line for line in open(pdbnm)]
        pdb = readPDB(pdb_lines, isTER=True) 
        outpdb = outpdblist[i]
        f = open(outpdb, "w")
        if 0 == i:
            orderatomlist = []
            for lines in pdb.prot:
                idf = lines[12:16] + lines[17:21] + lines[22:26]
                if idf in carryatomlist:
                    orderatomlist.append(idf)
                    f.write(lines + "\n")
        else:
            for atm in orderatomlist:
                f.write(pdbdict[i][atm] + "\n")
       
        f.write("TER\n")
      
        for lines in pdb.lig:
            f.write(lines + "\n")
        f.close()
       
    ## Print striped atom in each PDB
    if debug:
        f = open("Pocket_matching/strippedAtoms.dat", "w")
        for i,pdbnm in enumerate(pdblist):
            f.write("\n" + pdbnm + "\n")
            pdb_lines = [line for line in open(pdbnm)]      
            pdb = readPDB(pdb_lines, isTER=True)
            for lines in pdb.prot:
                idf = lines[12:16] + lines[17:21] + lines[22:26]
                if idf not in orderatomlist:
                    f.write(idf + "\n")
                   
        f.close()
        
    return outpdblist       



def write_matchpocket(matchpocketidx, gsim, nsnap, sortidx, npocket_list):
    """Write match pocket
    
    
    """
    scdict = {}
    
    for i in range(nsnap):
        fn = "AS" + str(i+1) + "/table_pockets.dat"
        for lines in fileinput.input(fn):
            if lines[0:5] == "#rank":
                continue
            else:
                line = lines.split()
                scdict[(i,int(line[0])-1)] = line[2:4]
        fileinput.close()
    
    #print os.getcwd() 
    with open('AS1/colors_table.txt','r') as f:
        colors = f.read().splitlines()
    f.close()

    with open('AS1/colors_chimera.txt','r') as f:
        colorschi = f.read().splitlines()
    f.close()


    colorchidict  = {}
    for i,j in zip(colors, colorschi):
        colorchidict[i] = j

    colors = colors * 20

    
    f = open("Result_match.txt", "w")
    
    head = "#Rank    Color     GroupSim\t" + \
        "\t".join(["Score\tOcc\t"] * nsnap) + "\n"
    
    f.write(head)
     
    colordict = {}
    
    for i,j in enumerate(sortidx):
        pktlist = matchpocketidx[j]
        line = "%-8d%-10s%8.2f\t\t" %(i+1, colors[i], gsim[j])
        tmpdict = {}
        for k in pktlist:
            idx1, idx2 = k
            tmpdict[idx1] = k
            colordict[k] = colors[i]
        for n in range(nsnap):
            if n in list(tmpdict.keys()):
                line = line + "\t".join(scdict[tmpdict[n]]) + "\t\t"
            else:
                line = line + "--\t--\t\t"
        
        f.write(line + "\n")
    f.close()   
    
    for i,j in enumerate(npocket_list):
        f = open("AS" + str(i+1) + "/colors_match.txt", "w")
        for k in range(j):
            f.write(colorchidict[colordict[(i,k)]]  + "\n")
        f.close()
        cmd1 = "sed -e 's/colors_chimera/colors_match/g' AS" + str(i+1) + "/AS_chimera.py > AS" + str(i+1) + "/AS_chimera_match.py"
        os.system(cmd1)
    
    return None
    


def plot_heatmap(pocket_atom_array,  metric="euclidean",\
                 figname = "match_clust_heatmap.png"):
    """Plot Heatmap of clustering
    
    This function is useful not only for this project and it is helpful to put 
    into other modules
    
    :param pocket_atom_array: the coord or the pairwised matrix in a list 
    :param metric: metric to calculated the distance, "euclidean" or "jaccard"
    :param figname: figure name
    :returns: figure 
    
    .. todo:: Should be put into other modules or plot modules
    
    """
    
    import scipy.spatial.distance as scidist 
    import scipy.cluster.hierarchy as hier
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import matplotlib.font_manager as font_manager
    
    # Font of figure
    path = '/Library/Frameworks/EPD64.framework/Versions/Current/lib/python2.7/\
site-packages/matplotlib/mpl-data/fonts/ttf/Helvetica.ttf'
    try:
        prop = font_manager.FontProperties(fname=path)
        mpl.rcParams['font.family'] = prop.get_name()
    except:
        pass
    
    # color
    mpl.rc('axes',edgecolor='k')
    colors = ['#1e90ff', '#ff69b4', '#32cd32', '#9370db', '#2f4f4f','#87ceeb',\
          '#ff7f00','#228b22','#b22222','#000000']
    hier.set_link_color_palette(colors)
    
    # create pdist matrix (stored as list and see numpy for details)
    if metric == 'euclidean':
        data_dist = scidist.pdist(pocket_atom_array)
    elif metric == 'jaccard':
        data_dist = scidist.pdist(pocket_atom_array, 'jaccard')
    else:
        print ('Error')

    # figure size
    fig = plt.figure(figsize=(10,10))
    
    # x ywidth height
    ax1 = fig.add_axes([0.05,0.1,0.2,0.6])
    Y = hier.linkage(data_dist, method="average")
    Z1 = hier.dendrogram(Y, orientation='right',leaf_font_size = 10) 
    # adding/removing the axes
    ax1.set_xticks([])
    ax1.set_yticks([])

    # Compute and plot second dendrogram.
    ax2 = fig.add_axes([0.27,0.72,0.6,0.2])
    Z2 = hier.dendrogram(Y)
    ax2.set_xticks([])
    ax2.set_yticks([])

    #Compute and plot the heatmap
    axmatrix = fig.add_axes([0.27,0.1,0.6,0.6])
    idx1 = Z1['leaves']
    idx2 = Z2['leaves']
    D = scidist.squareform(data_dist)
    D = D[idx1,:]
    D = D[:,idx2]
    im = axmatrix.matshow(D, aspect='auto', origin='lower', cmap=cm.RdBu)
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])

    # Plot colorbar.
    axcolor = fig.add_axes([0.89,0.1,0.03,0.6])
    plt.colorbar(im, cax=axcolor)    
    plt.savefig(figname, dpi = 300,bbox_inches='tight')

    return None

def plot_dendrogram(Y, figname="match_dendrogram.png"):
    import scipy.cluster.hierarchy as hier
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.font_manager as font_manager
    
    # Font of figure
    path = '/Library/Frameworks/EPD64.framework/Versions/Current/lib/python2.7/\
site-packages/matplotlib/mpl-data/fonts/ttf/Helvetica.ttf'
    try:
        prop = font_manager.FontProperties(fname=path)
        mpl.rcParams['font.family'] = prop.get_name()
    except:
        pass
    
    # color
    mpl.rc('axes',edgecolor='k')
    colors = ['#1e90ff', '#ff69b4', '#32cd32', '#9370db', '#2f4f4f','#87ceeb',\
          '#ff7f00','#228b22','#b22222','#000000']
    hier.set_link_color_palette(colors)
    
    fig = plt.figure(figsize=(10,5))
    Z = hier.dendrogram(Y, leaf_font_size=10,leaf_rotation =90, color_threshold = 0.75)

    plt.ylim(0.0, 1.0)
    plt.savefig(figname, dpi = 300,bbox_inches='tight')
    
    
    

def fill_pocket(natom, pocket):
    """Convert Pocket into binary atom array
    
    :param natom: Number of atoms
    :param pocket: Instance of pockets
    :returns: ndarray([len(pocekt), natom])
    :rtype: ndarray
    
    .. note:: If atom j in pocket i, assign (i,j) = 1.
    
    """
    
    npocket = len(pocket)
    pocket_atom_array = np.zeros([npocket, natom])
    for i in range(npocket):
        for j in pocket[i].atoms:
            pocket_atom_array[i,j] = 1.
    return pocket_atom_array


def calc_simsa(a, b):
    """Calculate similarity based on two SA array
    
    """
    atf = (a > 0.)
    btf = (b > 0.)
    
    abtf = atf & btf
    
    anew = a * abtf
    bnew = b * abtf
    

    return 1 - sum(abs(a-b)) / (sum(a) + sum(b))
    #return (sum(anew) + sum(bnew)) / (sum(a) + sum(b))
    

def calc_psimsa(pocket_atom_sa_array):
    """Calculate pairwised similarity based on array
    
    """
    nsnap, natom = pocket_atom_sa_array.shape
    
    psimsa = []
    
    for i in range(0, nsnap-1,1):
        for j in range(i+1, nsnap,1):
            psimsa.append(calc_simsa(pocket_atom_sa_array[i], \
                pocket_atom_sa_array[j]))
    return 1. - np.array(psimsa)



def calc_intraIdx(npocket_list):
    """Get the intra pocket index to assign to large number
    
    
    """
    n = sum(npocket_list) - 1
    a = [0] + list(range(n,1,-1))
    a = np.cumsum(np.array(a))
    
    bnum = np.array(npocket_list) - 1
    bidx = []
    
    for i in bnum:
        while i >= 0:
            bidx.append(i)
            i = i - 1
    bidx = bidx[:-1]
    
    newidx = []
    
    for i,j in zip(a, bidx):
        for k in range(0,j,1):
            newidx.append(i+k)
    return newidx


def calc_gsim(matchpocketidx, ftdict, saidx = 5):
    """
    
    """
    
    gsim = []
    
    for i,j in enumerate(matchpocketidx):
        if len(j) == 1:
            gsim.append(0.0)
        else:
            nsnap = len(j)

            for k,idx in enumerate(j):
                if 0 == k:
                    atomarr = [ftdict[idx][saidx]]
                else:
                    atomarr = np.append(atomarr, [ftdict[idx][saidx]], 0)
            gsim.append( np.round(1- np.average(calc_psimsa(atomarr)),3))
    return gsim 


def genMatchPocket(pdblist, SAwt = False, metric = 'jaccard',cmethod = 'distance',\
                      intraIgnore = False, isplot = False, isdebug = False):
    """
    Generate Pocket Matching for a list of PDBs
    
    :param pdblist: list of PDBs
    
    """


    pdblist = checkPDB(pdblist)
    numsnap = len(pdblist)
    
    ftdict = {}
    surfdict = {}
    snapidx = []
    npocket_list = []
    
    pdb_lines = [line for line in open(pdblist[0])]
    tmppdb = readPDB(pdb_lines, isTER = True)
    natom = len(tmppdb.prot)
    
    #olddir = os.getcwd()
    os.chdir("Pocket_matching")
    olddir = os.getcwd()

    for i,pdbnm in enumerate(pdblist):

        os.chdir("AS" + str(i+1))

        #adding this to copy existing option and parameter files into the directory to apply to AlphaSpace FCTM
        #this should be fixed to be generalizable in case the user has changed the name of these files from the default
        if os.path.isfile('../../AS_option.txt'):
            os.system("cp ../../AS_option.txt .")
        if os.path.isfile('../../AS_param.txt'):
            os.system("cp ../../AS_param.txt .")
        pdb_lines = [line for line in open(pdbnm)]
        pdb, ss = runAlphaSpace(pdb_lines, default = False, isWrite=True)
        ssp = ss.sub_pockets
        
        npocket_list.append(len(ssp))
        x = fill_pocket(natom, ssp) ## pocket atom array 
        surfdict[(i)] = pdb.asa_prot
        y = x * np.array(pdb.asa_prot) ## pocket atom array with sa
       
        for j,p in enumerate(ssp):
            snapidx.append((i,j))
            ftdict[(i,j)] = [x[j], list(p.atm_cntr()), \
                list(p.vert_cntr()), p.space, p.score, y[j]]
        os.chdir(olddir)
    
    #hard coding this to be true for now
    SAwt = True
    if SAwt:
        saidx = 5 ## use ftdict[?][5]
    else:
        saidx = 0
    
    pocket_atom_array = []

    for i,j in enumerate(snapidx):
        #print ([ftdict[j][saidx]])
        if 0 == i:
            pocket_atom_array = [ftdict[j][saidx]]
        else:
            pocket_atom_array = np.append(pocket_atom_array, [ftdict[j][saidx]], 0)

    if [] == pocket_atom_array:
        print("No pockets to match")
        sys.exit()


    ## Clustering       
    ## dendrogram tree
    if SAwt:
        data_dist = calc_psimsa(pocket_atom_array)
    elif metric == 'euclidean':
        data_dist = scidist.pdist(pocket_atom_array)
    elif metric == 'jaccard':
        data_dist = scidist.pdist(pocket_atom_array, "jaccard")
    else:
        print ("Error: wrong clustering method")
    

    intraIgnore = True
    
    if intraIgnore:
        data_dist[calc_intraIdx(npocket_list)] = 10.
    
    zmat = hier.linkage(data_dist, method="average")
    
    if isplot:
        plot_dendrogram(zmat)


    # method of cut the tree
    if cmethod == 'distance':
        #clust_dist =  0.7*max(zmat[:,2])
        clust_dist = 0.75
        cluster = hier.fcluster(zmat,clust_dist,criterion='distance')
    elif cmethod == 'maxclust':
        clust_num = int(np.average(npocket_list))   
        cluster = hier.fcluster(zmat,clust_num,criterion='maxclust')
    else:
        print ("Error: wrong cutoff method")
    
    #print cluster
    #print hier.leaders(zmat, cluster)
    
    if isplot:
        #np.savetxt("AtomArray_dataforplot.dat",pocket_atom_array)
        plot_heatmap(pocket_atom_array, figname ='match_clust_heatmap.png',\
                     metric=metric)
    
    nmatchpocket = max(cluster)   # number of pocket by MD clustering
    matchpocketidx = [[] for i in range(nmatchpocket)]  # empty list
    
    for i,j in zip(snapidx, cluster):
        matchpocketidx[j-1].append(i)
    
    gsim = calc_gsim(matchpocketidx, ftdict)

    sortidx =  sorted(list(range(len(gsim))), key=lambda x:gsim[x], reverse=True)
    
    write_matchpocket(matchpocketidx, gsim, numsnap, sortidx,  npocket_list)
    return None



def is_tool(name):
    try:
        devnull = open(os.devnull)
        subprocess.Popen([name],stdout=devnull, stderr=devnull).communicate()
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            return False
    return True


def calc_cluster_vol(coord):
    np_coord = np.array(coord)
    rad = [1.8 for i in range(len(coord))]
    np_rad = np.array(rad)
    if np_coord.size:
        #return AS_fcn.estVolume(p_coord,p_rad,5000)
        return estVolume(np_coord,np_rad,10000)
    else:
        return 0.0


def get_LFC_coord(lig_pdb,clust_dist):
    lig_coords = []
    for line in lig_pdb:
        lig_coords.append(pdbinfo.coord(line))
    zmat = hier.linkage(lig_coords, method='average')
    #zmat = hier.linkage(lig_coords, method='complete')
    cluster = hier.fcluster(zmat,clust_dist,criterion='distance')
    #cluster = hier.fcluster(zmat,3.5,criterion='distance')

    clusters = [ [] for i in range(max(cluster)) ]
    fill_p_verts(clusters,cluster)

    frag_centroids = []

    for count,i in enumerate(clusters):
        frag_coords = []
        for idx in i:
            frag_coords.append(lig_coords[idx])
        frag_centroid = np.average(frag_coords,axis=0)

        frag_centroids.append(frag_centroid)

    return frag_centroids














