# AlphaSpace v1.0 - python classes


import os

import numpy as np
from AlphaSpaceLite import AS_fcn
from AlphaSpaceLite import pdbinfo


class readPDB:
    """
    This is a class to read PDB and save into 2 parts (receptor and binder) with two options:
    Options one: isTER = False, it will read binder as anything with different chain ID or "HETATM"
    Options two: isTER = True, it will read binder as anything after "TER"

    if isReverse = True, it will swap the receptor and the binder

    water, protons, Na+, and Cl- will be stripped, 
    if multiple rotameric states present, state "A" will be used   
    """

    def __init__(self, pdb_lines: object, isTER: object = True, isReverse: object = False, useNaccess: object = False, getFace: object = False) -> object:
        # import pdbinfo
        # import tm_fcns
        self.pdb_lines = pdb_lines
        self.isTER = isTER
        self.isReverse = isReverse

        if self.isTER:
            self.prot = []
            self.lig = []
            isProt = True
            for lines in pdb_lines:
                if lines[0:3] == "TER":
                    isProt = False
                    continue
                if pdbinfo.isAtom(lines) == 1 and pdbinfo.isHydrogen(lines) == 0 and pdbinfo.isWater(lines) == 0 \
                        and pdbinfo.isNa(lines) == 0 and pdbinfo.isCl(lines) == 0 and (pdbinfo.roti(lines) == ' ' or \
                                                                                                   pdbinfo.roti(
                                                                                                           lines) == 'A'):
                    if isProt:
                        self.prot.append(lines.rstrip())
                    else:
                        self.lig.append(lines.rstrip())
            self.cmplx = self.prot + self.lig
            if isProt:
                print("******")
                print("Notice: option 'use_TER = True' ... however no 'TER' line was found in input PDB file")
                print(">>All atoms are being treated as 'receptor'")
                print(
                    ">>Without 'binder' atoms, AlphaSpace cannot screen for pockets by interface or by binder contact")
                print("******\n")

        else:

            self.cmplx = [lines.rstrip() for lines in self.pdb_lines \
                          if
                          pdbinfo.isAtom(lines) == 1 and pdbinfo.isHydrogen(lines) == 0 and pdbinfo.isWater(lines) == 0 \
                          and pdbinfo.isNa(lines) == 0 and (pdbinfo.roti(lines) == ' ' or pdbinfo.roti(lines) == 'A')]
            self.chain_1 = pdbinfo.chid(self.cmplx[0])
            self.prot = [lines for lines in self.cmplx if pdbinfo.isATOM(lines) == 1 \
                         and pdbinfo.chid(lines) == self.chain_1]
            self.lig = [lines for lines in self.cmplx if pdbinfo.isHETATM(lines) == 1 \
                        or pdbinfo.chid(lines) != self.chain_1]
            self.cmplx = self.prot + self.lig

        if self.isReverse:
            tmp_prot = self.prot
            self.prot = self.lig
            self.lig = tmp_prot
            self.cmplx = self.prot + self.lig

        # later add output table for the interface SA info as another option
        if useNaccess:
            # self.asa_prot = tm_fcns.calc_asa(self.prot)
            self.asa_prot,self.asa_radius = AS_fcn.calc_asa2(self.prot)
            self.asa_radius2 = np.array(self.asa_radius) ** 2

            # add asa for binder
            self.asa_lig,self.asa_lig_radius = AS_fcn.calc_asa2(self.lig)
            self.asa_lig_radius2 = np.array(self.asa_lig_radius) ** 2

            # calculate and save various interface features
            if getFace:
                # self.asa_desolv, self.face_idx, self.face_SA, self.perc_face_np, self.np_face_SA = self.get_face_SAs_cut()

                # adding binder side face idx/asa lists to fcn "get_face_idx"
                self.face_idx,self.face_asa,self.face_lig_idx,self.face_lig_asa = self.get_face_idx()

        # EDIT 11/12/15 : fix for "use_naccess = False"
        else:
            self.asa_prot = [1.0 for i in range(len(self.prot))]
            self.asa_radius = [0.0 for i in range(len(self.prot))]
            self.asa_radius2 = [0.0 for i in range(len(self.prot))]

            self.face_asa = [0.0 for i in range(len(self.prot))]
            self.face_idx = [i for i in range(len(self.prot))]
            self.face_lig_asa = [0.0 for i in range(len(self.lig))]
            self.face_lig_idx = [i for i in range(len(self.lig))]

            coord_prot = [pdbinfo.coord(l) for l in self.prot]
            coord_prot = np.array(coord_prot)
            coord_lig = [pdbinfo.coord(l) for l in self.lig]
            coord_lig = np.array(coord_lig)

            if getFace:
                self.face_idx = []
                self.face_lig_idx = []
                for i,c_p in enumerate(coord_prot):
                    for j,c_l in enumerate(coord_lig):

                        if np.linalg.norm(c_p - c_l) <= 5.0:
                            self.face_asa[i] = 1.0
                            # adding this (and commenting out "break") to get lig side face
                            self.face_lig_asa[j] = 1.0
                            # break
                for i,val in enumerate(self.face_asa):
                    if val > 0.0:
                        self.face_idx.append(i)
                for i,val in enumerate(self.face_lig_asa):
                    if val > 0.0:
                        self.face_lig_idx.append(i)

    def get_lig_vol(self,cntct_only=False,cntct_dist=4.5):
        """
        Function: to calculate the volume of the binder
        Options: 'cntct_only' will evaluate only binder atoms within 'cntct_dist' of the receptor
        Note: the function in 'return' uses random sampling to estimate the volume ('10000' = # of sample points) 
        """
        lig_coord = [pdbinfo.coord(l) for l in self.lig]
        lig_coord = np.array(lig_coord)
        lig = self.lig
        if cntct_only:
            prot_coord = [pdbinfo.coord(l) for l in self.prot]
            prot_coord = np.array(prot_coord)
            lig_coord_tmp = []
            lig_tmp = []
            for i,atm_l in enumerate(lig_coord):
                for atm_p in prot_coord:
                    if np.linalg.norm(atm_l - atm_p) <= cntct_dist:
                        lig_coord_tmp.append(atm_l)
                        lig_tmp.append(self.lig[i])
                        break
            lig = lig_tmp
            lig_coord = np.array(lig_coord_tmp)

        lig_asa,lig_rad = AS_fcn.calc_asa2(lig)
        lig_rad = np.array(lig_rad)
        return AS_fcn.estVolume(lig_coord,lig_rad,20000)

    def get_lig_vol_alt(self,cntct_only=False,cntct_dist=1.6,prot_coord=[]):
        """
        Function: alternate binder volume function, here 'cntct_only' ignores the receptor, but
                  is passed the alpha-atom list and only lig atoms in cntct w alpha-atoms are evaluated
        """
        lig_coord = [pdbinfo.coord(l) for l in self.lig]
        lig_coord = np.array(lig_coord)
        lig = self.lig
        if cntct_only:
            # prot_coord = [pdbinfo.coord(l) for l in self.prot]
            prot_coord = np.array(prot_coord)
            lig_coord_tmp = []
            lig_tmp = []
            for i,atm_l in enumerate(lig_coord):
                for atm_p in prot_coord:
                    if np.linalg.norm(atm_l - atm_p) <= cntct_dist:
                        lig_coord_tmp.append(atm_l)
                        lig_tmp.append(self.lig[i])
                        break
            lig = lig_tmp
            lig_coord = np.array(lig_coord_tmp)

        lig_asa,lig_rad = AS_fcn.calc_asa2(lig)
        lig_rad = np.array(lig_rad)
        return AS_fcn.estVolume(lig_coord,lig_rad,10000)

    def get_asa(self):
        self.asa_prot = AS_fcn.calc_asa(self.prot)
        return self.asa_prot

    def get_surf_idx(self):
        """
        Function: return a list of the receptor surface atoms' indices
        """
        idx_list = []
        for i,j in enumerate(self.asa_prot):
            if j > 0.0:
                idx_list.append(i)
        return idx_list

    def get_face_idx(self):
        face_idx = []
        face_asa = []
        face_lig_idx = []
        face_lig_asa = []
        if self.lig:
            asa_cmplx = AS_fcn.calc_asa(self.cmplx)

            for idx,val in enumerate(self.asa_prot):
                # check to see if each atomistic ASA has been reduced in the complex
                if val > asa_cmplx[idx]:
                    face_asa.append(val - asa_cmplx[idx])
                    face_idx.append(idx)
                else:
                    face_asa.append(0.0)

            # 4-11-16: added this for the binder side face info
            for idx,val in enumerate(self.asa_lig):
                # check to see if each atomistic ASA has been reduced in the complex
                if val > asa_cmplx[len(self.asa_prot) + idx]:
                    face_lig_asa.append(val - asa_cmplx[len(self.asa_prot) + idx])
                    face_lig_idx.append(idx)
                else:
                    face_lig_asa.append(0.0)

        return face_idx,face_asa,face_lig_idx,face_lig_asa

    def get_face(self,asa2):
        """
        parameters: asa1: list of apo receptor asa's, asa2: list of complex asa's
        return: index list of atoms at the interface (ie: change in asa upon complexation is > 0)
        """
        idx_list = []
        for i,j in enumerate(self.asa_prot):
            diff = j - asa2[i]
            if diff > 0.0:
                idx_list.append(i)
        return idx_list

    def get_face_lig(self):
        """
        function: get the binder side indices that are at the face
        """
        asa_cmplx = np.array(AS_fcn.calc_asa(self.cmplx)[len(self.prot):])
        asa_lig = np.array(AS_fcn.calc_asa(self.lig))
        idx_list = []
        for i,j in enumerate(self.asa_lig):
            diff = j - asa_cmplx[i]
            if diff > 0.0:
                idx_list.append(i)
        return idx_list

    def get_face_SAs(self):
        """
        Function: calculate total face SA, pnp for face SA, nonpolar face SA
        """
        # if no binder, then no face SA
        if len(self.lig) == 0:
            return [0.0,'n/a',0.0]
        else:
            face_ASA = 0.0
            face_ASA_np = 0.0

            asa_cmplx = AS_fcn.calc_asa(self.cmplx)

            asa_prot = np.array(self.asa_prot)
            asa_cmplx = np.array(asa_cmplx[0:len(self.asa_prot)])
            diff_asa = asa_prot - asa_cmplx

            for i,val in enumerate(diff_asa):
                # check to see if each atomistic ASA has been reduced in the complex
                face_ASA += val
                if pdbinfo.atmn(self.prot[i]).strip().startswith('N') or \
                        pdbinfo.atmn(pdb[i]).strip().startswith('O') or \
                                pdbinfo.atmn(pdb[i]).strip() == 'SG':
                    continue
                else:
                    face_ASA_np += val

            face_perc_np = face_ASA_np / face_ASA

            return [face_ASA,face_perc_np,face_ASA_np]

    def get_lig_asas(self):
        """
        Function: get all ASA data for the binder broken down by tot,face,non-face
        """
        np_sa_lig = 0.0
        p_sa_lig = 0.0
        tot_sa_lig = 0.0
        np_sa_face = 0.0
        p_sa_face = 0.0
        tot_sa_face = 0.0
        np_sa_nonf = 0.0
        p_sa_nonf = 0.0
        tot_sa_nonf = 0.0

        asa_cmplx = np.array(AS_fcn.calc_asa(self.cmplx)[len(self.prot):])
        asa_lig = np.array(AS_fcn.calc_asa(self.lig))

        # print(len(asa_lig))
        # print(len(asa_cmplx))
        # quit()

        asa_face = asa_lig - asa_cmplx

        pnp_lig = 0.0
        pnp_face = 0.0
        pnp_nonf = 0.0

        # tot_lig_sa = np.sum(asa_lig)

        for i,atm in enumerate(self.lig):
            tot_sa_nonf += asa_cmplx[i]
            tot_sa_face += asa_face[i]
            tot_sa_lig += asa_lig[i]
            if pdbinfo.atmn(atm).strip().startswith('C'):
                np_sa_nonf += asa_cmplx[i]
                np_sa_face += asa_face[i]
                np_sa_lig += asa_lig[i]
            # else:
            #     p_sa_nonf += asa_cmplx[i]
            #     p_sa_face += asa_face[i]
            #     p_sa_lig += asa_lig[i]

            p_sa_nonf = tot_sa_nonf - np_sa_nonf
            p_sa_face = tot_sa_face - np_sa_face
            p_sa_lig = tot_sa_lig - np_sa_lig

        if tot_sa_lig > 0.0:
            pnp_lig = np_sa_lig / tot_sa_lig
        else:
            pnp_lig = 0.0
        if tot_sa_face > 0.0:
            pnp_face = np_sa_face / tot_sa_face
        else:
            pnp_face = 0.0
        if tot_sa_nonf > 0.0:
            pnp_nonf = np_sa_nonf / tot_sa_nonf
        else:
            pnp_nonf = 0.0

        perc_cov = tot_sa_face / tot_sa_lig
        return [perc_cov,tot_sa_lig,pnp_lig,np_sa_lig,p_sa_lig,tot_sa_face,pnp_face,np_sa_face, \
                p_sa_face,tot_sa_nonf,pnp_nonf,np_sa_nonf,p_sa_nonf]

    def get_sum_SA(self):
        """
        Function: get total interface SA (pocket SA + binder SA)
        """
        asa_cmplx = np.array(AS_fcn.calc_asa(self.cmplx))
        asa_lig = np.array(AS_fcn.calc_asa(self.lig))
        asa_prot = np.array(self.asa_prot)

        # asa_apo_cmplx = np.array()

        asa_apo_cmplx = np.concatenate((asa_prot,asa_lig))

        diff_asa = asa_apo_cmplx - asa_cmplx
        sum_SA = np.sum(diff_asa)
        pock_SA = np.sum(diff_asa[:len(self.prot)])
        lig_SA = np.sum(diff_asa[len(self.prot):])
        lig_totSA = np.sum(asa_lig)
        perc_cov = lig_SA / lig_totSA
        struct_score = perc_cov * sum_SA
        return [sum_SA,pock_SA,lig_SA,perc_cov,struct_score]


class Snapshot:
    """
    Class: to contain, for a single snapshot or pdb structure, the atom_xyz
           information for the receptor and the binder, the data from the Voronoi
           tessellation, and the list of Pockets
    """

    def __init__(self,pdb,vertices,simplices,radii,sub_pockets=[],idx=None):
        self.prot_pdb = pdb.prot
        self.lig_pdb = pdb.lig
        if hasattr(pdb,'asa_prot'):
            self.asa_prot = pdb.asa_prot
            self.asa_radius = pdb.asa_radius
            self.asa_radius2 = pdb.asa_radius2

        if hasattr(pdb,'face_idx'):
            self.face_idx = pdb.face_idx
            self.face_asa = pdb.face_asa

        # 4-11-16: adding this for binder side
        if hasattr(pdb,'asa_lig'):
            self.asa_lig = pdb.asa_prot
            self.asa_lig_radius = pdb.asa_lig_radius
            self.asa_lig_radius2 = pdb.asa_lig_radius2

        if hasattr(pdb,'face_idx'):
            self.face_idx = pdb.face_idx
            self.face_asa = pdb.face_asa
        if hasattr(pdb,'face_lig_idx'):
            self.face_lig_idx = pdb.face_lig_idx
            self.face_lig_asa = pdb.face_lig_asa
        self.vertices = vertices
        self.simplices = simplices
        self.radii = radii

        self.sub_pockets = sub_pockets
        for p in self.sub_pockets:
            # link each pocket to its residue Snapshot
            p.set_ss(self)

        self.p_asa_prot = []
        self.cntct_arr = []

        self.lig_resid_resn = {}

        # beta-cluster is a list of beta-atoms with vert and atoms referencing
        # these new beta vert, simp, and rad lists
        self.beta_clusters = []
        # self.vert_beta = []
        # self.simp_beta = []
        # self.rad_beta = []

        # self.beta_cluster_frag = []
        # self.vert_beta_frag = []
        # self.simp_beta_frag = []
        # self.rad_beta_frag = []

        self.idx = idx

    # #This the original "fill_beta" before going from full face betas to fragment-centric betas
    # def fill_beta(self,clust_cutoff):
    #     import scipy.cluster.hierarchy as hier
    #     import AS_fcn

    #     vert_beta = []
    #     simp_beta = []
    #     rad_beta = []
    #     for p in self.sub_pockets:
    #         for v in p.vertices_coords:
    #             vert_beta.append(self.vertices_coords[v])
    #             simp_beta.append(self.simplex_idx_list[v])
    #             rad_beta.append(self.radii[v])

    #     self.vert_beta = vert_beta
    #     self.simp_beta = simp_beta
    #     self.rad_beta = rad_beta

    #     if len(vert_beta) > 1:

    #         #do the clustering  
    #         zmat = hier.linkage(vert_beta, method="complete")
    #         cluster = hier.fcluster(zmat,clust_cutoff,criterion='distance')

    #     elif len(vert_beta) == 1:
    #         cluster = [1]

    #     else: 
    #         return 

    #     b_verts = [ [] for i in range(max(cluster)) ]
    #     b_atms = [ [] for i in range(max(cluster)) ]

    #     AS_fcn.fill_p_verts(b_verts,cluster)
    #     AS_fcn.fill_p_atms(b_atms,b_verts,simp_beta)

    #     beta_cluster = []
    #     for i,b in enumerate(b_verts):
    #         beta_cluster.append(Beta_atom(b,b_atms[i],self))

    #     self.beta_cluster = beta_cluster

    #     return


    def fill_beta_frag(self,clust_cutoff):
        import scipy.cluster.hierarchy as hier
        from retired.AlphaSpace import AS_fcn

        for p in self.sub_pockets:

            vert_beta = []
            simp_beta = []
            rad_beta = []

            for v in p.vertices:
                vert_beta.append(self.vertices[v])
                simp_beta.append(self.simplices[v])
                rad_beta.append(self.radii[v])

            # self.vert_beta_frag.append(vert_beta)
            # self.simp_beta_frag.append(simp_beta)
            # self.rad_beta_frags.append(rad_beta)

            # self.beta_clusters.append(Beta_clust(vert_beta,simp_beta,rad_beta))
            curr_beta_clust = Beta_clust(vert_beta,simp_beta,rad_beta)

            if len(vert_beta) > 1:

                # do the clustering
                zmat = hier.linkage(vert_beta,method="complete")
                cluster = hier.fcluster(zmat,clust_cutoff,criterion='distance')

            elif len(vert_beta) == 1:
                cluster = [1]

            else:
                return

            b_verts = [[] for i in range(max(cluster))]
            b_atms = [[] for i in range(max(cluster))]

            AS_fcn.fill_p_verts(b_verts, cluster)
            AS_fcn.fill_p_atms(b_atms, b_verts, simp_beta)

            beta_cluster = []
            for i,b in enumerate(b_verts):
                curr_beta_clust.beta_atoms.append(Beta_atom(b,b_atms[i],self,curr_beta_clust))

            self.beta_clusters.append(curr_beta_clust)

        return

    # #This the original "write_beta" before going from full face betas to fragment-centric betas
    # def write_beta(self,class_cutoff,dir_name=".",hit_dist=1.6,high_cutoff=20.0,mid_cutoff=10.0,get_beta_vol=False):

    #     import re

    #     combine_atoms = []
    #     for b in self.beta_cluster:
    #         combine_atoms += [atm for atm in b.atoms if atm not in combine_atoms]
    #     combine_atoms.sort()

    #     #store trajectory to calculate centroids for beta-cluster and all pocket atoms
    #     atom_coord = []
    #     beta_coord = []
    #     beta_coord_occ = []

    #     f_head = open(dir_name + '/pockets/beta.info','w')
    #     header = ['beta-score','#all_atm', '#pol_atm', '#chg_atm', '#pos_atm', '#neg_atm']
    #     print(*header,sep='\t',file=f_head)
    #     print("\n\nalso: beta_atom resid is 1 (all-nonpolar), 2 (polar atom, does not include backbone), 3 (charged atom)",file=f_head)
    #     print("\nbeta.vol file (if exists) lists 1. full beta-cluster volume; 2. occupied beta-cluster volume",file=f_head)
    #     f_head.close()
    #     f = open(dir_name + '/pockets/beta.pdb','w')
    #     for a_idx in combine_atoms:
    #         f.write(self.prot_pdb[a_idx][0:60])

    #         atom_coord.append(pdbinfo.coord(self.prot_pdb[a_idx]))

    #         #write atomic ACSA 
    #         f.write('{0:6.2f}'.format(self.asa_prot[a_idx]) + '\n')
    #     pos_atoms = ['NE','NH1','NH2','NZ']
    #     asp_atoms = ['OD1','OD2']
    #     glu_atoms = ['OE1','OE2']
    #     hip_atoms = ['ND1', 'NE2']
    #     for i,b in enumerate(self.beta_cluster):

    #         #collect all the classifications
    #         if b.get_occ_status(hit_dist):
    #             occ_str = "BAO"
    #         else:
    #             occ_str = "BAU"
    #         #EDIT 11-17-15: changing this to be space ('1') instead of score ('0') in order to take account of polar space too...
    #         beta_score = b.beta_score_space()[1]
    #         if beta_score >= high_cutoff:
    #             chain = 'H'
    #         elif beta_score >= mid_cutoff:
    #             chain = 'M'
    #         else:
    #             chain = 'L'

    #         b_center = b.vert_cntr()
    #         beta_coord.append(b_center)
    #         if occ_str == "BAO":
    #             beta_coord_occ.append(b_center)
    #         atm_cnt = 0
    #         pol_cnt = 0
    #         chg_cnt = 0
    #         pos_cnt = 0
    #         neg_cnt = 0
    #         for a in b.atoms:
    #             line = self.prot_pdb[a]
    #             if np.linalg.norm(np.array(pdbinfo.coord(line)) - b_center) <= class_cutoff:
    #                 atm_cnt += 1
    #                 #if pdbinfo.atmn(line).strip().startswith('N') or \
    #                 # pdbinfo.atmn(line).strip().startswith('O') or \
    #                 # pdbinfo.atmn(line).strip() == 'SG':

    #                 #changing this to not classify back bone atoms as polar****
    #                 if re.match("(N\w+)", pdbinfo.atmn(line).strip()) or \
    #                  re.match("(O\w+)", pdbinfo.atmn(line).strip()) or \
    #                  pdbinfo.atmn(line).strip() == 'SG':


    #                     pol_cnt += 1
    #                     if pdbinfo.atmn(line).strip() in pos_atoms or (pdbinfo.atmn(line).strip() in hip_atoms \
    #                      and pdbinfo.resn(line).strip() == 'HIP'):
    #                         chg_cnt += 1
    #                         pos_cnt += 1
    #                     if (pdbinfo.atmn(line).strip() in asp_atoms and pdbinfo.resn(line).strip() == 'ASP') \
    #                      or  (pdbinfo.atmn(line).strip() in glu_atoms and pdbinfo.resn(line).strip() == 'GLU'):
    #                         chg_cnt += 1
    #                         neg_cnt += 1

    #         if chg_cnt > 0:
    #             res_id = '3'
    #         #note setting pol threshhold to "1" to make it a little more selective (could also try not counting backbone atoms)
    #         #so I did set pol to 0 and remove backbone atoms for now, later will only remove if secondary structure is sheet or helix
    #         elif pol_cnt > 0:
    #             res_id = '2'
    #         elif atm_cnt > 0:
    #             res_id = '1'
    #         else:
    #             res_id = '0' 


    #         #f.write(b.vert_cntr)
    #         f.write('ATOM   ' + str(i).rjust(4) + '  ' + occ_str + ' BAC ')
    #         f.write(chain)
    #         f.write('   ' + res_id + '    ')

    #         for j in b_center:
    #             print('{0:8.3f}'.format(j),end="",file=f)
    #         #print('{0:12.2f}\n'.format(self.rad_beta[v_idx]),end="",file=f)

    #         #scr = b.beta_score_space()[0]
    #         #f.write('{0:8.3f}'.format(scr),end="",file=f)
    #         f.write('{0:6.1f}'.format(beta_score) + '  ')

    #         bins = [atm_cnt,pol_cnt,chg_cnt,pos_cnt,neg_cnt]
    #         print(*bins,sep='   ',file=f)
    #         #f.write("\n")

    #     if get_beta_vol:
    #         beta_vol = AS_fcn.calc_cluster_vol(beta_coord)
    #         beta_vol_occ = AS_fcn.calc_cluster_vol(beta_coord_occ)
    #         f_vol = open(dir_name + '/pockets/beta.vol','w')
    #         vol_data = [beta_vol,beta_vol_occ]
    #         print(*vol_data,sep='\t',file=f_vol)
    #         f_vol.close()

    #     #adding this "if" to match what I see in v1.1.4, don't remeber how these diverged?
    #     if atom_coord and beta_coord:

    #         atom_centroid = np.average(atom_coord,axis=0)
    #         beta_centroid = np.average(beta_coord,axis=0)

    #         f.write('ATOM    100  PBC PBC X   0    ')
    #         for j in atom_centroid:
    #                 print('{0:8.3f}'.format(j),end="",file=f)
    #         f.write('\n')
    #         f.write('ATOM    101  BCC BCC X   0    ')
    #         for j in beta_centroid:
    #                 print('{0:8.3f}'.format(j),end="",file=f)
    #     f.close()
    #     return


    def write_beta_frag(self,class_cutoff,dir_name=".",hit_dist=1.6,core_p_cutoff=100.0,aux_p_cutoff=30.0, \
                        high_cutoff=20.0,mid_cutoff=10.0,get_beta_vol=False):

        import re

        # write the beta information file
        f_head = open(dir_name + '/pockets/beta.info','w')
        header = ['beta-np-space','beta-pol-space','beta-total-space','#all_atm','#pol_atm','#chg_atm','#pos_atm',
                  '#neg_atm']
        print(*header,sep='\t',file=f_head)
        print(
            "\n\nalso: beta_atom resid is 1 (all-nonpolar), 2 (polar atom, does not include backbone), 3 (charged atom)",
            file=f_head)
        print("\nbeta.vol file (if exists) lists 1. full beta-cluster volume; 2. occupied beta-cluster volume",
              file=f_head)
        f_head.close()

        f_total_beta = open(dir_name + '/pockets/beta.pdb','w')
        tot_combine_atoms = []
        for index,b_c in enumerate(self.beta_clusters):
            combine_atoms = []
            for b in b_c.beta_atoms:
                combine_atoms += [atm for atm in b.atoms if atm not in combine_atoms]
                tot_combine_atoms += [atm for atm in b.atoms if atm not in tot_combine_atoms]
            combine_atoms.sort()

            # store trajectory to calculate centroids for beta-cluster and all pocket atoms
            atom_coord = []
            beta_coord = []
            beta_coord_occ = []

            f = open(dir_name + '/pockets/' + str(index).zfill(3) + '_beta.pdb','w')
            for a_idx in combine_atoms:
                f.write(self.prot_pdb[a_idx][0:60])

                atom_coord.append(pdbinfo.coord(self.prot_pdb[a_idx]))

                # write atomic ACSA
                f.write('{0:6.2f}'.format(self.asa_prot[a_idx]) + '\n')

            pos_atoms = ['NE','NH1','NH2','NZ']
            asp_atoms = ['OD1','OD2']
            glu_atoms = ['OE1','OE2']
            hip_atoms = ['ND1','NE2']
            for i,b in enumerate(b_c.beta_atoms):

                # collect all the classifications
                if b.get_occ_status(hit_dist):
                    occ_str = "BAO"
                else:
                    occ_str = "BAU"
                # EDIT 11-17-15: changing this to be space ('1') instead of score ('0') in order to take account of polar space too...
                # EDIT 2-1-15: changing back because now I want score... and I will treat polar separately... not sure about this...
                beta_data = b.beta_score_space()

                beta_space = beta_data[2]
                beta_pol_space = beta_data[1]
                beta_score = beta_data[0]
                if beta_score >= high_cutoff:
                    chain = 'H'
                elif beta_score >= mid_cutoff:
                    chain = 'M'
                else:
                    chain = 'L'

                b_center = b.vert_cntr()
                beta_coord.append(b_center)
                if occ_str == "BAO":
                    beta_coord_occ.append(b_center)
                atm_cnt = 0
                pol_cnt = 0
                chg_cnt = 0
                pos_cnt = 0
                neg_cnt = 0
                for a in b.atoms:
                    line = self.prot_pdb[a]
                    if np.linalg.norm(np.array(pdbinfo.coord(line)) - b_center) <= class_cutoff:
                        atm_cnt += 1
                        # if pdbinfo.atmn(line).strip().startswith('N') or \
                        # pdbinfo.atmn(line).strip().startswith('O') or \
                        # pdbinfo.atmn(line).strip() == 'SG':

                        # changing this to not classify back bone atoms as polar****
                        if re.match("(N\w+)",pdbinfo.atmn(line).strip()) or \
                                re.match("(O\w+)",pdbinfo.atmn(line).strip()) or \
                                        pdbinfo.atmn(line).strip() == 'SG':

                            pol_cnt += 1
                            if pdbinfo.atmn(line).strip() in pos_atoms or (pdbinfo.atmn(line).strip() in hip_atoms \
                                                                                   and pdbinfo.resn(
                                        line).strip() == 'HIP'):
                                chg_cnt += 1
                                pos_cnt += 1
                            if (pdbinfo.atmn(line).strip() in asp_atoms and pdbinfo.resn(line).strip() == 'ASP') \
                                    or (pdbinfo.atmn(line).strip() in glu_atoms and pdbinfo.resn(
                                            line).strip() == 'GLU'):
                                chg_cnt += 1
                                neg_cnt += 1

                if chg_cnt > 0:
                    res_id = '3'
                # note setting pol threshhold to "1" to make it a little more selective (could also try not counting backbone atoms)
                # so I did set pol to 0 and remove backbone atoms for now, later will only remove if secondary structure is sheet or helix
                elif pol_cnt > 0:
                    res_id = '2'
                elif atm_cnt > 0:
                    res_id = '1'
                else:
                    res_id = '0'


                    # f.write(b.vert_cntr)
                f.write('ATOM   ' + str(i).rjust(4) + '  ' + occ_str + ' BAC ')
                f.write(chain)
                f.write('   ' + res_id + '    ')

                # write to full beta face
                f_total_beta.write('ATOM   ' + str(i).rjust(4) + '  ' + occ_str + ' BAC ')
                f_total_beta.write(chain)
                f_total_beta.write('   ' + res_id + '    ')

                for j in b_center:
                    print('{0:8.3f}'.format(j),end="",file=f)
                    print('{0:8.3f}'.format(j),end="",file=f_total_beta)
                # print('{0:12.2f}\n'.format(self.rad_beta[v_idx]),end="",file=f)

                # scr = b.beta_score_space()[0]
                # f.write('{0:8.3f}'.format(scr),end="",file=f)
                f.write('{0:6.1f}'.format(beta_score))
                f_total_beta.write('{0:6.1f}'.format(beta_score))

                f.write('{0:6.1f}'.format(beta_pol_space))
                f_total_beta.write('{0:6.1f}'.format(beta_pol_space))

                f.write('{0:6.1f}'.format(beta_space) + '  ')
                f_total_beta.write('{0:6.1f}'.format(beta_space) + '  ')

                bins = [atm_cnt,pol_cnt,chg_cnt,pos_cnt,neg_cnt]
                print(*bins,sep='   ',file=f)
                print(*bins,sep='   ',file=f_total_beta)
                # f.write("\n")

            if get_beta_vol:
                beta_vol = AS_fcn.calc_cluster_vol(beta_coord)
                beta_vol_occ = AS_fcn.calc_cluster_vol(beta_coord_occ)
                f_vol = open(dir_name + '/pockets/beta.vol','w')
                vol_data = [beta_vol,beta_vol_occ]
                print(*vol_data,sep='\t',file=f_vol)
                f_vol.close()

            # adding this "if" to match what I see in v1.1.4, don't remeber how these diverged?
            if atom_coord and beta_coord:

                if self.sub_pockets[index].score >= core_p_cutoff:
                    chain_BCC = 'H'
                elif self.sub_pockets[index].score >= aux_p_cutoff:
                    chain_BCC = 'M'
                else:
                    chain_BCC = 'L'

                atom_centroid = np.average(atom_coord,axis=0)
                beta_centroid = np.average(beta_coord,axis=0)

                f.write('ATOM    100  PBC PBC ' + chain_BCC + '   0    ')
                for j in atom_centroid:
                    print('{0:8.3f}'.format(j),end="",file=f)
                f.write('\n')
                f.write('ATOM    101  BCC BCC ' + chain_BCC + '   0    ')
                for j in beta_centroid:
                    print('{0:8.3f}'.format(j),end="",file=f)
            f.close()

        tot_combine_atoms_filter = []
        tot_combine_atoms_filter += [atm for atm in tot_combine_atoms if atm not in tot_combine_atoms_filter]
        tot_combine_atoms_filter.sort()
        for a_idx in tot_combine_atoms_filter:
            f_total_beta.write(self.prot_pdb[a_idx][0:60])

            # atom_coord.append(pdbinfo.coord(self.prot_pdb[a_idx]))

            # write atomic ACSA
            f_total_beta.write('{0:6.2f}'.format(self.asa_prot[a_idx]) + '\n')
        f_total_beta.close()

        return

    def write_cntct_NEW_2(self,cntct_dict,e_cntct_dict,dir_name=".",write_header=True):

        f = open(dir_name + '/table_cntctscore_new_2.dat','w')
        header = ["resID","total_score","np_score","pol_score","e_score","unocc_space"]
        if write_header:
            print(*header,sep='\t',file=f)

        face_total_score = 0.0
        face_np_score = 0.0
        face_pol_score = 0.0
        face_e_score = 0.0
        face_unocc_space = 0.0

        # make non-redundant list of all resID keys that have any cntct score
        keys = [key for key in cntct_dict]
        for key in e_cntct_dict:
            if key not in keys:
                keys.append(key)
        for key in sorted(keys):
            total_score = 0.0
            if key in cntct_dict:
                resid = key + '-' + cntct_dict[key]["resn"]

                total_score += cntct_dict[key]["match_occ_space"]
                np_score = "{0}".format(int(round(cntct_dict[key]["np_match_occ_space"],0)))
                pol_score = "{0}".format(int(round(cntct_dict[key]["pol_match_occ_space"],0)))
                unocc_space = "{0}".format(int(round(cntct_dict[key]["unocc_space"],0)))

                face_np_score += cntct_dict[key]["np_match_occ_space"]
                face_pol_score += cntct_dict[key]["pol_match_occ_space"]
                face_unocc_space += cntct_dict[key]["unocc_space"]

            else:
                np_score = 0
                pol_score = 0
                unocc_space = 0

            if key in e_cntct_dict:
                resid = key + '-' + e_cntct_dict[key]["resn"]

                total_score += e_cntct_dict[key]["match_occ_space"]
                e_score = "{0}".format(int(round(e_cntct_dict[key]["match_occ_space"],0)))
                face_e_score += e_cntct_dict[key]["match_occ_space"]

            else:
                e_score = 0

            if total_score == 0.0:
                continue

            tot_score = "{0}".format(int(round(total_score,0)))
            face_total_score += total_score

            data = [resid,tot_score,np_score,pol_score,e_score,unocc_space]
            print(*data,sep='\t',file=f)

        resid = "FACE"
        tot_score = "{0}".format(int(round(face_total_score,0)))
        np_score = "{0}".format(int(round(face_np_score,0)))
        pol_score = "{0}".format(int(round(face_pol_score,0)))
        e_score = "{0}".format(int(round(face_e_score,0)))
        unocc_space = "{0}".format(int(round(face_unocc_space,0)))

        data = [resid,tot_score,np_score,pol_score,e_score,unocc_space]
        print(*data,sep='\t',file=f)

        return

    def write_cntct_NEW(self,cntct_dict,dir_name=".",write_header=True):
        """
        function to print out the cntct score data for all residues in "lig_resids" as a table
        """
        f = open(dir_name + '/table_cntctscore_new.dat','w')

        header = ["resID","Pockets","total_space","%_nonpolar","occ_space","%_occupied"]

        if write_header:
            print(*header,sep='\t',file=f)

        face_space = 0.0
        face_occ_space = 0.0
        face_np_space = 0.0

        for key in sorted(cntct_dict.keys()):
            resid = key + '-' + cntct_dict[key]["resn"]
            pockets = ','.join(map(str,cntct_dict[key]["pockets"]))

            # score = "{0}".format(int(round(cntct_dict[key]["match_occ_space"],0)))

            if cntct_dict[key]["space"] > 0.0:
                space = "{0}".format(int(round(cntct_dict[key]["space"],0)))
                perc_occ = "{0}".format(
                    int(round(cntct_dict[key]["occ_space"] / cntct_dict[key]["space"],2) * 100.0)) + '%'
                perc_np = "{0}".format(
                    int(round(cntct_dict[key]["np_space"] / cntct_dict[key]["space"],2) * 100.0)) + '%'
                occ_space = "{0}".format(int(round(cntct_dict[key]["occ_space"],0)))

                face_space += cntct_dict[key]["space"]
                face_np_space += cntct_dict[key]["np_space"]
                face_occ_space += cntct_dict[key]["occ_space"]
            else:
                space = "0"
                perc_occ = '0%'
                perc_np = '0%'
                occ_space = "0"


                # for line in cntct_arr:
                #     resid = line[0]
                #     space = "{0}".format(int(round(line[2],0)))
                #     perc_np = "{0}".format(int(round(line[4],2) * 100.0)) + '%'
                #     occ_space = "{0}".format(int(round(line[5],0)))
                #     match_occ_space = "{0}".format(int(round(line[7],0)))
                #     perc_occ = "{0}".format(int(round(line[6],2) * 100.0)) + '%'
                #     perc_match_occ = "{0}".format(int(round(line[8],2) * 100.0)) + '%'

            # data = [resid,space,perc_np,occ_space,match_occ_space,perc_occ,perc_match_occ]

            data = [resid,pockets,space,perc_np,occ_space,perc_occ]

            print(*data,sep='\t',file=f)

        resid = "FACE"
        space = "{0}".format(int(round(face_space,0)))
        pockets = "n/a"
        if face_space > 0.0:
            perc_np = "{0}".format(int(round(face_np_space / face_space,2) * 100.0)) + '%'
        else:
            perc_np = "0.0%"
        occ_space = "{0}".format(int(round(face_occ_space,0)))
        if face_space > 0.0:
            perc_occ = "{0}".format(int(round(face_occ_space / face_space,2) * 100.0)) + '%'
        else:
            perc_occ = "0.0%"

        data = [resid,pockets,space,perc_np,occ_space,perc_occ]
        print(*data,sep='\t',file=f)

        return

    def write_cntct(self,cntct_arr,dir_name=".",write_header=True):
        """
        function to print out the cntct score data for all residues in "lig_resids" as a table
        """
        f = open(dir_name + '/table_cntctscore.dat','w')

        header = ["resID","total-score","c-p-count","c-score","c-pnp","c-space","e-p-count","e-score","e-pnp","e-space"]
        if write_header:
            print(*header,sep='\t',file=f)
        for line in cntct_arr:
            resid = line[0]
            tot_score = "{0}".format(int(round(line[1],0)))
            c_p_count = line[2]
            c_score = "{0}".format(int(round(line[3],0)))
            c_pnp = "{0}".format(int(round(line[4],2) * 100.0)) + '%'
            c_space = "{0}".format(int(round(line[5],0)))
            e_p_count = line[6]
            e_score = "{0}".format(int(round(line[7],0)))
            e_pnp = "{0}".format(int(round(line[8],2) * 100.0)) + '%'
            e_space = "{0}".format(int(round(line[9],0)))

            data = [resid,tot_score,c_p_count,c_score,c_pnp,c_space,e_p_count,e_score,e_pnp,e_space]
            print(*data,sep='\t',file=f)

    def get_e_cntct_NEW(self,pdb_lines,option_dict,param_dict):
        cntct_data = {}
        data_template = {"occ_space": 0.0,"match_occ_space": 0.0,"np_match_occ_space": 0.0,"pol_match_occ_space": 0.0}
        skip_reverse = False

        # make sure the "binder" has at least 5 atoms in order to run the reverse surface mapping
        if len(self.lig_pdb) >= 5:

            # if already did "reverse" mapping, then do regular for inverse mapping
            # if did regular, now do reverse for inverse mapping
            if option_dict['do_reverse']:
                envelope_do_reverse = False
            else:
                envelope_do_reverse = True

            pdb_lig = readPDB(pdb_lines,isTER=option_dict['use_TER'],isReverse=envelope_do_reverse, \
                              useNaccess=option_dict['use_naccess'],getFace=option_dict['get_face'])
            env_face_idx = self.face_idx

            # perform the Voronoi tessellation
            cutvert_lig,cutsimp_lig,cutrad_lig = AS_fcn.genVoro(pdb_lig,min_r=param_dict['min_r'],
                                                                max_r=param_dict['max_r'])

            if len(cutvert_lig) >= 1:

                # log the start of mapping
                AS_fcn.log("\nmapping the interface...\n",option_dict,to_screen=True,to_file=False)

                # perform fragmnet-centric topographical mapping (FCTM)
                ss_lig = AS_fcn.genPocket_wVoro(pdb_lig,cutvert_lig,cutsimp_lig,cutrad_lig,
                                                is_use_naccess=option_dict['use_naccess'], \
                                                is_get_face=option_dict['get_face'],clust_dist=param_dict['clust_dist'], \
                                                min_num_sph=param_dict['min_num_alph'],
                                                screen_by_face_atoms=option_dict['screen_by_face'], \
                                                screen_face_perc=option_dict['screen_face_perc'],
                                                hit_dist=param_dict['hit_dist'], \
                                                screen_by_score=option_dict['screen_by_score'],
                                                min_score=option_dict['min_score'], \
                                                screen_out_subsurf=option_dict['screen_out_subsurf'],
                                                subsurf_perc=option_dict['max_desolv_perc'], \
                                                is_use_exact_ACSA=option_dict['use_exact_ACSA'],
                                                screen_by_lig_contact=True, \
                                                expand_around_cntct=False,is_use_NP_wt=option_dict['use_NP_wt'],
                                                clust_method='average', \
                                                screen_by_resid=False,resid_list=[],screen_by_seam=False,
                                                expand_around_seam=False)
            else:
                skip_reverse = True
                print("\nSkipping mapping of the 'binder' surface due to too few binder voronoi vertices_coords...\n")

        else:
            skip_reverse = True
            print("\nSkipping mapping of the 'binder' surface due to too few binder atoms...\n")

        if skip_reverse:
            return cntct_data

        pol_face_pdb = [self.prot_pdb[pdb_idx] for pdb_idx in self.face_idx if
                        not pdbinfo.atmn(self.prot_pdb[pdb_idx]).strip().startswith('C') and \
                        not pdbinfo.atmn(self.prot_pdb[pdb_idx]).strip() == 'SD']

        pol_coord_list = [pdbinfo.coord(line) for line in pol_face_pdb]
        pol_coord_list = np.array(pol_coord_list)

        np_face_pdb = [self.prot_pdb[pdb_idx] for pdb_idx in self.face_idx if
                       pdbinfo.atmn(self.prot_pdb[pdb_idx]).strip().startswith('C') or \
                       pdbinfo.atmn(self.prot_pdb[pdb_idx]).strip() == 'SD']

        np_coord_list = [pdbinfo.coord(line) for line in np_face_pdb]
        np_coord_list = np.array(np_coord_list)

        for p in ss_lig.sub_pockets:

            for v in p.vertices:
                pol_hit = False
                np_hit = False
                for coord in pol_coord_list:
                    if (np.linalg.norm(coord - ss_lig.vertices[v]) <= param_dict['hit_dist']):
                        pol_hit = True
                        break
                for coord in np_coord_list:
                    if (np.linalg.norm(coord - ss_lig.vertices[v]) <= param_dict['hit_dist']):
                        np_hit = True
                        break

                if pol_hit or np_hit:
                    # calculate the alpha-space
                    points = []
                    for j in ss_lig.simplices[v]:
                        points.append(pdbinfo.coord(ss_lig.prot_pdb[j]))
                    points = np.array(points)
                    vol = AS_fcn.calc_simp_vol(points)

                    curr_v_tot_sa = 0.0
                    for a_idx in ss_lig.simplices[v]:
                        curr_v_tot_sa += ss_lig.asa_prot[a_idx]
                    for a_idx in ss_lig.simplices[v]:
                        curr_resid = ss_lig.prot_pdb[a_idx][23:26].strip()
                        curr_atmn = pdbinfo.atmn(ss_lig.prot_pdb[a_idx]).strip()
                        curr_resid = str(curr_resid).zfill(3)
                        if curr_resid not in cntct_data:
                            cntct_data[curr_resid] = dict(data_template)
                            cntct_data[curr_resid]["resn"] = pdbinfo.resn(ss_lig.prot_pdb[a_idx])
                        curr_space = (ss_lig.asa_prot[a_idx] / curr_v_tot_sa) * vol
                        cntct_data[curr_resid]["occ_space"] += curr_space

                        if curr_atmn.startswith('C') or curr_atmn == 'SD':
                            if np_hit:
                                cntct_data[curr_resid]["match_occ_space"] += curr_space
                                cntct_data[curr_resid]["np_match_occ_space"] += curr_space
                        else:
                            if pol_hit:
                                cntct_data[curr_resid]["match_occ_space"] += curr_space
                                cntct_data[curr_resid]["pol_match_occ_space"] += curr_space

        # print(cntct_data)
        # exit()
        return cntct_data

    def get_cntct_NEW(self,option_dict,param_dict):
        """
        function to return the 2D list containing the cntct score features for all residues 
        listed in lig_resids. This is updated to define all space closest to a residue as a
        new kind of "pocket" associated with that residue. And output will include contact 
        as occupied space, but also "potential" contact by unoccupied space. Instead of
        "contact score" I can start by calling it occupied pocket score and unoccupied 
        pocket score OR contact alpha-space and unoccupied alpha-space
        (and also they are broken down by npol vs pol).
        """

        # empty data template for holding all cntct data by residue:

        # [overall cntct score, total space, total np space, %np, total occ space, % of total, match occ space, % of total, \
        # np match occ space, % of match, pol match occ space, % of match, mismatch occ space, % of total, np mismatch occ space,\
        # % of mismatch, pol mismatch occ space, % of mismatch, unocc space, % of total, np unocc space. % of unocc, pol unocc space,\
        # % of unocc, total envelope space, occ envelope space, match occ envelope space, mismatch occ envelope space]
        #
        # data_template = [0.0] * 27

        # initialize a dict to hold dicts for all cntct data
        cntct_data = {}
        # res_names = {}

        # for res in lig_resids:
        #     cntct_data[res] = {"space":0.0,"np_space":0.0,"occ_space":0.0,"match_occ_space":0.0,"np_match_occ_space":0.0,\
        #      "pol_match_occ_space":0.0,"mismatch_occ_space":0.0,"np_mismatch_occ_space":0.0,"pol_mismatch_occ_space":0.0,\
        #      "unocc_space":0.0,"np_unocc_space":0.0,"pol_unocc_space":0.0}

        data_template = {"space":                  0.0,"np_space": 0.0,"occ_space": 0.0,"match_occ_space": 0.0,
                         "np_match_occ_space":     0.0, \
                         "pol_match_occ_space":    0.0,"mismatch_occ_space": 0.0,"np_mismatch_occ_space": 0.0,
                         "pol_mismatch_occ_space": 0.0, \
                         "unocc_space":            0.0,"np_unocc_space": 0.0,"pol_unocc_space": 0.0}

        # make sure the "binder" has at least 5 atoms in order to run the reverse surface mapping
        if not self.lig_pdb:
            return cntct_data

        for i,p in enumerate(self.sub_pockets):
            for j,v in enumerate(p.vertices):
                curr_hits = []
                np_hit = False
                pol_hit = False
                low_dist = 1000.0
                low_resid = ''
                # low_resn = 'NA'
                for la_idx in self.face_lig_idx:
                    curr_dist = np.linalg.norm(
                        np.array(pdbinfo.coord(self.lig_pdb[la_idx])) - np.array(self.vertices[v]))
                    if curr_dist < low_dist:
                        low_dist = curr_dist
                        low_resid = self.lig_pdb[la_idx][23:26].strip()
                        # low_resid = pdbinfo.resi(self.lig_pdb[la_idx]).strip()
                        low_resn = pdbinfo.resn(self.lig_pdb[la_idx])
                    if curr_dist <= param_dict['hit_dist']:
                        curr_hits.append(la_idx)
                low_resid = str(low_resid).zfill(3)
                if low_resid not in cntct_data:
                    cntct_data[low_resid] = dict(data_template)
                    # cntct_data[low_resid] = {"space":0.0,"np_space":0.0,"occ_space":0.0,"match_occ_space":0.0,"np_match_occ_space":0.0\
                    #  "pol_match_occ_space":0.0,"mismatch_occ_space":0.0,"np_mismatch_occ_space":0.0,"pol_mismatch_occ_space":0.0,\
                    #   "unocc_space":0.0,"np_unocc_space":0.0,"pol_unocc_space":0.0}
                    # res_names[low_resid] = low_resn
                    cntct_data[low_resid]["resn"] = low_resn
                    cntct_data[low_resid]["pockets"] = []
                    cntct_data[low_resid]["AAC"] = []
                    cntct_data[low_resid]["AAC-spaces"] = []
                    cntct_data[low_resid]["pATOM"] = []

                if i + 1 not in cntct_data[low_resid]["pockets"]:
                    cntct_data[low_resid]["pockets"].append(i + 1)

                if low_resid:

                    # adding for pocket collapse onto residues
                    cntct_data[low_resid]["AAC"].append(v)
                    cntct_data[low_resid]["AAC-spaces"].append((p.aac_np_space[j],p.aac_pol_space[j],p.aac_space[j]))
                    for k in self.simplices[v]:
                        if k not in cntct_data[low_resid]["pATOM"]:
                            cntct_data[low_resid]["pATOM"].append(k)

                    cntct_data[low_resid]["space"] += p.aac_space[j]
                    cntct_data[low_resid]["np_space"] += p.aac_np_space[j]

                    if low_dist <= param_dict['hit_dist']:

                        cntct_data[low_resid]["occ_space"] += p.aac_space[j]
                        for la_hit in curr_hits:
                            if pdbinfo.atmn(self.lig_pdb[la_hit]).strip().startswith('N') or \
                                    pdbinfo.atmn(self.lig_pdb[la_hit]).strip().startswith('O') or \
                                            pdbinfo.atmn(self.prot_pdb[j]).strip() == 'SG':
                                pol_hit = True
                            else:
                                np_hit = True
                        if np_hit:
                            cntct_data[low_resid]["match_occ_space"] += p.aac_np_space[j]
                            cntct_data[low_resid]["np_match_occ_space"] += p.aac_np_space[j]
                        else:
                            cntct_data[low_resid]["mismatch_occ_space"] += p.aac_np_space[j]
                            cntct_data[low_resid]["np_mismatch_occ_space"] += p.aac_np_space[j]
                        if pol_hit:
                            cntct_data[low_resid]["match_occ_space"] += p.aac_pol_space[j]
                            cntct_data[low_resid]["pol_match_occ_space"] += p.aac_pol_space[j]
                        else:
                            cntct_data[low_resid]["mismatch_occ_space"] += p.aac_pol_space[j]
                            cntct_data[low_resid]["pol_mismatch_occ_space"] += p.aac_pol_space[j]
                    else:
                        cntct_data[low_resid]["unocc_space"] += p.aac_space[j]
                        cntct_data[low_resid]["np_unocc_space"] += p.aac_np_space[j]
                        cntct_data[low_resid]["pol_unocc_space"] += p.aac_pol_space[j]

        # #~~~~~NOTE: STILL NEED TO ADD ENVELOPE PART!!
        # #cntct_data_sorted = sorted(cntct_data.items(), key=operator.itemgetter(0))
        # #cntct_output = []


        # # for i in cntct_data_sorted:
        # #     curr_output = [str(i[0]).strip() + res_names[i[0]].strip(),str(i[1][6]+i[1][26]),str(i[1][1]),str(i[1][2]),str(i[1][2]/i[1][1]),\
        # #      str(i[1][4]),str(i[1][4]/i[1][1]),str(i[1][6]),str(i[1][6]/i[1][1]),str(i[1][8]),str(i[1][8]/[1][6]),str(i[1][10]),str(i[1][10]/i[1][6]),\
        # #      str(i[1][12]),str(i[1][12]/i[1][1]),str(i[1][14]),str(i[1][14]/i[1][12]),str(i[1][16]),str(i[1][16]/i[1][12]),\
        # #      str(i[1][18]),str(i[1][18]/i[1][1]),str(i[1][20]),str(i[1][20]/i[1][18]),str(i[1][22]),str(i[1][22]/i[1][18])]


        #     for i in cntct_data_sorted: 

        #         if i[1][6] == 0.0:
        #             curr_output = [str(i[0]).strip() + res_names[i[0]].strip(),i[1][6]+i[1][26],i[1][1],i[1][2],i[1][2]/i[1][1],\
        #              i[1][4],i[1][4]/i[1][1],i[1][6],i[1][6]/i[1][1],i[1][8],i[1][8]/i[1][6],i[1][10],i[1][10]/i[1][6],\
        #              i[1][12],i[1][12]/i[1][1],i[1][14],i[1][14]/i[1][12],i[1][16],i[1][16]/i[1][12],\
        #              i[1][18],i[1][18]/i[1][1],i[1][20],i[1][20]/i[1][18],i[1][22],i[1][22]/i[1][18]]


        #         else:

        #             curr_output = [str(i[0]).strip() + res_names[i[0]].strip(),i[1][6]+i[1][26],i[1][1],i[1][2],i[1][2]/i[1][1],\
        #              i[1][4],i[1][4]/i[1][1],i[1][6],i[1][6]/i[1][1],i[1][8],i[1][8]/i[1][6],i[1][10],i[1][10]/i[1][6],\
        #              i[1][12],i[1][12]/i[1][1],i[1][14],i[1][14]/i[1][12],i[1][16],i[1][16]/i[1][12],\
        #              i[1][18],i[1][18]/i[1][1],i[1][20],i[1][20]/i[1][18],i[1][22],i[1][22]/i[1][18]]


        #         # curr_output = [str(i[0]).strip() + str(res_names[i[0]].strip(),str(i[1][6]+i[26]),str(i[1][1]),str(i[1][2]),str(i[1][2]/i[1][1]),\
        #         #  str(i[4]),str(i[4]/i[1]),str(i[6]),str(i[6]/i[1]),str(i[8]),str(i[8]/i[6]),str(i[10]),str(i[10]/i[6]),\
        #         #  str(i[12]),str(i[12]/i[1]),str(i[14]),str(i[14]/i[12]),str(i[16]),str(i[16]/i[12]),\
        #         #  str(i[18]),str(i[18]/i[1]),str(i[20]),str(i[20]/i[18]),str(i[22]),str(i[22]/i[18])]

        #         cntct_output.append(curr_output)
        #     #now get the total face data
        #     tmp_matrix = [i[1] for i in cntct_data_sorted]

        #     i = []
        #     for j in range(len[tmp_matrix[0]]):

        #         i.append(sum([j[i] for j in tmp_matrix]))
        #     curr_output = ["FACE",i[6]+i[26],i[1],i[2],i[2]/i[1],\
        #      i[4],i[4]/i[1],i[6],i[6]/i[1],i[8],i[8]/i[6],i[10],i[10]/i[6],\
        #      i[12],i[12]/i[1],i[14],i[14]/i[12],i[16],i[16]/i[12],\
        #      i[18],i[18]/i[1],i[20],i[20]/i[18],i[22],i[22]/i[18]]

        # # curr_output = ["FACE",str(i[6]+i[26]),str(i[1]),str(i[2]),str(i[2]/i[1]),\
        # #  str(i[4]),str(i[4]/i[1]),str(i[6]),str(i[6]/i[1]),str(i[8]),str(i[8]/i[6]),str(i[10]),str(i[10]/i[6]),\
        # #  str(i[12]),str(i[12]/i[1]),str(i[14]),str(i[14]/i[12]),str(i[16]),str(i[16]/i[12]),\
        # #  str(i[18]),str(i[18]/i[1]),str(i[20]),str(i[20]/i[18]),str(i[22]),str(i[22]/i[18])]

        # #~~~~~~~NOTE: need to count the # of distinct pockets associated with each residues




        # cntct_output.append(curr_output)

        # returning dict with cntct data by resID

        return cntct_data

    def get_cntct(self,pdb_lines,lig_resids,option_dict,param_dict):
        """
        function to return the 2D list containing the cntct score features for all residues 
        listed in lig_resids.
        """
        cntct_arr = []
        counted_v = []
        cntct_p_list = []
        c_space_total = 0.0
        c_score_total = 0.0

        counted_lig_v = []
        e_space_total = 0.0
        e_score_total = 0.0

        p_cntct_list = []
        lig_p_cntct_list = []

        skip_reverse = False

        # make sure the "binder" has at least 5 atoms in order to run the reverse surface mapping
        if len(self.lig_pdb) >= 5:

            # if already did "reverse" mapping, then do regular for inverse mapping
            # if did regular, now do reverse for inverse mapping
            if option_dict['do_reverse']:
                envelope_do_reverse = False
            else:
                envelope_do_reverse = True

            pdb_lig = readPDB(pdb_lines,isTER=option_dict['use_TER'],isReverse=envelope_do_reverse, \
                              useNaccess=option_dict['use_naccess'],getFace=option_dict['get_face'])
            env_face_idx = self.face_idx

            # perform the Voronoi tessellation
            cutvert_lig,cutsimp_lig,cutrad_lig = AS_fcn.genVoro(pdb_lig,min_r=param_dict['min_r'],
                                                                max_r=param_dict['max_r'])

            if len(cutvert_lig) >= 1:

                # log the start of mapping
                AS_fcn.log("\nmapping the interface...\n",option_dict,to_screen=True,to_file=False)

                # perform fragmnet-centric topographical mapping (FCTM)
                ss_lig = AS_fcn.genPocket_wVoro(pdb_lig,cutvert_lig,cutsimp_lig,cutrad_lig,
                                                is_use_naccess=option_dict['use_naccess'], \
                                                is_get_face=option_dict['get_face'],clust_dist=param_dict['clust_dist'], \
                                                min_num_sph=param_dict['min_num_alph'],
                                                screen_by_face_atoms=option_dict['screen_by_face'], \
                                                screen_face_perc=option_dict['screen_face_perc'],
                                                hit_dist=param_dict['hit_dist'], \
                                                screen_by_score=option_dict['screen_by_score'],
                                                min_score=option_dict['min_score'], \
                                                screen_out_subsurf=option_dict['screen_out_subsurf'],
                                                subsurf_perc=option_dict['max_desolv_perc'], \
                                                is_use_exact_ACSA=option_dict['use_exact_ACSA'],
                                                screen_by_lig_contact=True, \
                                                expand_around_cntct=False,is_use_NP_wt=option_dict['use_NP_wt'],
                                                clust_method='average', \
                                                screen_by_resid=False,resid_list=[],screen_by_seam=False,
                                                expand_around_seam=False)
            else:
                skip_reverse = True
                print("\nSkipping mapping of the 'binder' surface due to too few binder voronoi vertices_coords...\n")

        else:
            skip_reverse = True
            print("\nSkipping mapping of the 'binder' surface due to too few binder atoms...\n")

        for i,res in enumerate(lig_resids):
            # get the trajectory for this residue
            # filter out polar atoms from the binder-contact score
            np_res_pdb = [line for line in self.lig_pdb if not pdbinfo.atmn(line).strip().startswith('N') and \
                          not pdbinfo.atmn(line).strip().startswith('O') and \
                          not pdbinfo.atmn(line).strip() == 'SG']

            coord_list = [pdbinfo.coord(line) for line in np_res_pdb if line[17:26].strip() == res]
            coord_list = np.array(coord_list)

            p_count = 0
            c_score = 0.0
            c_space = 0.0

            # cntct_atom_list_all = []

            for i2,p in enumerate(self.sub_pockets):
                p_found = False
                for v in p.vertices:
                    for coord in coord_list:
                        if np.linalg.norm(coord - self.vertices[v]) <= param_dict['hit_dist']:

                            # cntct_atom_list_all = cntct_atom_list_all + list(ss.simplex_idx_list[v])

                            # keep track of all cntct pockets for this snapshot
                            # if p not in cntct_p_list:
                            #    cntct_p_list.append(p)

                            # keep count of how many pockets are contacted by this residue
                            if not p_found:
                                p_found = True
                                p_count += 1

                            if i2 not in p_cntct_list:
                                p_cntct_list.append(i2)

                            # calculate the percent np SA for weighting of the alpha-space
                            np_sa = 0.0
                            tot_sa = 0.0
                            for j in self.simplices[v]:
                                tot_sa += p.asa_prot[j]
                                if pdbinfo.atmn(self.prot_pdb[j]).strip().startswith('N') or \
                                        pdbinfo.atmn(self.prot_pdb[j]).strip().startswith('O') or \
                                                pdbinfo.atmn(self.prot_pdb[j]).strip() == 'SG':
                                    continue
                                else:
                                    np_sa += p.asa_prot[j]
                            if tot_sa > 0.0:
                                pnpsa = np_sa / tot_sa
                            else:
                                pnpsa = 0.0

                            # calculate the alpha-space and np-weighted alpha-space
                            points = []
                            for j in self.simplices[v]:
                                points.append(pdbinfo.coord(self.prot_pdb[j]))
                            points = np.array(points)
                            vol = AS_fcn.calc_simp_vol(points)
                            np_vol = vol * pnpsa

                            c_space += vol
                            c_score += np_vol

                            # keep track of total interface space and score
                            if v not in counted_v:
                                #     c_space_total += vol
                                #     c_score_total += np_vol
                                counted_v.append(v)

                            # putting these outside the loop now, so that total pocket-binding cntct score
                            # will be exactly the sum of all residues' cntct scores (this is more interpretable)
                            c_space_total += vol
                            c_score_total += np_vol

                            break

            # NEW: first get the indices for the atoms in the current residue
            res_atm_idx_list = [idx for idx,line in enumerate(self.lig_pdb) if line[17:26].strip() == res]

            # then figure out which alpha-atoms are defined by these residue atoms and in cntct with the receptor (face)

            e_score = 0.0
            e_space = 0.0

            lig_p_count = 0

            if not skip_reverse:

                # for i2,aac in enumerate(cutvert_lig):
                for i2,lig_p in enumerate(ss_lig.sub_pockets):
                    lig_p_found = False

                    for lig_v in lig_p.vertices:

                        # found_lig_v = False
                        for atm_idx in res_atm_idx_list:

                            # this to use in case we want to screen out polar atoms from contributing to envelope contact score
                            # if pdbinfo.atmn(self.lig_pdb[atm_idx]).strip().startswith('N') or pdbinfo.atmn(self.lig_pdb[atm_idx]).strip().startswith('O') \
                            # or pdbinfo.atmn(self.lig_pdb[atm_idx]).strip().startswith('SG'):
                            #    continue



                            # pol_res_pdb = [line for line in res_pdb if pdbinfo.atmn(line).strip().startswith('N') or \
                            # pdbinfo.atmn(line).strip().startswith('O') or \
                            # pdbinfo.atmn(line).strip() == 'SG']

                            # if atm_idx in cutsimp_lig[i2]:
                            if atm_idx in cutsimp_lig[lig_v]:
                                for env_face_i in env_face_idx:
                                    # check if this vertex, which has a non-polar binder atom in its simp, is "occupied" by any receptor atoms
                                    # if np.linalg.norm(pdbinfo.coord(self.prot_pdb[env_face_i]) - aac) <= param_dict['hit_dist']:
                                    if np.linalg.norm(pdbinfo.coord(self.prot_pdb[env_face_i]) - ss_lig.vertices[
                                        lig_v]) <= param_dict['hit_dist']:

                                        # keep count of how many pockets are contacted by this residue
                                        if not lig_p_found:
                                            lig_p_found = True
                                            lig_p_count += 1

                                        if i2 not in lig_p_cntct_list:
                                            lig_p_cntct_list.append(i2)

                                        # TEST
                                        # curr_sa = pdb_lig.face_asa[atm_idx]
                                        curr_sa = lig_p.asa_prot[atm_idx]

                                        # np_sa = 0.0
                                        tot_sa = 0.0
                                        for j in cutsimp_lig[lig_v]:
                                            # tot_sa += pdb_lig.face_asa[j]
                                            tot_sa += lig_p.asa_prot[j]
                                            # if pdbinfo.atmn(pdb_lig.prot[j]).strip().startswith('N') or \
                                            #  pdbinfo.atmn(pdb_lig.prot[j]).strip().startswith('O') or \
                                            #  pdbinfo.atmn(pdb_lig.prot[j]).strip() == 'SG':
                                            #    continue
                                            # else:
                                            #    # np_sa += pdb_lig.face_asa[j]
                                            #    np_sa += lig_p.asa_prot[j]
                                        # if tot_sa > 0.0:
                                        #    pnpsa = np_sa / tot_sa
                                        # else:
                                        #    pnpsa = 0.0

                                        # TEST
                                        curr_perc = curr_sa / tot_sa

                                        points = []
                                        for j in cutsimp_lig[lig_v]:
                                            points.append(pdbinfo.coord(pdb_lig.prot[j]))
                                        points = np.array(points)
                                        vol = AS_fcn.calc_simp_vol(points)
                                        # np_vol = vol * pnpsa

                                        # e_space += vol
                                        # e_score += np_vol

                                        # TEST
                                        if not pdbinfo.atmn(self.lig_pdb[atm_idx]).strip().startswith(
                                                'N') and not pdbinfo.atmn(self.lig_pdb[atm_idx]).strip().startswith('O') \
                                                and not pdbinfo.atmn(self.lig_pdb[atm_idx]).strip().startswith('SG'):
                                            e_score += vol * curr_perc
                                            e_score_total += vol * curr_perc
                                        e_space += vol * curr_perc
                                        e_space_total += vol * curr_perc

                                        # found_lig_v = True

                                        # commenting this for TEST
                                        # if i2 not in counted_lig_v:
                                        # e_space_total += vol
                                        # #e_score_total += np_vol

                                        # #TEST
                                        # e_score_total += vol * curr_perc

                                        # counted_lig_v.append(i2)

                                        break

                                        # if found_lig_v:
                                        #     break

            resID = res.split()[0] + res.split()[-1]

            if c_space and e_space:

                features = [resID,c_score + e_score,p_count,c_score,c_score / c_space,c_space,lig_p_count,e_score,
                            e_score / e_space,e_space]
                # print(*features,sep='\t',file=out_files[i]) 
                # for i3,feat in enumerate(features[1:]):
                #     #print(feat)
                #     append_list[i][header[i3+1]].append(feat)

            elif c_space:

                features = [resID,c_score + e_score,p_count,c_score,c_score / c_space,c_space,lig_p_count,e_score,0.0,
                            e_space]
                # print(*features,sep='\t',file=out_files[i]) 
                # for i3,feat in enumerate(features[1:]):
                #     #print(feat)
                #     append_list[i][header[i3+1]].append(feat)

            elif e_space:

                features = [resID,c_score + e_score,p_count,c_score,0.0,c_space,lig_p_count,e_score,e_score / e_space,
                            e_space]
                # print(*features,sep='\t',file=out_files[i]) 
                # for i3,feat in enumerate(features[1:]):
                #     #print(feat)
                #     append_list[i][header[i3+1]].append(feat)

            else:
                features = [resID,c_score + e_score,p_count,c_score,0.0,c_space,lig_p_count,e_score,0.0,e_space]

            cntct_arr.append(features)

        if c_space_total:
            c_pnp_total = c_score_total / c_space_total
        else:
            c_pnp_total = 0.0
        if e_space_total:
            e_pnp_total = e_score_total / e_space_total
        else:
            e_pnp_total = 0.0

        face_features = ["face",c_score_total + e_score_total,len(p_cntct_list),c_score_total,c_pnp_total,c_space_total, \
                         len(lig_p_cntct_list),e_score_total,e_pnp_total,e_space_total]
        cntct_arr.append(face_features)

        self.cntct_arr = cntct_arr

        return cntct_arr

    def pocket_asa(self,idx):
        """
        Function to use asa list and pocket index to retreive the surface area of the pocket
        """
        return self.pockets[idx].asa(self.p_asa_prot)

    def write_resid_pockets_2(self,cntct_dict,core_cutoff=100.0,aux_cutoff=30.0,dir_name="."):

        # remove any existing pockets directory
        if os.path.isdir(dir_name + '/resid_pockets_2'):
            os.system("rm -r " + dir_name + "/resid_pockets_2")
        os.system("mkdir " + dir_name + "/resid_pockets_2")

        sorted_cntct_dict = sorted(iter(cntct_dict.items()),key=lambda x_y: x_y[1]['np_space'],reverse=True)

        for idx,pocket in enumerate(sorted_cntct_dict):
            f = open(dir_name + '/resid_pockets_2/' + str(idx).zfill(3) + '.pdb','w')

            ####Now I need to get the ACC position... and write this out
            ####just make a atom_xyz list as I go through and then find the average of that
            ####trajectory = []
            # for i in self.vertices_coords:
            #    trajectory.append(self.ss.vertices_coords[i])
            # center = (np.sum(trajectory,axis=0)) / len(self.vertices_coords)

            for a_idx in pocket[1]["pATOM"]:
                f.write(self.prot_pdb[a_idx][0:21])
                f.write('X')
                # f.write(self.prot_pdb[a_idx][22:60])
                f.write(self.prot_pdb[a_idx][22:54])

                f.write('\n')

            coords = []

            for i,v_idx in enumerate(pocket[1]["AAC"]):

                coords.append(self.vertices[v_idx])

                f.write('ATOM   ' + str(i).rjust(4) + '  AAO AAC ')
                f.write('X')
                f.write('   0    ')

                for j in self.vertices[v_idx]:
                    print('{0:8.3f}'.format(j),end="",file=f)
                for j in pocket[1]["AAC-spaces"][i]:
                    print('{0:6.1f}'.format(j),end="",file=f)
                print('{0:6.2f}\n'.format(self.radii[v_idx]),end="",file=f)

                # center = (np.sum(trajectory,axis=0)) / len(pocket[1]["AAC"])
            center = (np.average(coords,axis=0))

            f.write('ATOM   ' + str(1).rjust(4) + '  ACC ACC ')
            f.write('X')
            f.write('   0    ')

            for j in center:
                print('{0:8.3f}'.format(j),end="",file=f)
            print('\n',end="",file=f)

            f.close()

        # for key in cntct_dict:
        #     f = open(dir_name + '/resid_pockets/' + key + '.pdb','w')

        #     for a_idx in cntct_dict[key]["pATOM"]:
        #         f.write(self.prot_pdb[a_idx][0:21])
        #         f.write('X')
        #         #f.write(self.prot_pdb[a_idx][22:60])
        #         f.write(self.prot_pdb[a_idx][22:54])

        #         f.write('\n')

        #     for i,v_idx in enumerate(cntct_dict[key]["AAC"]):
        #         f.write('ATOM   ' + str(i).rjust(4) + '  AAO AAC ')
        #         f.write('X')
        #         f.write('   0    ')

        #         for j in self.vertices_coords[v_idx]:
        #             print('{0:8.3f}'.format(j),end="",file=f)
        #         print('{0:12.2f}\n'.format(self.radii[v_idx]),end="",file=f)  

        #     f.close()

        return

    def write_resid_pockets(self,cntct_dict,core_cutoff=100.0,aux_cutoff=30.0,dir_name="."):

        # remove any existing pockets directory
        if os.path.isdir(dir_name + '/resid_pockets'):
            os.system("rm -r " + dir_name + "/resid_pockets")
        os.system("mkdir " + dir_name + "/resid_pockets")

        for key in cntct_dict:
            f = open(dir_name + '/resid_pockets/' + key + '.pdb','w')

            for a_idx in cntct_dict[key]["pATOM"]:
                f.write(self.prot_pdb[a_idx][0:21])
                f.write('X')
                # f.write(self.prot_pdb[a_idx][22:60])
                f.write(self.prot_pdb[a_idx][22:54])

                f.write('\n')

            for i,v_idx in enumerate(cntct_dict[key]["AAC"]):
                f.write('ATOM   ' + str(i).rjust(4) + '  AAO AAC ')
                f.write('X')
                f.write('   0    ')

                for j in self.vertices[v_idx]:
                    print('{0:8.3f}'.format(j),end="",file=f)
                print('{0:12.2f}\n'.format(self.radii[v_idx]),end="",file=f)

            f.close()

        return

    def write_pockets(self,core_cutoff=100.0,aux_cutoff=30.0,dir_name="."):
        """
        Function: to write out PDB files of all pockets for visualization
        Arguments: high: threshold value to qualify as high scoring
                   mid: threshold value to qualify as medium scoring
                   dir_name: directory to write output
        """
        # remove any existing pockets directory
        if os.path.isdir(dir_name + '/pockets'):
            os.system("rm -r " + dir_name + "/pockets")
        os.system("mkdir " + dir_name + "/pockets")

        for idx,p in enumerate(self.sub_pockets):
            f = open(dir_name + '/pockets/' + str(idx).zfill(3) + '.pdb','w')
            if p.score >= core_cutoff:
                chain = 'H'
            elif p.score >= aux_cutoff:
                chain = 'M'
            else:
                chain = 'L'
            for a_idx in p.atoms:
                f.write(self.prot_pdb[a_idx][0:21])
                f.write(chain)
                # f.write(self.prot_pdb[a_idx][22:60])
                f.write(self.prot_pdb[a_idx][22:54])

                # write atomic space
                f.write('{0:6.2f}'.format(p.atomic_space[a_idx]))

                # write atomic ACSA
                f.write('{0:6.2f}'.format(p.asa_prot[a_idx]) + '\n')
            for i,v_idx in enumerate(p.vertices):
                if p.vert_occ_stat[i] == True:
                    f.write('ATOM   ' + str(i).rjust(4) + '  AAO AAC ')
                    f.write(chain)
                    f.write('   0    ')


                else:
                    f.write('ATOM   ' + str(i).rjust(4) + '  AAU AAC ')
                    f.write(chain)
                    f.write('   0    ')

                for j in self.vertices[v_idx]:
                    print('{0:8.3f}'.format(j),end="",file=f)
                print('{0:12.2f}\n'.format(self.radii[v_idx]),end="",file=f)

            vert_centroid = p.vert_cntr()
            atm_centroid = p.atm_cntr()

            f.write('ATOM    100  PAC PAC ')
            f.write(chain)
            f.write('   0    ')
            for i in atm_centroid:
                print('{0:8.3f}'.format(i),end="",file=f)
            f.write('\n')
            f.write('ATOM    101  ACC ACC ')
            f.write(chain)
            f.write('   0    ')
            for i in vert_centroid:
                print('{0:8.3f}'.format(i),end="",file=f)
            f.close()

    def set_p_asa_prot(self):
        """
        Function: Set the snapshot level atomistic ACSA (alpha-atom contact SA)
                  to be used for non-polar weighting if NOT using exact (pocket-level) ACSA 
        Method: Write to a file the prot_pdb and then all alpha-atoms of all pockets in the snapshot.
                Calculate the SASA of this receptor-multiple pocket complex and take the difference
                from the apo receptor to get the total ACSA at the atomistic level.
        """
        holo_pdb = self.prot_pdb[:]
        for p in self.sub_pockets:
            for i,v_idx in enumerate(p.vertices):
                ln = 'ATOM   1000  AAC AAC X'
                ln = ln + str(i).rjust(4) + '    '

                for j in self.vertices[v_idx]:
                    ln = ln + '{0:8.3f}'.format(j)
                ln = ln + '{0:12.2f}'.format(self.radii[v_idx])
                holo_pdb.append(ln)

        asa_holo = AS_fcn.calc_asa(holo_pdb)[0:len(self.prot_pdb)]

        # make into np.arrays and take the difference to get the actual pocket SA
        numpy_asa = np.array(self.asa_prot)
        numpy_asa_holo = np.array(asa_holo)
        numpy_p_asa = np.subtract(numpy_asa,numpy_asa_holo)

        self.p_asa_prot = numpy_p_asa

        return

    def get_clusters_percSA(self):
        """
        Function: evaluate the occluded SA and npSA (SA weighted by npsa of interacting pocket) 
                  for each alpha-cluster a return in an array and also the percent of SA occluded
        To do still: as is, this method is too expensive. it recalculates the receptor-cluster
                     SASA for every cluster. we should calculate this once, and then be smarter
                     about finding the alpha-spheres from this list of all alpha-sphere SASA.
                     need to know the exact number of pockets and the length for each such that
                     we can move through them and get the individual cluster SASAs.
        """
        clusters_percSA = []
        clusters_npSA = []
        for p in self.sub_pockets:
            cluster_pdb = []
            for i,v_idx in enumerate(p.vertices):
                ln = 'ATOM   1000  AAC AAC X'
                ln = ln + str(i).rjust(4) + '    '
                for j in self.vertices[v_idx]:
                    ln = ln + '{0:8.3f}'.format(j)
                ln = ln + '{0:12.2f}'.format(self.radii[v_idx])
                cluster_pdb.append(ln)
            # get the pdb for the receptor with a single pocket
            cluster_prot_pdb = cluster_pdb[:]
            cluster_prot_pdb.extend(self.prot_pdb)
            cluster_gas_asa = np.array(AS_fcn.calc_asa(cluster_pdb))
            cluster_complex_asa = np.array(AS_fcn.calc_asa(cluster_prot_pdb)[0:len(cluster_gas_asa)])

            # cluster_covered_asa = np.subtract(cluster_gas_asa,cluster_complex_asa)
            amount_gas = np.sum(cluster_gas_asa)
            amount_cov = amount_gas - np.sum(cluster_complex_asa)
            perc_cov = amount_cov / amount_gas

            clusters_percSA.append(perc_cov)
        return (clusters_percSA)

    def write_pdb(self,dir_name='.',clust_dist='4.0'):
        """
        Function: write out the pdb files for the receptor, the binder, and the interface atoms
        """

        if os.path.isdir(dir_name + '/pdb_out'):
            os.system("rm -r " + dir_name + "/pdb_out")
        os.system("mkdir " + dir_name + "/pdb_out")

        f = open(dir_name + '/pdb_out/prot.pdb','w')
        for line in self.prot_pdb:
            f.write("%s\n" % line)
        f.close()

        f = open(dir_name + '/pdb_out/lig.pdb','w')
        f2 = open(dir_name + '/pdb_out/lig_LFC.pdb','w')
        for line in self.lig_pdb:
            f.write("%s\n" % line)
            f2.write("%s\n" % line)
        # to get/write LFCs
        f.close()

        f2.write("TER\n")

        if self.lig_pdb:

            LFCs = AS_fcn.get_LFC_coord(self.lig_pdb,clust_dist)

            for count,i in enumerate(LFCs):
                f2.write('ATOM    100  LFC LFC X' + str(count).rjust(4) + '    ')
                for j in i:
                    print('{0:8.3f}'.format(j),end="",file=f2)
                f2.write('\n')

            f2.close()

        # lig_coords = []
        # for line in self.lig_pdb:
        #     lig_coords.append(pdbinfo.coord(line))
        # #zmat = hier.linkage(lig_coords, method='average')
        # zmat = hier.linkage(lig_coords, method='complete')
        # #cluster = hier.fcluster(zmat,clust_dist,criterion='distance')
        # cluster = hier.fcluster(zmat,3.5,criterion='distance')

        # clusters = [ [] for i in range(max(cluster)) ]
        # AS_fcn.fill_p_verts(clusters,cluster)

        # for count,i in enumerate(clusters):
        #     frag_coords = []
        #     for idx in i:
        #         frag_coords.append(lig_coords[idx])
        #     frag_centroid = np.average(frag_coords,axis=0)

        #     f.write('ATOM    100  LFC LFC X' + str(count).rjust(4) + '    ')
        #     for j in frag_centroid:
        #         print('{0:8.3f}'.format(j),end="",file=f)
        #     f.write('\n')


        f = open(dir_name + '/pdb_out/face.pdb','w')
        if hasattr(self,'face_idx'):
            for i in self.face_idx:
                f.write("%s\n" % self.prot_pdb[i])
        # if no "face" defined, use the list of all pocket atoms
        elif self.sub_pockets:
            all_p_atoms = []
            for p in self.sub_pockets:
                for a in p.atoms:
                    if a not in all_p_atoms:
                        f.write("%s\n" % self.prot_pdb[a])
                    all_p_atoms.append(a)
        f.close()

        return

    def write_tetras(self,dir_name='.'):
        """
        Function: write out the file for visualizing the alpha-space tetrahedrons with Chimera
        """
        if os.path.isdir(dir_name + '/draw'):
            os.system("rm -r " + dir_name + "/draw")
        os.system("mkdir " + dir_name + "/draw")
        colors = [lines.rstrip() for lines in open('colors_chimera.txt','r')]
        colors = colors * 20

        f = open(dir_name + '/draw/alpha_spaces.bild','w')

        for idx,p in enumerate(self.sub_pockets):
            f.write(".color " + colors[idx] + " \n")
            for v in p.vertices:
                fourpoint = self.simplices[v]
                for i in range(3):
                    for j in range(i + 1,4,1):
                        f.write(".cylinder ")
                        for k in pdbinfo.coord(self.prot_pdb[fourpoint[i]]):
                            f.write(str(k) + " ")
                        for k in pdbinfo.coord(self.prot_pdb[fourpoint[j]]):
                            f.write(str(k) + " ")
                        f.write(" 0.1 \n")
        f.close()
        return

    def write_score_table(self,dir_name='.',write_header=True):
        """
        Function: write out the basic "score table" (includes: pocket rank, color,
        score, space, percent occupied, and percent nonpolar)
        """
        f = open(dir_name + '/table_pockets.dat','w')

        header = ['#rank','color','score','% occupied','space','% non-polar']
        if write_header:
            print(*header,sep='\t',file=f)

        colors = [lines.rstrip() for lines in open('colors_table.txt','r')]
        colors = colors * 20

        for i,p in enumerate(self.sub_pockets):
            rank = i + 1
            color = colors[i]
            score = "{0}".format(int(round(p.score,0)))
            perc_occ = "{0}".format(int(round(p.perc_occ,2) * 100.0)) + '%'
            space = "{0}".format(int(round(p.space,0)))

            if p.pnp == 'n/a':
                pnp = 'n/a'
            else:
                pnp = "{0}".format(int(round(p.pnp,2) * 100.0)) + '%'

            data = [rank,color,score,perc_occ,space,pnp]
            print(*data,sep='\t\t',file=f)

        return

    def get_total_vol(self,occ_only=False):
        """
        Function: return the total alpha-cluster volume for all pockets in "sub_pockets"
        Options: is "occ_only = True", filters out all alpha-atoms NOT in binder contact
        """
        vert = []
        if occ_only:
            for p in self.sub_pockets:
                vert = vert + p.vert_occ
        else:
            for p in self.sub_pockets:
                vert = vert + p.vertices

        p_coord = []
        p_rad = []
        for v in vert:
            p_coord.append(self.vertices[v])
            p_rad.append(1.8)
        p_rad = np.array(p_rad)
        p_coord = np.array(p_coord)
        if p_coord.size:
            # return AS_fcn.estVolume(p_coord,p_rad,5000)
            return AS_fcn.estVolume(p_coord,p_rad,20000)
        else:
            return 0.0


class Pocket:
    """
    Class: main pocket class to hold pocket atom indices, alpha-cluster indices, and 
           pocket features used for evaluation  
    """

    def __init__(self,vertices,atoms):

        self.ss = None
        self.ss_idx = None
        self.vertices = vertices
        self.atoms = atoms

        self.asa_prot = []

        self.score = 0.0
        self.space = 0.0
        self.pnp = 0.0
        self.occ_score = 0.0
        # occ_cntct_score counted for np binder only
        self.occ_cntct_score = 0.0
        self.occ_space = 0.0
        self.occ_pnp = 0.0
        self.unocc_score = 0.0
        # unocc_cntct_score includes polar binder cntct and empty alpha-spaces
        self.unocc_cntct_score = 0.0
        self.unocc_space = 0.0
        self.unocc_pnp = 0.0

        self.perc_occ = 0.0

        self.vert_occ = []
        self.vert_occ_stat = []

        # dict to store space by pocket atom (keys are idx in prot_atm array)
        self.atomic_space = {}

        # lists to save space, np space, and pol space by alpha atom (vertice)
        self.aac_space = []
        self.aac_np_space = []
        self.aac_pol_space = []

        # list of relative simplex atom SA percents for each vertice
        self.vert_atomSA_percents = []

    def set_ss(self,ss):
        """
        Function: to associate the residue Snapshot with its Pocket
        """
        self.ss = ss
        return

    def set_ss_idx(self,ss_idx):
        """
        set the residue Snapshot index
        :param ss_idx: int
        """
        self.ss_idx = ss_idx

    def vert_cntr(self):
        """
        Function: calculate and return the centroid of the alpha-cluster
        """
        coords = []
        for i in self.vertices:
            coords.append(self.ss.vertices[i])
        center = (np.sum(coords,axis=0)) / len(self.vertices)
        return center

    def atm_cntr(self):
        """
        Function: calculate and return the centroid of the atoms from the pocket
        """
        coords = []
        for atm in self.atoms:
            coords.append(pdbinfo.coord(self.ss.prot_pdb[atm]))
        center = (np.sum(coords,axis=0)) / len(self.atoms)
        return center

    def is_contact(self,hit_dist=1.6):
        """
        Function to determine if the pocket is a binder contact pocket
        """
        # get binder atoms trajectory
        lig_coord = [pdbinfo.coord(l) for l in self.ss.lig_pdb]
        lig_coord = np.array(lig_coord)
        status = False
        # in case later I want to add a "count" threshold
        atm_cnt = 0
        for v in self.vertices:
            for a in lig_coord:
                if np.linalg.norm(a - self.ss.vertices[v]) <= hit_dist:
                    atm_cnt += 1
                    status = True
                    break
            if status == True:
                break
        return status

    def get_exact_asa_prot(self):
        """
        Function: To calculate the "exact" SA for a single pocket
        """
        # import AS_fcn
        holo_pdb = self.ss.prot_pdb[:]

        for i,v_idx in enumerate(self.vertices):
            ln = 'ATOM   1000  AAC AAC X'
            ln = ln + str(i).rjust(4) + '    '

            for j in self.ss.vertices[v_idx]:
                ln = ln + '{0:8.3f}'.format(j)
            ln = ln + '{0:12.2f}'.format(self.ss.radii[v_idx])
            holo_pdb.append(ln)

        asa_holo = AS_fcn.calc_asa(holo_pdb)[0:len(self.ss.prot_pdb)]

        # make into np.arrays and take the difference to get the actual pocket SA
        numpy_asa = np.array(self.ss.asa_prot)
        numpy_asa_holo = np.array(asa_holo)
        numpy_p_asa = np.subtract(numpy_asa,numpy_asa_holo)

        return numpy_p_asa

    def get_cluster_vol(self):
        """
        Function: get the alpha-cluster volume 
        """
        vert = self.ss.vertices
        p_coord = []
        p_rad = []
        for v in self.vertices:
            p_coord.append(vert[v])
            p_rad.append(1.8)
        p_rad = np.array(p_rad)
        p_coord = np.array(p_coord)
        # return AS_fcn.estVolume(p_coord,p_rad,5000)
        return AS_fcn.estVolume(p_coord,p_rad,10000)

    def get_cntct_vol(self):
        """
        Function: get the occupied only alpha-cluster volume
        """
        vert = self.ss.vertices
        p_coord = []
        p_rad = []
        for v in self.vert_occ:
            p_coord.append(vert[v])
            p_rad.append(1.8)
        p_rad = np.array(p_rad)
        p_coord = np.array(p_coord)
        if p_coord.size:
            # return AS_fcn.estVolume(p_coord,p_rad,5000)
            return AS_fcn.estVolume(p_coord,p_rad,10000)
        else:
            return 0.0

    def set_space_score_occ(self,hit_dist):
        """
        Function: set pocket features (when using non-polar weighting): 
        score, alpha-space, percent occupied, percent non-polar
        """
        vert = self.ss.vertices
        simp = self.ss.simplices
        pdb = self.ss.prot_pdb
        lig_pdb = self.ss.lig_pdb
        asa = self.asa_prot

        # get binder atoms trajectory
        lig_coord = [pdbinfo.coord(l) for l in lig_pdb]
        lig_coord = np.array(lig_coord)
        lig_atmn = [pdbinfo.atmn(l).strip() for l in lig_pdb]

        # pocket features
        space = 0.0
        np_space = 0.0  # score
        occ_space = 0.0
        occ_np_space = 0.0
        occ_score = 0.0
        unocc_space = 0.0
        unocc_np_space = 0.0
        unocc_score = 0.0

        pnpsas = self.pnpsa_by_asph()

        vert_occ_stat = []
        vert_cntct = []

        # open dict to store alpha-space by pocket atom for feature set
        atomic_space = {}

        aac_space = []
        aac_np_space = []
        aac_pol_space = []

        # cycle through the alpha-spheres
        for i,v in enumerate(self.vertices):
            hit = False
            np_hit = False

            # ~~~NOTE: I could probably remove this and put it into the new contact score fcn, because now it's redundant!!~~~~
            # check for binder contact
            for i2,a in enumerate(lig_coord):
                if np.linalg.norm(a - vert[v]) <= hit_dist:
                    hit = True
                    # check if cntct binder atom is non-polar
                    if not lig_atmn[i2].startswith('N') and not lig_atmn[i2].startswith('O') \
                            and not lig_atmn[i2] == 'SG':
                        np_hit = True
                        break

            # assign the space and np_space for the vertex
            points = []
            for j in simp[v]:
                points.append(pdbinfo.coord(pdb[j]))
            points = np.array(points)
            vol = AS_fcn.calc_simp_vol(points)
            np_vol = vol * pnpsas[i]
            space += vol
            np_space += np_vol

            self.aac_space.append(vol)
            self.aac_np_space.append(np_vol)
            self.aac_pol_space.append(vol - np_vol)

            # new for feature set...
            for i2,j in enumerate(simp[v]):
                if j not in atomic_space:
                    atomic_space[j] = vol * self.vert_atomSA_percents[i][i2]
                    # atomic_space[j] = vol / 4.0
                else:
                    atomic_space[j] += vol * self.vert_atomSA_percents[i][i2]
                    # atomic_space[j] += vol / 4.0

            if hit:
                vert_cntct.append(v)
                vert_occ_stat.append(True)
                occ_space += vol
                occ_np_space += np_vol

            else:
                unocc_space += vol
                unocc_np_space += np_vol
                vert_occ_stat.append(False)

            if hit and np_hit:
                occ_score += np_vol
            else:
                unocc_score += np_vol

        self.space = space
        self.score = np_space
        self.occ_space = occ_space
        self.occ_score = occ_np_space
        self.occ_cntct_score = occ_score
        self.unocc_space = unocc_space
        self.unocc_score = unocc_np_space
        self.unocc_cntct_score = unocc_score

        self.perc_occ = (occ_space / space if space else 0.0)

        self.vert_occ_stat = vert_occ_stat
        self.vert_occ = vert_cntct

        self.atomic_space = atomic_space

        if space >= 0.5:
            self.pnp = round(np_space,0) / round(space,0)
        else:
            self.pnp = (round(np_space,0) / space if space else 0.0)

        if occ_space > 0.0:
            self.occ_pnp = occ_np_space / occ_space
        else:
            self.occ_pnp = 'n/a'
        if unocc_space > 0.0:
            self.unocc_pnp = unocc_np_space / unocc_space
        else:
            self.unocc_pnp = 'n/a'

        return

    def set_space_as_score_occ(self,hit_dist):
        """
        Function: set pocket features (if NOT using non-polar weighting):
        score (same as alpha-space), alpha-space, percent occupied, percent non-polar (is N/A)
        """
        vert = self.ss.vertices
        simp = self.ss.simplices
        pdb = self.ss.prot_pdb
        lig_pdb = self.ss.lig_pdb

        # get binder atoms trajectory
        lig_coord = [pdbinfo.coord(l) for l in lig_pdb]
        lig_coord = np.array(lig_coord)

        # pocket features
        space = 0.0
        occ_space = 0.0
        unocc_space = 0.0

        vert_occ_stat = []
        vert_cntct = []

        # cycle through the alpha-spheres
        for i,v in enumerate(self.vertices):
            hit = False

            # check for binder contact
            for a in lig_coord:
                if np.linalg.norm(a - vert[v]) <= hit_dist:
                    hit = True
                    break

            # assign the space and np_space for the vertex
            points = []
            for j in simp[v]:
                points.append(pdbinfo.coord(pdb[j]))
            points = np.array(points)
            vol = AS_fcn.calc_simp_vol(points)
            space += vol
            if hit:
                vert_cntct.append(v)
                occ_space += vol
                vert_occ_stat.append(True)
            else:
                unocc_space += vol
                vert_occ_stat.append(False)

        self.space = space
        self.score = space
        self.occ_space = occ_space
        self.occ_score = occ_space
        self.occ_cntct_score = occ_space
        self.unocc_space = unocc_space
        self.unocc_score = unocc_space
        self.unocc_cntct_score = unocc_space

        self.perc_occ = (occ_space / space if space else 0.0)

        self.vert_occ_stat = vert_occ_stat
        self.vert_occ = vert_cntct

        self.pnp = 'n/a'
        self.occ_pnp = 'n/a'
        self.unocc_pnp = 'n/a'

        return

    def pnpsa_by_asph(self):
        """
        Function: return a list of pnpsa values by alpha-sphere 
        """
        asph_pnpsa = []
        for v in self.vertices:
            np_sa = 0.0
            tot_sa = 0.0
            for i in self.ss.simplices[v]:
                tot_sa += self.asa_prot[i]
                if pdbinfo.atmn(self.ss.prot_pdb[i]).strip().startswith('N') or \
                        pdbinfo.atmn(self.ss.prot_pdb[i]).strip().startswith('O') or \
                                pdbinfo.atmn(self.ss.prot_pdb[i]).strip() == 'SG':
                    continue
                else:
                    np_sa += self.asa_prot[i]

            if tot_sa > 0.0:
                asph_pnpsa.append(np_sa / tot_sa)
            else:
                asph_pnpsa.append(0.0)

            # adding this list to pocket to store the relative percent SA of the 4 simplex atoms for each vertex
            curr_atomSA_percent = []
            for i in self.ss.simplices[v]:
                curr_atomSA_percent.append(self.asa_prot[i] / tot_sa)
            self.vert_atomSA_percents.append(curr_atomSA_percent)

        return asph_pnpsa

    def perc_clust_cov(self):
        """
        Function: return the percent of the alpha-cluster SA that is occluded by the receptor
        """
        cluster_pdb = []
        for i,v_idx in enumerate(self.ss.vertices):
            ln = 'ATOM   1000  AAC AAC X'
            ln = ln + str(i).rjust(4) + '    '
            for j in self.ss.vertices[v_idx]:
                ln = ln + '{0:8.3f}'.format(j)
            ln = ln + '{0:12.2f}'.format(self.ss.radii[v_idx])
            cluster_pdb.append(ln)
        # get the pdb for the receptor with a single pocket
        cluster_prot_pdb = cluster_pdb[:]
        cluster_prot_pdb.extend(self.ss.prot_pdb)
        cluster_gas_asa = np.array(AS_fcn.calc_asa(cluster_pdb))
        cluster_complex_asa = np.array(AS_fcn.calc_asa(cluster_prot_pdb)[0:len(cluster_gas_asa)])

        # cluster_covered_asa = np.subtract(cluster_gas_asa,cluster_complex_asa)
        amount_gas = np.sum(cluster_gas_asa)
        amount_cov = amount_gas - np.sum(cluster_complex_asa)
        perc_cov = amount_cov / amount_gas

        return perc_cov


class Beta_clust:
    """
    Class: beta-cluster class will hold list of Beta-atoms and associated lists for 
    vert, simp, and radii
    Funtion: will consolidate the information needed to look up vert, simp, and radii,
    for each distinct fragment-centric beta-cluster
    """

    def __init__(self,vert_beta,simp_beta,rad_beta):
        self.vert_beta = vert_beta
        self.simp_beta = simp_beta
        self.rad_beta = rad_beta
        self.beta_atoms = []


class Beta_atom:
    """
    Class: beta-atom class will hold list of alpha-atoms and contact pocket atoms
    Functions: to get the total space and score and the classification vector
    """

    def __init__(self,vertices,atoms,ss,b_clust):

        self.ss = ss
        self.b_clust = b_clust
        self.vertices = vertices
        self.atoms = atoms

        self.asa_prot = []

        self.occ_status = None

        self.space = None
        self.score = None
        self.pol_space = None

    def get_occ_status(self,hit_dist=1.6):

        vert = self.b_clust.vert_beta
        lig_pdb = self.ss.lig_pdb

        # get binder atoms trajectory
        lig_coord = [pdbinfo.coord(l) for l in lig_pdb]
        lig_coord = np.array(lig_coord)
        # lig_atmn = [pdbinfo.atmn(l).strip() for l in lig_pdb]

        beta_center = self.vert_cntr()
        hit = False
        for i in lig_coord:
            if np.linalg.norm(i - beta_center) <= hit_dist:
                hit = True
                break
        if hit:
            self.occ_status = True
        else:
            self.occ_status = False
        return hit

    def vert_cntr(self):
        """
        Function: calculate and return the centroid of the beta-atom
        """
        coords = []
        for i in self.vertices:
            coords.append(self.b_clust.vert_beta[i])
        center = (np.sum(coords,axis=0)) / len(self.vertices)
        return center

    def beta_score_space(self):
        """
        Function: set/return beta-score and beta-space features 
        score, alpha-space, percent occupied, percent non-polar
        """
        import numpy as np

        vert = self.b_clust.vert_beta
        simp = self.b_clust.simp_beta
        pdb = self.ss.prot_pdb
        lig_pdb = self.ss.lig_pdb
        asa = self.asa_prot

        # #get binder atoms trajectory
        # lig_coord = [pdbinfo.coord(l) for l in lig_pdb]
        # lig_coord = np.array(lig_coord)
        # lig_atmn = [pdbinfo.atmn(l).strip() for l in lig_pdb]

        # pocket features
        space = 0.0
        np_space = 0.0  # score
        pol_space = 0.0

        pnpsas = self.pnpsa_by_asph()

        # cycle through the alpha-spheres
        for i,v in enumerate(self.vertices):
            # hit = False
            # np_hit = False

            # #check for binder contact
            # for i2,a in enumerate(lig_coord):
            #     if np.linalg.norm(a - vert[v]) <= hit_dist:
            #         hit = True
            #         #check if cntct binder atom is non-polar
            #         if not lig_atmn[i2].startswith('N') and not lig_atmn[i2].startswith('O') \
            #          and not lig_atmn[i2] == 'SG':
            #             np_hit = True 
            #             break

            # assign the space and np_space for the vertex
            points = []
            for j in simp[v]:
                points.append(pdbinfo.coord(pdb[j]))
            points = np.array(points)
            vol = AS_fcn.calc_simp_vol(points)
            np_vol = vol * pnpsas[i]
            pol_vol = vol - np_vol
            space += vol
            np_space += np_vol
            pol_space += pol_vol
            # if hit:
            #     vert_cntct.append(v)
            #     vert_occ_stat.append(True)
            #     occ_space += vol
            #     occ_np_space += np_vol

            # else:
            #     unocc_space += vol
            #     unocc_np_space += np_vol
            #     vert_occ_stat.append(False)

            # if hit and np_hit:
            #         occ_score += np_vol
            # else:
            #     unocc_score += np_vol

        self.space = space
        self.score = np_space
        self.pol_space = pol_space
        # self.occ_space = occ_space
        # self.occ_score = occ_np_space
        # self.occ_cntct_score = occ_score
        # self.unocc_space = unocc_space
        # self.unocc_score = unocc_np_space
        # self.unocc_cntct_score = unocc_score

        # self.perc_occ = (occ_space / space if space else 0.0)

        # self.vert_occ_stat = vert_occ_stat
        # self.vert_occ = vert_cntct

        # if space >= 0.5:
        #     self.pnp = round(np_space,0) / round(space,0)
        # else:
        #     self.pnp = (round(np_space,0) / space if space else 0.0)

        # if occ_space > 0.0:
        #     self.occ_pnp = occ_np_space / occ_space
        # else:
        #     self.occ_pnp = 'n/a' 
        # if unocc_space > 0.0:
        #     self.unocc_pnp = unocc_np_space / unocc_space
        # else:
        #     self.unocc_pnp = 'n/a'

        return [np_space,pol_space,space]

    def pnpsa_by_asph(self):
        """
        Function: return a list of pnpsa values by alpha-sphere 
        """
        asph_pnpsa = []
        for v in self.vertices:
            np_sa = 0.0
            tot_sa = 0.0
            for i in self.b_clust.simp_beta[v]:
                tot_sa += self.ss.asa_prot[i]
                if pdbinfo.atmn(self.ss.prot_pdb[i]).strip().startswith('N') or \
                        pdbinfo.atmn(self.ss.prot_pdb[i]).strip().startswith('O') or \
                                pdbinfo.atmn(self.ss.prot_pdb[i]).strip() == 'SG':
                    continue
                else:
                    np_sa += self.ss.asa_prot[i]

            if tot_sa > 0.0:
                asph_pnpsa.append(np_sa / tot_sa)
            else:
                asph_pnpsa.append(0.0)
        return asph_pnpsa

# class Beta_cluster:
#     """
#     Class: beta-cluster class holds list of beta-atoms
#     Function: to write out the pdb file 
#     """
#     def __init__(self,beta_atoms):

#         self.ss = None
#         self.beta_atoms = beta_atoms

#         self.asa_prot = []
