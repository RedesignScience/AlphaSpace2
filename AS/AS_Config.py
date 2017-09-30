from colour import Color


class AS_Config(object):
    def __init__(self):
        self.min_r = 3.2
        self.max_r = 5.4

        self.pocket_cluster_distance = 4.7

        self.contact_threshold = 3.0

        self.screen_by_ligand_contact = True


        self.screen_by_face = True
        ### select pockets that include atoms from the interface atom list

        self.screen_face_perc = 0.5
        ### minimum percentage of atoms that must be in the interface atom list



        self.screen_by_lig_cntct = False
        ### select only pockets in "contact" with the ligand
        expand_around_cntct = False
        ### expand selection to include all pockets overlapping the "contact" pockets

        self.screen_by_seam = False
        ### select pocket that contain atoms from >1 protein chain
        expand_around_seam = False
        ### exapnd selection to include all pockets overlapping the "seam" pockets


        self.screen_by_score = False
        self.min_score = 100
        ### select only pockets above a minimum score



        self.screen_by_perc_rank = False
        self.min_perc_rank = 0.9
        ### select only pockets above of minimum percentage ranking



        self.screen_out_subsurf = False
        self.max_desolv_perc = 0.95
        ### screen out pocket if too much of the alpha-cluster ASA is desolvated by the pocket
        ### (i.e. pocket is beneath the surface) (more expensive calculation)



        self.screen_by_res = False
        ### select all pockets containing any atom from any residue within the residue pdb file



        self.pocket_communities = True
        self.tight_option = True
        ### calculate and output the pocket communities
        ### "tight" option uses additional 8.5A core clustering and aux/minor pocket distance cutoff



        self.contact_score = True
        ### calculate contact score features by residue in "lig_resid_file"

        self.beta_cluster = True
        ### recluster the alpha-clusters as a single beta_cluster

        self.beta_cluster_cutoff = 1.5

        self.get_beta_vol = False
        ###calculate and write the beta-cluster volume (total and occupied)

        self.color_table = ['green', 'yellow', 'pink', 'orange', 'blue', 'purple', 'tan', 'olive', 'lime', 'gold', 'aqua', 'rosybrown', 'coral']

    def default(self):
        self.__init__()



    def get_color_RBG(self,color_idx):
        return Color(self.color_table[color_idx%len(self.color_table)]).rgb


if __name__ == '__main__':
    param = AS_Config()
    newcolor = []
    for color in param.color_table:
        try:
            C = Color(color)
            newcolor.append(color)
            print(C)
        except:
            print(color, 'not available')

    print(newcolor)


