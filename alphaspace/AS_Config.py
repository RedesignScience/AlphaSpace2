from colour import Color


class AS_Config(object):
    def __new__(cls, *args, **kwargs):
        # Alpha Sphere screening upper and lower limits
        cls.min_r = 3.2
        cls.max_r = 5.4

        # clustering distance cutoff atom alpha atom to pocket
        cls.pocket_cluster_distance = 4.7

        # Contact cutoff threshold
        cls.contact_threshold = 3.0


        # Screening options
        ### select only pockets in "contact" with the ligand
        cls.screen_by_lig_cntct = True

        ### expand selection to include all pockets overlapping the "contact" pockets
        cls.expand_around_cntct = False

        ### select pockets that include atoms from the interface atom list
        cls.screen_by_face = True

        ### minimum percentage of atoms that must be in the interface atom list
        cls.screen_face_perc = 0.5



        cls.screen_by_seam = False
        ### select pocket that contain atoms from >1 protein chain
        cls.expand_around_seam = False
        ### exapnd selection to include all pockets overlapping the "seam" pockets


        cls.screen_by_score = False
        cls.min_score = 100
        ### select only pockets above a minimum score



        cls.screen_by_perc_rank = False
        cls.min_perc_rank = 0.9
        ### select only pockets above of minimum percentage ranking


        cls.screen_out_subsurf = False
        cls.max_desolv_perc = 0.95
        ### screen out pocket if too much of the alpha-cluster ASA is desolvated by the pocket
        ### (i.e. pocket is beneath the surface) (more expensive calculation)

        cls.screen_by_res = False
        ### select all pockets containing any atom from any residue within the residue pdb file


        cls.pocket_communities = True
        cls.tight_option = True
        ### calculate and output the pocket communities
        ### "tight" option uses additional 8.5A core clustering and aux/minor pocket distance cutoff

        cls.contact_score = True
        ### calculate contact score features by residue

        cls.beta_cluster = True
        ### recluster the alpha-clusters as a single beta_cluster

        cls.beta_cluster_cutoff = 1.5

        cls.get_beta_vol = False
        ###calculate and write the beta-cluster volume (total and occupied)

        cls.color_table = ['green', 'yellow', 'pink', 'orange', 'blue', 'purple', 'tan', 'olive', 'lime', 'gold',
                            'aqua', 'rosybrown', 'coral']

        cls.dpocket_cluster_cutoff = 0.75

    def default(self):
        self.__init__()




if __name__ == '__main__':
    config = AS_Config()
    newcolor = []
    for color in config.color_table:
        C = Color(color)
        newcolor.append(color)
        print(C.rgb)
