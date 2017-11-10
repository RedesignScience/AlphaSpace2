import configparser,os,sys

_COLOR_DICT = dict(green=[0.0, 0.5019607843137255, 1.6718652606120004e-16], yellow=[1.0, 0.9999999999999998, 0.0],
                   pink=[1.0, 0.7529411764705882, 0.7960784313725489], orange=[1.0, 0.6470588235294115, 0.0],
                   blue=[0.0, 0.0, 1.0], purple=[0.5019607843137255, 0.0, 0.5019607843137254],
                   tan=[0.8235294117647057, 0.7058823529411764, 0.5490196078431374],
                   olive=[0.5019607843137255, 0.5019607843137254, 0.0], lime=[0.0, 1.0, 3.3306690738754696e-16],
                   gold=[1.0, 0.8431372549019605, 0.0], aqua=[0.0, 0.9999999999999998, 1.0],
                   rosybrown=[0.7372549019607844, 0.5607843137254902, 0.5607843137254902],
                   coral=[0.9999999999999999, 0.49803921568627463, 0.31372549019607854])

_COLOR_IDX = {0: "green",
             1: "yellow",
             2: "pink",
             3: "orange",
             4: "blue",
             5: "purple",
             6: "tan",
             7: "olive",
             8: "lime",
             9: "gold",
             10: "aqua",
             11: "rosybrown",
             12: "coral"}

_DEFAULT_CONFIG_FILE_PATH = "./config.ini"


class AS_Config(object):
    def __new__(cls,config_path = None):
        cls.config = configparser.ConfigParser()
        cls.config_path = config_path if config_path else _DEFAULT_CONFIG_FILE_PATH
        cls.load_config(cls,cls.config_path)
        return cls

    def color(self,idx):
        assert type(idx) == int
        rounded_idx = idx%len(_COLOR_IDX)
        return _COLOR_DICT[_COLOR_IDX[rounded_idx]]

    def default(self):
        self.__init__()

    def load_config(self,path):

        self.config_path = path
        self.config.read(self.config_path)

        self.output_dir =           self.config.get('options', 'output_dir')
        self.output_to_screen =     self.config.getboolean('options', 'output_to_screen')
        self.lig_resid_file =       self.config.get('options', 'lig_resid_file')
        self.do_reverse =           self.config.getboolean('options', 'do_reverse')
        self.screen_by_face =       self.config.getboolean('options', 'screen_by_face')
        self.screen_face_perc =     self.config.getfloat('options', 'screen_face_perc')
        self.screen_by_lig_cntct =  self.config.getboolean('options', 'screen_by_lig_cntct')
        self.expand_around_cntct =  self.config.getboolean('options', 'expand_around_cntct')
        self.screen_by_score =      self.config.getboolean('options', 'screen_by_score')
        self.min_score =            self.config.getfloat('options', 'min_score')
        self.screen_by_perc_rank =  self.config.getboolean('options', 'screen_by_perc_rank')
        self.min_perc_rank =        self.config.getfloat('options', 'min_perc_rank')
        self.screen_by_res =        self.config.getboolean('options', 'screen_by_res')
        self.res_file =             self.config.get('options', 'res_file')
        self.pocket_communities =   self.config.getboolean('options', 'pocket_communities')
        self.contact_score =        self.config.getboolean('options', 'contact_score')
        self.beta_cluster =         self.config.getboolean('options', 'beta_cluster')
        self.get_beta_vol =         self.config.getboolean('options', 'get_beta_vol')
        self.min_score =            self.config.getfloat('options', 'min_score')
        self.screen_by_perc_rank =  self.config.getboolean('options', 'screen_by_perc_rank')
        self.min_perc_rank =        self.config.getfloat('options', 'min_perc_rank')

        self.min_r =                    self.config.getfloat('parameters','min_r')
        self.max_r =                    self.config.getfloat('parameters','max_r')
        self.clust_dist =               self.config.getfloat('parameters','clust_dist')
        self.min_num_alph =             self.config.getint('parameters','min_num_alph')
        self.hit_dist =                 self.config.get('parameters','hit_dist')
        self.core_cutoff =              self.config.getfloat('parameters','core_cutoff')
        self.aux_cutoff =               self.config.getfloat('parameters','aux_cutoff')
        self.tight_communities_cutoff=  self.config.getfloat('parameters','tight_communities_cutoff')
        self.beta_clust_cutoff =        self.config.getfloat('parameters','beta_clust_cutoff')
        self.beta_class_cutoff =        self.config.getfloat('parameters','beta_class_cutoff')
        self.beta_high_cutoff =         self.config.getfloat('parameters','beta_high_cutoff')
        self.beta_mid_cutoff =          self.config.getfloat('parameters','beta_mid_cutoff')
        self.lig_clust_dist =           self.config.getfloat('parameters','lig_clust_dist')
        self.contact_threshold =        self.config.getfloat('parameters','contact_threshold')

    def write_config(self,path):
        self.config.set('options', 'output_dir', self.output_dir)
        self.config.set('options', 'output_to_screen', self.output_to_screen)
        self.config.set('options', 'lig_resid_file', self.lig_resid_file)
        self.config.set('options', 'do_reverse', self.do_reverse)
        self.config.set('options', 'screen_by_face', self.screen_by_face)
        self.config.set('options', 'screen_face_perc', self.screen_face_perc)
        self.config.set('options', 'screen_by_lig_cntct', self.screen_by_lig_cntct)
        self.config.set('options', 'expand_around_cntct', self.expand_around_cntct)
        self.config.set('options', 'screen_by_score', self.screen_by_score)
        self.config.set('options', 'min_score', self.min_score)
        self.config.set('options', 'screen_by_perc_rank', self.screen_by_perc_rank)
        self.config.set('options', 'min_perc_rank', self.min_perc_rank)
        self.config.set('options', 'screen_by_res', self.screen_by_res)
        self.config.set('options', 'res_file', self.res_file)
        self.config.set('options', 'pocket_communities', self.pocket_communities)
        self.config.set('options', 'contact_score', self.contact_score)
        self.config.set('options', 'beta_cluster', self.beta_cluster)
        self.config.set('options', 'get_beta_vol', self.get_beta_vol)
        self.config.set('options', 'min_score', self.min_score)
        self.config.set('options', 'screen_by_perc_rank', self.screen_by_perc_rank)
        self.config.set('options', 'min_perc_rank', self.min_perc_rank)


        with open(path,'w') as handle:
            self.config.write(handle)

