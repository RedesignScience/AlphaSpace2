"""
This is the default configuration container for AlphaSpace.

The AS_Config is loaded automatically when initializing AS_Universe, you can modify the content by
changing the attributes or loading from a new config.ini file.
"""




import configparser
import os

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

_DEFAULT_CONFIG_FILE_PATH = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'config.ini')


# noinspection PyAttributeOutsideInit
class AS_Config(object):
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
    # def __new__(cls, config_path=None):
    #     cls.config = configparser.ConfigParser()
    #     cls.config_path = config_path if config_path else _DEFAULT_CONFIG_FILE_PATH
    #     cls.load_config(cls, cls.config_path)
    #     return cls

    def __init__(self,config_path=None):
        self.config = configparser.ConfigParser()
        self.config_path = config_path if config_path else _DEFAULT_CONFIG_FILE_PATH
        self.load_config(self.config_path)


    @staticmethod
    def color(idx):
        """
        Get the RGB of a color corresponding to a given index.
        :param idx: int
        :return: list[r,g,b]
        """
        rounded_idx = int(idx) % len(_COLOR_IDX)
        return _COLOR_DICT[_COLOR_IDX[rounded_idx]]

    @staticmethod
    def color_name(idx):
        """
        Get the name of a color corresponding to a given index.
        :param idx: int
        :return: str
        """
        rounded_idx = int(idx) % len(_COLOR_IDX)
        return _COLOR_IDX[rounded_idx]

    def default(self):
        """
        Reset to default
        :return: AS_Config
        """
        self.__init__()
        return self

    def load_config(self, path):
        """
        Load the configuration from the default config.ini file.
        :param path: str
        :return: bool, if loading was successful
        """
        try:
            self.config_path = path
            self.config.read(self.config_path)

            self.output_dir = self.config.get('options', 'output_dir')
            self.output_to_screen = self.config.getboolean('options', 'output_to_screen')
            self.lig_resid_file = self.config.get('options', 'lig_resid_file')
            self.do_reverse = self.config.getboolean('options', 'do_reverse')
            self.screen_by_face = self.config.getboolean('options', 'screen_by_face')
            self.screen_face_perc = self.config.getfloat('options', 'screen_face_perc')
            self.screen_by_lig_cntct = self.config.getboolean('options', 'screen_by_lig_cntct')
            self.expand_around_cntct = self.config.getboolean('options', 'expand_around_cntct')
            self.screen_by_space = self.config.getboolean('options', 'screen_by_space')
            self.min_space = self.config.getfloat('options', 'min_space')
            self.screen_by_perc_rank = self.config.getboolean('options', 'screen_by_perc_rank')
            self.min_perc_rank = self.config.getfloat('options', 'min_perc_rank')
            self.screen_by_res = self.config.getboolean('options', 'screen_by_res')
            self.res_file = self.config.get('options', 'res_file')
            self.pocket_communities = self.config.getboolean('options', 'pocket_communities')
            self.contact_space = self.config.getboolean('options', 'contact_space')
            self.beta_cluster = self.config.getboolean('options', 'beta_cluster')
            self.get_beta_vol = self.config.getboolean('options', 'get_beta_vol')
            self.min_space = self.config.getfloat('options', 'min_space')
            self.screen_by_perc_rank = self.config.getboolean('options', 'screen_by_perc_rank')
            self.min_perc_rank = self.config.getfloat('options', 'min_perc_rank')
            self.use_asa = self.config.getboolean('options', 'use_asa')
            self.screen_out_subsurf = self.config.getboolean('options', 'screen_out_subsurf')
            self.max_desolv_perc = self.config.getfloat('options', 'max_desolv_perc')

            self.min_r = self.config.getfloat('parameters', 'min_r')
            self.max_r = self.config.getfloat('parameters', 'max_r')
            self.clust_dist = self.config.getfloat('parameters', 'clust_dist')
            self.min_num_alph = self.config.getint('parameters', 'min_num_alph')
            self.hit_dist = self.config.getfloat('parameters', 'hit_dist')
            self.core_cutoff = self.config.getfloat('parameters', 'core_cutoff')
            self.aux_cutoff = self.config.getfloat('parameters', 'aux_cutoff')
            self.tight_communities_cutoff = self.config.getfloat('parameters', 'tight_communities_cutoff')
            self.beta_clust_cutoff = self.config.getfloat('parameters', 'beta_clust_cutoff')
            self.beta_class_cutoff = self.config.getfloat('parameters', 'beta_class_cutoff')
            self.beta_high_cutoff = self.config.getfloat('parameters', 'beta_high_cutoff')
            self.beta_mid_cutoff = self.config.getfloat('parameters', 'beta_mid_cutoff')
            self.lig_clust_dist = self.config.getfloat('parameters', 'lig_clust_dist')
            self.contact_threshold = self.config.getfloat('parameters', 'contact_threshold')
            self.probe_radius = self.config.getfloat('parameters', 'probe_radius')
            self.n_sphere_points = self.config.getfloat('parameters', 'n_sphere_points')
            self.dpocket_cluster_cutoff = self.config.getfloat('parameters', 'dpocket_cluster_cutoff')
            return True
        except:
            print(IOError("Cannot load config {}".format(path)))
            return False

    def write_config(self, path):
        """
        Write the current configuration to a file.
        :param path: str, file location
        :return: bool, if the writing operation is successful
        """
        try:
            self.config.set('options', 'output_dir', self.output_dir)
            self.config.set('options', 'output_to_screen', self.output_to_screen)
            self.config.set('options', 'lig_resid_file', self.lig_resid_file)
            self.config.set('options', 'do_reverse', self.do_reverse)
            self.config.set('options', 'screen_by_face', self.screen_by_face)
            self.config.set('options', 'screen_face_perc', self.screen_face_perc)
            self.config.set('options', 'screen_by_lig_cntct', self.screen_by_lig_cntct)
            self.config.set('options', 'expand_around_cntct', self.expand_around_cntct)
            self.config.set('options', 'screen_by_space', self.screen_by_space)
            self.config.set('options', 'min_space', self.min_space)
            self.config.set('options', 'screen_by_perc_rank', self.screen_by_perc_rank)
            self.config.set('options', 'min_perc_rank', self.min_perc_rank)
            self.config.set('options', 'screen_by_res', self.screen_by_res)
            self.config.set('options', 'res_file', self.res_file)
            self.config.set('options', 'pocket_communities', self.pocket_communities)
            self.config.set('options', 'contact_space', self.contact_space)
            self.config.set('options', 'beta_cluster', self.beta_cluster)
            self.config.set('options', 'get_beta_vol', self.get_beta_vol)
            self.config.set('options', 'min_space', self.min_space)
            self.config.set('options', 'screen_by_perc_rank', self.screen_by_perc_rank)
            self.config.set('options', 'min_perc_rank', self.min_perc_rank)
            self.config.set('options', 'use_asa', self.use_asa)
            self.config.set('options', 'screen_out_subsurf', self.screen_out_subsurf)
            self.config.set('options', 'max_desolv_perc', self.max_desolv_perc)

            self.config.set('parameters', 'min_r', self.min_r)
            self.config.set('parameters', 'max_r', self.max_r)
            self.config.set('parameters', 'clust_dist', self.clust_dist)
            self.config.set('parameters', 'min_num_alph', self.min_num_alph)
            self.config.set('parameters', 'hit_dist', self.hit_dist)
            self.config.set('parameters', 'core_cutoff', self.core_cutoff)
            self.config.set('parameters', 'aux_cutoff', self.aux_cutoff)
            self.config.set('parameters', 'tight_communities_cutoff', self.tight_communities_cutoff)
            self.config.set('parameters', 'beta_clust_cutoff', self.beta_clust_cutoff)
            self.config.set('parameters', 'beta_class_cutoff', self.beta_class_cutoff)
            self.config.set('parameters', 'beta_high_cutoff', self.beta_high_cutoff)
            self.config.set('parameters', 'beta_mid_cutoff', self.beta_mid_cutoff)
            self.config.set('parameters', 'lig_clust_dist', self.lig_clust_dist)
            self.config.set('parameters', 'contact_threshold', self.contact_threshold)
            self.config.set('parameters', 'probe_radius', self.probe_radius)
            self.config.set('parameters', 'n_sphere_points', self.n_sphere_points)
            self.config.set('parameters', 'dpocket_cluster_cutoff', self.dpocket_cluster_cutoff)
            with open(path, 'w') as handle:
                self.config.write(handle)
            return True
        except:
            print("Failed to write configuration to {}".format(path))
