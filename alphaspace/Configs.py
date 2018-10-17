"""
This is the default configuration container for AlphaSpace.

The AS_Config is loaded automatically when initializing Trajectory, you can modify the content by
changing the attributes or loading from a new config.ini file.
"""

import configparser
import os

_COLOR_DICT = {'blue': [0.22265625, 0.44921875, 0.671875],
               'coral': [0.99609375, 0.25, 0.25],
               'coregreen': [0.3046875, 0.6328125, 0.3046875],
               'dkpurple': [0.5390625, 0.39453125, 0.63671875],
               'foam': [0.48828125, 0.7578125, 0.75],
               'gold': [0.8984375, 0.7578125, 0.0],
               'lime': [0.74609375, 0.99609375, 0.0],
               'ltblue': [0.44140625, 0.671875, 0.8359375],
               'ltgreen': [0.49609375, 0.7578125, 0.48828125],
               'orange': [0.99609375, 0.5078125, 0.1015625],
               'peach': [0.8671875, 0.59375, 0.4765625],
               'peri': [0.48828125, 0.49609375, 0.7578125],
               'pink': [0.8671875, 0.4765625, 0.5546875],
               'rasp': [0.54296875, 0.0, 0.26953125],
               'teal': [0.0, 0.64453125, 0.64453125]}
_COLOR_IDX = {0: "teal",
              1: "gold",
              2: "coregreen",
              3: "lime",
              4: "peri",
              5: "orange",
              6: "dkpurple",
              7: "ltgreen",
              8: "coral",
              9: "pink",
              10: "blue",
              11: "ltblue",
              12: "peach", }

_DEFAULT_CONFIG_FILE_PATH = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'config.ini')


# noinspection PyAttributeOutsideInit
class AS_Config(object):
    _COLOR_DICT = {'blue': [0.22265625, 0.44921875, 0.671875],
                   'coral': [0.99609375, 0.25, 0.25],
                   'coregreen': [0.3046875, 0.6328125, 0.3046875],
                   'dkpurple': [0.5390625, 0.39453125, 0.63671875],
                   'foam': [0.48828125, 0.7578125, 0.75],
                   'gold': [0.8984375, 0.7578125, 0.0],
                   'lime': [0.74609375, 0.99609375, 0.0],
                   'ltblue': [0.44140625, 0.671875, 0.8359375],
                   'ltgreen': [0.49609375, 0.7578125, 0.48828125],
                   'orange': [0.99609375, 0.5078125, 0.1015625],
                   'peach': [0.8671875, 0.59375, 0.4765625],
                   'peri': [0.48828125, 0.49609375, 0.7578125],
                   'pink': [0.8671875, 0.4765625, 0.5546875],
                   'rasp': [0.54296875, 0.0, 0.26953125],
                   'teal': [0.0, 0.64453125, 0.64453125]}

    _COLOR_IDX = {0: "teal",
                  1: "gold",
                  2: "coregreen",
                  3: "lime",
                  4: "peri",
                  5: "orange",
                  6: "dkpurpl",
                  7: "ltgreen",
                  8: "coral",
                  9: "pink",
                  10: "blue",
                  11: "ltblue",
                  12: "peach", }

    # def __new__(cls, config_path=None):
    #     cls.config = configparser.ConfigParser()
    #     cls.config_path = config_path if config_path else _DEFAULT_CONFIG_FILE_PATH
    #     cls.load_config(cls, cls.config_path)
    #     return cls

    def __init__(self, config_path=None):
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
            self.hdbscan_min_samples = self.config.getfloat('parameters', 'hdbscan_min_samples')

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
