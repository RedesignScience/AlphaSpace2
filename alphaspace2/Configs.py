"""
This is the default configuration container for AlphaSpace.

The AS_Config is loaded automatically when initializing Trajectory, you can modify the content by
changing the attributes or loading from a new config.ini file.
"""


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
