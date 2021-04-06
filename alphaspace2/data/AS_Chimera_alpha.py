# AlphaSpace v1.0 - visualization script

import chimera
import sys
import os
from chimera import runCommand as rc  # use 'rc' as shorthand for runCommand
from chimera import replyobj  # for emitting status messages
import glob

pock_dir = 'pockets'

rc("open pdb_out/prot.pdb")
rc("open pdb_out/lig.pdb")
pdb_cnt = 2

for file in glob.glob(pock_dir + "/*_alpha.pdb"):
    rc("open " + file)
pock_cnt = len(glob.glob(pock_dir + "/*_alpha.pdb"))

models_cnt = pock_cnt + pdb_cnt

rc("background solid white")
rc("windowsize 1200 900")
try:
    rc("~longbond")
except:
    pass

rc("colordef my_teal #00a5a5")
rc("colordef my_gold #e6c200")
rc("colordef my_coregreen #4ea24e")
rc("colordef my_lime #bfff00")
rc("colordef my_peri #7d7fc2")
# rc("colordef my_orange #ff8c1a")
rc("colordef my_orange #ff821a")
rc("colordef my_dkpurple #8a65a3")
# rc("colordef my_dkpurple #9748e5")
# rc("colordef my_dkpurple #b148e5")
rc("colordef my_ltgreen #7fc27d")
rc("colordef my_coral #FF4040")

rc("colordef my_pink #de7a8e")
rc("colordef my_blue #3973ac")
# rc("colordef my_blue #336699")
rc("colordef my_ltblue #71acd6")

rc("colordef my_peach #de987a")
rc("colordef my_rasp #8b0045")

rc("colordef my_foam #7dc2c0")

with open('colors_chimera.txt', 'r') as f:
    colors = f.read().splitlines()

colors = colors

rc("~display #0-2")
rc("vdwdefine 1.0 #2")
rc("color white,a #2")

rc("rep stick #0,1")
rc("color white,r,a #1,0")
rc("color byhet,a #1,0")
rc("rib #0,1")

rc("~bond #2-" + str(models_cnt - 1))
rc("rep sphere #2-" + str(models_cnt - 1))

try:
    rc("ribscale 'thin ribbon' #0")
except:
    pass

try:
    rc("ribscale 'all 0.1' #1")
except:
    pass

rc("sel #1 & ~@C,N,O,CA")
rc("sel :AAC za<1.6 & sel")
rc("sel up")
rc("disp sel")
rc("disp :AAC")

rc("rep bs :AAC")
rc("vdwdefine 1.4 :AAC")

rc("display :AAC")
rc("~display :PAC")
rc("color white,a :ACC")
rc("transparency 65,a :ACC")
rc("vdwdefine 1.3 :ACC")
rc("display :ACC")

for i in range(pdb_cnt, models_cnt + 1):
    rc("color " + colors[(i - pdb_cnt) % len(colors)] + ",a #" + str(i) + " & ~:PAC,ACC")
    rc("surfcat surf_p" + str(i) + " #" + str(i) + " & ~:PAC,AAC,ACC")
    rc("surface surf_p" + str(i))
    rc("color " + colors[(i - pdb_cnt) % len(colors)] + ",s #" + str(i))
    rc("~display #" + str(i) + " & ~:PAC,AAC")
    rc("color magenta,a #" + str(i) + ":PAC")
    rc("vdwdefine 1 #" + str(i) + ":PAC")
    rc("transparency 50,a #" + str(i) + " & ~:PAC,AAC")

rc("lighting reflectivity 0")

rc("set shadows")

rc("scene 1 save")

rc("~surf")

rc("~disp :ACC")

for i in range(pdb_cnt, models_cnt + 1):
    rc("display #" + str(i) + " & ~:PAC,AAC,ACC")
    rc("disp #" + str(i) + " zr<0.2 & #0")

rc("scene 2 save")

for i in range(pdb_cnt, models_cnt + 1):
    rc("~display #" + str(i) + " & ~:PAC,AAC,ACC")

rc("color dim grey,a #0")
rc("color light grey,r #0")
rc("color byhet,a #0")

rc("surfcat surf_0 #0")
rc("surface surf_0")
rc("color white,s #0")

rc("lighting reflectivity 100")
rc("lighting reflectivity 0")

rc("color rosy brown,s :.L za<0.2 & #0")
rc("color rosy brown,a :AAC.L")
rc("color my_blue,s :.M za<0.2 & #0")
rc("color my_blue,a :AAC.M")
rc("color my_coregreen,s :.H za<0.2 & #0")
rc("color forest green,a :AAC.H")
rc("disp :AAC,ACC")

rc("scene 3 save")

rc("color light grey,s #0")
rc("~disp :AAC")
rc("color rosy brown,a :ACC.L")
rc("color my_blue,a :ACC.M")
rc("color my_coregreen,a :ACC.H")
rc("vdwdefine 1.1 :ACC")

rc("scene 4 save")

for i in range(pdb_cnt, models_cnt + 1):
    rc("color " + colors[(i - pdb_cnt) % len(colors)] + ",a #" + str(i) + ":ACC")

rc("scene 5 save")

rc("color white,a :ACC")
rc("transparency 65,a :ACC")
rc("vdwdefine 1.3 :ACC")

rc("color white,a @AAO")
rc("color my_coral,a @AAU")
rc("sel @AAO")
rc("sel up")
rc("color my_lime sel & ~@AAO")
rc("disp :AAC")

rc("scene 6 save")

rc("color white,s #0")

for i in range(models_cnt, pdb_cnt - 1, -1):
    rc("sel #" + str(i) + " za<0.2 & #0")
    rc("color " + colors[(i - pdb_cnt) % len(colors)] + ",s sel")
for i in range(pdb_cnt, models_cnt + 1):
    rc("color " + colors[(i - pdb_cnt) % len(colors)] + ",a #" + str(i) + ":AAC")
    rc("disp :AAC")

rc("scene 7 save")

rc("~sel")
