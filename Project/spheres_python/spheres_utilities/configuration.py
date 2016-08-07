__author__ = 'Idan'
import time
import re

#Default values
N_spheres = 108
n_s = 3
L=1.0
DIAMETER = 0
VOLUME_RATIO = 1.1

debug_mode = False
ZERO_AVG_MOM = True
MASS=1.0
SPEED_SCALE =1.0

speedScale = SPEED_SCALE

INF = 1.7976931348623157*10.**( 308 - 10)# infinity =(1/(10^(-14)))

#paths
RESULTS = "results.txt"
DEBUG = "res/logs/logs_" +re.sub(":","",time.strftime("%j:%H:%M:%S")) + ".txt"

tag=  ""#"_" + re.sub(":","",time.strftime("%H:%M:%S"))
STA_RES = "res/STATUS/status_results"+tag+".json"
COL_RES = "collision_results"+".csv"
POS_RES = "positions_results"+".csv"
VEL_RES = "velocities_results"+".csv"
INIT_RES = "init_results"+".csv"

IMG = "res/img/"

