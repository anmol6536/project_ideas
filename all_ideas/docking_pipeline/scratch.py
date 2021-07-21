import os
from pandas import read_csv
from update_results import update_refs as ur
from update_results import update_good_peptides as ugp
from similiar_peptides import create_similar_peptides as csp
from FPSim2 import FPSim2Engine

fp_filename='similaritydb.h5'
sql_file='similarity.db'
fpe = FPSim2Engine(fp_filename, in_memory_fps=False)
res = csp(fpe=fpe,
     similar_peptides_for = ['AAAAA', 'KKKKK'],
     sql_file=sql_file, n=5,
     engine='sqlite')
print(res)
# res = [i for i in os.listdir('docking_files') if 'results' in i]
# project_folder = os.path.abspath('./docking_files')
# print(project_folder)
# for res_f in res:
#     ur('docking_files/master.csv', f"{project_folder}/{res_f}")
#     ugp('docking_files/good_pep.csv', f"{project_folder}/{res_f}", -10)

