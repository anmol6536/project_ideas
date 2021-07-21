from pandas import read_csv
from random_peptides import create_random_peptides as crp
from argparse import ArgumentParser
from icm_input import create_input_icm_docking as ciid
from FPSim2 import FPSim2Engine
from run_docking import run_docking as rd
from update_results import update_refs as ur
from update_results import update_good_peptides as ugp
from util import three_to_one_letter_aa
from os import getenv, listdir, path
from dotenv import load_dotenv
from similiar_peptides import create_similar_peptides as csp
load_dotenv('./.env')
working_dir = path.abspath('.')


def montecarlo_peptides(sim_search=False, search_pep=None):
    icm_location = getenv('icm_location')
    docking_script = getenv('docking_script')
    project_folder = getenv('project_folder')
    results_processing_script = getenv('results_processing_sctipt')
    master_peptides_file = f"{project_folder}/master.csv"
    good_peptides_file = f"{project_folder}/good_pep.csv"
    fp_filename = getenv('fp_filename')
    sql_file = getenv('sql_file')
    input_file = getenv('input_file')
    project_name = getenv('project_name')
    score = int(getenv('score_filter'))
    thoroughness=float(getenv('thoroughness'))
    mnCalls = int(getenv('mnCalls'))
    calls_completed = 0

    fpe = FPSim2Engine(fp_filename, in_memory_fps=False)  # initialize similarity engine

    peptides_tested_df = read_csv(master_peptides_file)  # check for already tested peptides
    peptides_tested = [*map(three_to_one_letter_aa, peptides_tested_df.sequence)]

    if sim_search:
        find_similar_peptides = True
        assert search_pep is not None, "Peptides must be provided with sim_search flasg"
        good_peptides = search_pep
    else:
        find_similar_peptides = False
        good_peptides = []

    while calls_completed <= mnCalls:
        if not find_similar_peptides:
            input_peptides = crp(fpe=fpe,
                                 peptides_tested=peptides_tested,
                                 sql_file=sql_file,
                                 n=10,
                                 engine='sqlite')
            print(input_peptides)
        else:
            if len(good_peptides) > 0:
                input_peptides = csp(fpe=fpe,
                                     similar_peptides_for=good_peptides,
                                     peptides_tested=peptides_tested,
                                     sql_file=sql_file, n=20,
                                     engine='sqlite')

        ciid(list_of_peptides=input_peptides,
             file_name=input_file,
             if_not_exists='create',
             if_exists='replace',
             project_folder=project_folder
             )

        rd(icm_location=icm_location,
            docking_script=docking_script,
            results_processing_script=results_processing_script,
            project_folder=project_folder,
            input_file=input_file,
            output_file='output',
            project_name=project_name,
            throughness=thoroughness,
            conformations=1)

        results_file = sorted([i for i in listdir(project_folder) if 'results' in i])[-1]
        print(results_file)
        ur(master_file=master_peptides_file,
           new_results=f"{project_folder}/{results_file}")  # update results
        ugp(results_file=good_peptides_file,
            new_results=f"{project_folder}/{results_file}",
            score=score)  # update good peptides

        good_peptides_df = read_csv(f"{project_folder}/{results_file}")
        good_peptides = [*good_peptides_df[good_peptides_df.Score <= score].sequence]
        if len(good_peptides) > 0:
            find_similar_peptides = True
            good_peptides = [*map(three_to_one_letter_aa, good_peptides)]
            score += score*0.2
        else:
            find_similar_peptides = False
            good_peptides = []
        peptides_tested_df = read_csv(master_peptides_file)
        peptides_tested = [*map(three_to_one_letter_aa, peptides_tested_df.sequence)]

        calls_completed += 1

    return


if __name__ == "__main__":
    parser = ArgumentParser(description="", epilog="")
    # parser.add_argument()
    montecarlo_peptides()
