from os import popen, getenv, path
from dotenv import load_dotenv
from argparse import ArgumentParser
from time import localtime
from pandas import read_csv, DataFrame
load_dotenv(verbose=True)


def _process_icmdocking_results(results_file, input_file, to_csv=False, csv_name=None):
    results = read_csv(results_file)
    index = []
    with open(input_file, 'r') as f:
        text = f.read()
    compound_lists = [i.split('\n') for i in text.split('ml a') if i]
    for ls in compound_lists:
        sequence = ' '.join(i for i in ls[1:]).replace('se ', '')
        index.append([int(ls[0]), sequence.rstrip().lstrip()])
    index = DataFrame(index, columns=['IX', 'sequence'])
    results = index.merge(results, on='IX')
    if to_csv:
        if csv_name:
            results.to_csv(csv_name, index=False)
            return
        else:
            raise ValueError
    return results


def run_docking(icm_location=None,
                docking_script=None,
                results_processing_script=None,
                project_folder=None,
                input_file='input.se',
                output_file='output',
                project_name='mydock',
                throughness=1.,
                conformations=1
                ) -> None:
    # Default Parameters if not provided
    if not icm_location:
        icm_location = getenv('icm_location')
    if not docking_script:
        docking_script = getenv('docking_script')

    if not path.exists(icm_location):
        raise FileNotFoundError('ICM File does not exist at the provided location')
    if not path.exists(docking_script):
        raise FileNotFoundError('docking_script does not exist at the provided location')

    res = popen(f'''{icm_location} {docking_script} \
                        -s -a thorough={throughness} \
                        input={project_folder}{input_file} \
                        -S \
                        confs={conformations} \
                        {project_folder}{project_name}''')

    # Create output tail
    tail = localtime()
    output_tail = f"""{tail.tm_year}_{str(tail.tm_mon).rjust(2,'0')}_{str(tail.tm_mday).rjust(2,'0')}_{str(tail.tm_hour).rjust(2,'0')}{str(tail.tm_min).rjust(2,'0')}{str(tail.tm_sec).rjust(2,'0')}"""

    with open(f'{project_folder}/log_{output_tail}.log', 'w+') as f:
        f.write(res.read())
    # todo: create a popen and process the docking results
    res = popen(f'''mv {path.abspath('.')}/{project_name}_{input_file.split(".")[0]}*.ob \
             {project_folder}{output_file}.ob
             
             mv {project_folder}{input_file} {project_folder}{input_file.split('.')[0]}_{output_tail}.se
             
             {icm_location} {results_processing_script} \
             {project_folder} {project_name} {output_file} {output_file}{output_tail}
          ''')
    print(res.read())
    _process_icmdocking_results(input_file=f"{project_folder}{input_file.split('.')[0]}_{output_tail}.se",
                                results_file=f'{project_folder}{output_file}{output_tail}.csv',
                                to_csv=True,
                                csv_name=f'{project_folder}results_{output_tail}.csv'
                                )

    return


if __name__ == "__main__":
    help_message = """Example Usage: python run_docking.py -l /Applications/MolsoftICM64.app/Contents/MacOS/icm64 \
                        -d /Applications/MolsoftICM64.app/Contents/Resources/icm/_dockScan \
                        -f /Users/anmol_gorakshakar/python/github/agn_docking/docking_files/ \
                        -i input.se \
                        -o output \
                        -p mydock \
                        -t 0.3 \
                        -c 1
                        """
    parser = ArgumentParser(description='Run ICM Molsoft Docking with limited parameters Optimized for VLS',
                            epilog=help_message
                            )
    parser.add_argument('-i', help='ICM input file', type=str)
    parser.add_argument('-l', help='ICM location', type=str)
    parser.add_argument('-d', help='ICM docking_script location', type=str)
    parser.add_argument('-o', help='name of the output file', type=str)
    parser.add_argument('-p', help='docking project name', type=str)
    parser.add_argument('-f', help='project folder location', type=str)
    parser.add_argument('-a', help='if exists append to existing file', type=str)
    parser.add_argument('-t', help='throughness of the peptide docking', type=float)
    parser.add_argument('-c', help='conformers per peptide to score and save', type=int)
    parser.add_argument('-r', help='ICM results_processing_script location', type=str)
    args = parser.parse_args()
    run_docking(icm_location=args.l,
                docking_script=args.d,
                project_folder=args.f,
                input_file=args.i,
                output_file=args.o,
                project_name=args.p,
                throughness=args.t,
                conformations=args.c,
                results_processing_script=args.r
                )
