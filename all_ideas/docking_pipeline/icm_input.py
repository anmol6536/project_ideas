from os.path import exists
from util import convert_one_to_tree_aa_code
from argparse import ArgumentParser


# todo: Append Mode does not work right now
def create_input_icm_docking(list_of_peptides,
                             file_name,
                             if_not_exists='create',
                             if_exists='append',
                             project_folder = '~/'
                             ) -> None:
    """

    Args:
        list_of_peptides:
        file_name:
        if_not_exists:
        if_exists:

    Returns:

    """
    three_letter_list = [*map(convert_one_to_tree_aa_code, list_of_peptides)]
    if not exists(file_name):
        if if_not_exists == 'create':
            open(file_name, 'x')
        else:
            raise FileNotFoundError

    if if_exists == 'append':
        mode = 'a'
    elif if_exists == 'replace':
        mode = 'w'
    else:
        raise ValueError
    with open(f'{project_folder}{file_name}', mode) as f:
        for index, peptide in enumerate(three_letter_list):
            f.write(f'ml a{index + 1}\n')
            f.write(f'se {peptide}\n')

    return


if __name__ == "__main__":
    help_message = """Example Usage: python icm_input.py \
                      -p 'AAAAA,GGGGG' \
                      -o input.se\
                      -c 'create' \
                      -a 'replace' \
                      -f /docking_files"""
    parser = ArgumentParser(description='converts the list of peptides to ICM input file for docking',
                            epilog=help_message
                            )
    parser.add_argument('-p', help='list of peptides to make the input file', type=str)
    parser.add_argument('-o', help='name of the output file', type=str)
    parser.add_argument('-c', help='if does not exist create a new file in working directory', type=str)
    parser.add_argument('-a', help='if exists append to existing file', type=str)
    parser.add_argument('-f', help='project folder', type=str)
    args = parser.parse_args()
    create_input_icm_docking(args.p.split(','),
                             args.o,
                             if_not_exists=args.c,
                             if_exists=args.a,
                             project_folder=args.f
                             )
