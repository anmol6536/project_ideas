from pandas import read_csv, concat
from argparse import ArgumentParser
from os.path import exists


def _create_master(path):
    if not exists(path):
        with open(path, 'w+') as f:
            f.write('IX,sequence,mol,L,NAME,Score,Natom,Nflex,Hbond,Hphob,VwInt,Eintl,Dsolv,SolEl,mfScore,'
                    'RTCNNscore,dTSsc,FILE,POS,CONF,RecConf,O')
            f.close()
    return


def update_refs(master_file, new_results):
    if not exists(master_file):
        _create_master(master_file)
    old_refs = read_csv(master_file)
    new_refs = read_csv(new_results)
    assert all(old_refs.columns == new_refs.columns), "Columns must match for updating reference files"
    old_refs = concat([old_refs, new_refs])
    old_refs.to_csv(master_file, index=False)
    return


def update_good_peptides(results_file, new_results, score=-30):
    if not exists(results_file):
        _create_master(results_file)
    old_refs = read_csv(results_file)
    new_refs = read_csv(new_results)
    assert all(old_refs.columns == new_refs.columns), "Columns must match for updating reference files"
    # filter new results
    new_refs = new_refs[new_refs.Score <= score]
    if not new_refs.empty:
        old_refs = concat([old_refs, new_refs])
        old_refs.to_csv(results_file, index=False)
    return


if __name__ == "__main__":
    desc = """
    Updates previously existing Master files with all tested peptides with new peptides
    """
    epi = ""
    parser = ArgumentParser(description=desc, epilog=epi)
    parser.add_argument('-m', '--master', help='Master file with all references')
    parser.add_argument('-n', '--new', help='New Docking Results')
    args = parser.parse_args()

    update_refs(master_file=args.master, new_results=args.new)
