from sqlite3 import connect, OperationalError
from progressbar import ProgressBar, Bar, Percentage
from sqlalchemy import create_engine
import os
from util import all_peptides_combinations, list_splitter, get_rdkit_smiles
from FPSim2.io import create_db_file
from dotenv import load_dotenv
from os import getenv
from argparse import ArgumentParser
from util import _create_connection

parser = ArgumentParser()
parser.add_argument('-f', help='sql_file_path for database creation', type=str)
parser.add_argument('-d', help='type of Database SQLite3 supported', type=str)
parser.add_argument('-t', help='sql_table_name | Default-five_pep_aa', type=str)
parser.add_argument('-n', help='length of peptides to create database | Default 5', type=int)
parser.add_argument('-e', help='action to conduct if table exists | Default replace', type=str)
parser.add_argument('-s', help='similarity database name', type=str)
args = parser.parse_args()
load_dotenv(verbose=True)


def create_sqlite_pep_database(sqlite_file_path: str,
                               if_exists: str = 'replace',
                               db: str = 'sqlite',
                               sql_table_name: str = 'five_aa_pep',
                               peptide_len: int = 5
                               ) -> None:
    """

    Args:
        sqlite_file_path:
        if_exists:
        db:
        sql_table_name:
        peptide_len:

    Returns:

    """
    con = _create_connection(file=sqlite_file_path)
    if db == 'psql':
        con = con.connect()
    try:
        if if_exists == 'replace':
            con.execute(f'DROP TABLE IF EXISTS {sql_table_name};')
    except OperationalError:
        os.remove(sqlite_file_path)
        open(sqlite_file_path, 'x')

    schema = f"""
                CREATE TABLE IF NOT EXISTS {sql_table_name}(
                            peptide_id INT,
                            smiles TEXT,
                            peptide CHAR(5),
                            UNIQUE (peptide_id, smiles, peptide),
                            PRIMARY KEY(peptide_id))
            """
    con.execute(schema)
    peptides = all_peptides_combinations(n=peptide_len)
    q = f'INSERT INTO {sql_table_name} (peptide_id, smiles, peptide) VALUES (?, ?, ?)'
    bar = ProgressBar(maxval=100,
                      widgets=[Bar('#', '[', ']'), ' ', Percentage()])

    bar.start()
    count = 0
    for lst in list_splitter(peptides, 1000):
        smiles = map(get_rdkit_smiles, lst)
        for sm, peptide in zip(smiles, lst):
            con.execute(q, [count, sm, peptide])
            count += 1
        percent_completed = (count / len(peptides)) * 100
        bar.update(percent_completed)
    bar.finish()
    con.commit()
    con.close()
    return


def create_similarity_database(sqlite_file_path=None,
                               sql_table_name=None,
                               engine='sqlite',
                               similarity_db='similarity_db.h5'
                               ):
    """

    Args:
        sqlite_file_path:
        sql_table_name:
        engine:
        similarity_db:

    Returns:

    """
    connection = _create_connection(file=sqlite_file_path,
                                    db=engine
                                    )
    if engine == 'postgres':
        c = connection.connect()
    elif engine == 'sqlite':
        c = connection
    else:
        return
    res = c.execute(f'''
                    SELECT smiles, peptide_id FROM {sql_table_name};
                    ''')
    create_db_file(res, similarity_db, 'Morgan', {'radius': 2, 'nBits': 2048})
    connection.close()
    print(f"[Info] --> Created similarity_database for table {sql_table_name} in file {sqlite_file_path}")
    return


create_sqlite_pep_database(args.f,
                           if_exists=args.e,
                           db=args.d,
                           sql_table_name=args.t,
                           peptide_len=args.n
                           )
create_similarity_database(sqlite_file_path=args.f,
                           similarity_db=args.s,
                           engine=args.d,
                           sql_table_name=args.t
                           )
