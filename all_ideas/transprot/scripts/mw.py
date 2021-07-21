import numpy as np
from requests import get
from json import loads
from scipy.stats import gmean
import bs4
from sqlite3 import connect
from os import getenv, path
from dotenv import load_dotenv
from sqlite3 import OperationalError
from Bio.SeqUtils import molecular_weight
from mygene import MyGeneInfo
import matplotlib.pyplot as plt

mg = MyGeneInfo()

def create_connection(sql_file):
    return connect(sql_file)


def _validate_mw_table(sql_file, create_file=False, create_table=False):
    if create_file:
        if not path.exists(sql_file):
            with open(sql_file, 'x') as f:
                print('sql_file_created')
                f.close()
    if path.exists(sql_file):
        con = create_connection(sql_file)
    else:
        raise FileNotFoundError('try the validation with create_file flag')

    if create_table:
        if_exists = 'replace'
        if if_exists == 'replace':
            drop_original = 'drop table if exists mw;'
            con.execute(drop_original)

        table_schema = """
                       create table mw(
                       gene varchar(20),
                       mw_da float,
                       primary key (gene),
                       unique (gene, mw_da)
                       );
                       """
        con.execute(table_schema)
        print('table created for molecular mw')

    query = "select * from sqlite_master where tbl_name='mw';"
    res = con.execute(query)
    res = res.fetchone()
    if not res:
        raise OperationalError('table mw not found')
    return


def get_uniprot(gene):
    uniprot_id = mg.query(gene.upper(), scopes='symbol', fields='uniprot', species='human')
    uniprot_id = [(i['_score'], i['uniprot']['Swiss-Prot']) for i in uniprot_id['hits'] if 'uniprot' in i.keys()]
    if uniprot_id:
        best_score = [i[0] for i in uniprot_id]
        best_match = [i[1] for i in uniprot_id if i[0] == max(best_score)]
        return best_match
    else:
        raise ValueError('No Uniprot Found')


def mw_finder(gene):
    sql_file = getenv('sqlite_database')
    connection = create_connection(sql_file)
    result = connection.execute(f"select * from mw where gene in ('{gene}');")
    result = result.fetchall()

    if result:
        return result

    else:
        uniprot_ids = get_uniprot(gene)
        sequences = []
        for uniprot_id in uniprot_ids:
            uniprot_info = get(f'https://www.uniprot.org/uniprot/{uniprot_id}.fasta')
            if uniprot_info.ok:
                sequence = ''.join(uniprot_info.text.split('\n')[1:])
                sequences.append(sequence)
        mass = [*map(molecular_weight, sequences, ['protein' for _ in range(len(sequences))])]
        # return masses
        if len(mass) == 1:
            connection.execute(f"insert into mw('gene', 'mw_da') VALUES ('{gene}', '{mass[0]}');")
            connection.commit()
            connection.close()

    return [(gene, mass[0])]


if __name__ == "__main__":
    load_dotenv()
    sql_file = getenv('sqlite_database')
    _validate_mw_table(sql_file,
                       create_file=False,
                       create_table=False)
    res = mw_finder('NTRK1')
    print(res)