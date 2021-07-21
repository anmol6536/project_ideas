from pandas import read_csv
from sqlite3 import connect
import numpy as np
from mygene import MyGeneInfo as mg
from dotenv import load_dotenv
from os import getenv
from sqlite3 import IntegrityError

load_dotenv()
data_folder = getenv('data_folder')
sql_database = getenv('sqlite_database')
connection = connect(sql_database)
mg = mg()

df = read_csv(f"{data_folder}/protein_expression.csv").dropna()
# filter for true proteins
df = df[[i for i in df if ('Adult' in i) | ('RefSeq' in i)]]
df.columns = [i.lower().replace(' ', '_') for i in df]
df = df[[True if 'NP_' in i else False for i in df.refseq_accession]]
df['refseq_accession'] = df['refseq_accession'].apply(lambda x: x.split('.')[0])

refseq_db = mg.querymany(df.refseq_accession.unique(),
                         scopes='refseq',
                         fields='symbol',
                         species='human',
                         as_dataframe=True)

drop_table = '''drop table if exists refseq_to_symbol'''
connection.execute(drop_table)

refseq_db = refseq_db.reset_index()[['query', 'symbol']]
refseq_table = """create table refseq_to_symbol(
               refseq varchar,
               gene varchar,
               primary key (gene)
               );"""
connection.execute(refseq_table)

for idx, values in refseq_db.iterrows():
    try:
        connection.execute(f"""insert into refseq_to_symbol(refseq, gene) 
                               values ('{values['query']}', '{values['symbol']}')""")
        connection.commit()
    except IntegrityError:
        pass

