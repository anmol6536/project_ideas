from sqlite3 import connect
from sqlalchemy import create_engine
from os import getenv
from random import choice
from rdkit.Chem import MolFromFASTA, MolToSmiles, RDKFingerprint, MolFromSmiles
from rdkit.DataStructs import CosineSimilarity

three_letter_code = {'A': 'ala',
                     'R': 'arg',
                     'N': 'asn',
                     'D': 'asp',
                     'C': 'cys',
                     'E': 'glu',
                     'Q': 'gln',
                     'G': 'gly',
                     'H': 'his',
                     'I': 'ile',
                     'L': 'leu',
                     'K': 'lys',
                     'M': 'met',
                     'F': 'phe',
                     'P': 'pro',
                     'S': 'ser',
                     'T': 'thr',
                     'W': 'trp',
                     'Y': 'tyr',
                     'V': 'val'}
one_letter_code = dict([(value, key) for key, value in three_letter_code.items()])
three_to_one_letter_aa = lambda x: ''.join(one_letter_code[i] for i in x.split(' ') if len(x.replace(' ', '')) > 1)

aa = ['A', 'R', 'N', 'D', 'C', 'E', 'Q',
      'G', 'H', 'I', 'L', 'K', 'M', 'F',
      'P', 'S', 'T', 'W', 'Y', 'V']


def _create_connection(file, db='sqlite', autocommit=False):
    """
    Create SQL Engine to connect to SQLite or PSQL Database
    Args:
        file: (string) abs path to sqlite database
        db: (string) sqlite | psql

    Returns:
        Live SQL Engine for sqlite or psql
    """
    if db == 'sqlite':
        connection = connect(file)
    elif db == 'psql':
        connection = create_engine(getenv('psql_agn_db'))
        if autocommit:
            connection.autocommit = True
    else:
        connection = 404
    return connection


def convert_one_to_tree_aa_code(peptide):
    if type(peptide) == list:
        holder = []
        for pep in peptide:
            try:
                assert type(pep) == str, 'please enter peptides as string'
                pep = ' '.join(three_letter_code[i.upper()] for i in pep)
                holder.append(pep)
            except Exception as e:
                if e.__class__ == AssertionError:
                    print(f'{pep} raised AssertionError')
                elif e.__class__ == KeyError:
                    print('peptide raised KeyError: Check peptide sequence')
                return e

        return holder

    elif type(peptide) == str:
        try:
            pep = ' '.join(three_letter_code[i.upper()] for i in peptide)
        except Exception as e:
            if e.__class__ == KeyError:
                print('peptide raised KeyError: Check peptide sequence')
            else:
                print(e)
            return
        return pep


random_peptide = lambda x: ''.join(choice([*three_letter_code.keys()]) for i in range(x))
fasta_to_smiles = lambda x: MolToSmiles(MolFromFASTA(x))


def process_icmdocking_results(results_file, input_file):
    results = pd.read_csv(results_file)
    index = []
    with open(input_file, 'r') as f:
        text = f.read()
    compound_lists = [i.split('\n') for i in a.split('ml a') if i]
    for ls in compound_lists:
        sequence = ' '.join(i for i in ls[1:]).replace('se ', '')
        index.append([int(ls[0]), sequence.rstrip().lstrip()])
    index = pd.DataFrame(index, columns=['IX', 'sequence'])
    results = index.merge(results, on='IX')
    return results


def all_peptides_combinations(n=1):
    aminoacids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M',
                  'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    current_aa_length = 0
    list_of_peptides = []
    while current_aa_length < n:
        if list_of_peptides:
            new_list_of_peptides = []
            for peptide in list_of_peptides:
                for aa in aminoacids:
                    new_peptide = peptide + aa
                    new_list_of_peptides.append(new_peptide)
            list_of_peptides = new_list_of_peptides
            current_aa_length += 1
        else:
            list_of_peptides = aminoacids
            current_aa_length += 1
    return list_of_peptides


def cosine_similarity(molecule_smiles=None, peptide_length=3, cl=None):
    peptides = all_peptides_combinations(peptide_length)

    rdkit_query_molecules = [*map(MolFromSmiles,
                                  molecule_smiles)]
    rdkit_query_fingerprint = [*map(RDKFingerprint,
                                    rdkit_query_molecules)]

    rdkit_peptide_molecules = [*map(MolFromFASTA,
                                    peptides)]
    rdkit_peptide_smiles = [*map(MolToSmiles,
                                 rdkit_peptide_molecules)]
    rdkit_peptide_fingerprint = [*map(RDKFingerprint,
                                      rdkit_peptide_molecules)]
    permutation = []
    for query_fp, mol_smiles, cluster in zip(rdkit_query_fingerprint,
                                             molecule_smiles,
                                             cl):
        for peptide_fp, peptide_smiles, peptide_symbol in zip(rdkit_peptide_fingerprint,
                                                              rdkit_peptide_smiles,
                                                              peptides):
            info = {
                'query_molecule': mol_smiles,
                'peptide_to_compare': peptide_smiles,
                'peptide_symbol': peptide_symbol,
                'similarity': CosineSimilarity(query_fp, peptide_fp),
                'cluster_no': cluster
            }
            permutation.append(info)
    return permutation


def list_splitter(ls=[], n=100):
    for i in range(0, len(ls), n):
        yield ls[i:i + n]


def get_rdkit_smiles(string=None, kind='peptide'):
    assert type(string) == str
    if string.upper() != string:
        string = string.upper()
    if kind == 'peptide' and string != None:
        try:
            return MolToSmiles(MolFromFASTA(string.upper()))
        except Exception as e:
            non_standard_aa = [i for i in string if i not in aa]
            print(non_standard_aa, 'is / are not a canonical aa')


def mutate_peptide(peptide, pos=None, mutate_with=None):
    aa = ['A', 'R', 'N', 'D', 'C', 'E', 'Q',
          'G', 'H', 'I', 'L', 'K', 'M', 'F',
          'P', 'S', 'T', 'W', 'Y', 'V']
    if not mutate_with:
        mutate_with = aa
    peptide = list(peptide)
    for i in peptide:
        if i.upper() not in aa:
            raise ValueError(f'{i.upper()} is not an amino acid')
    mutated_peptides = []
    for i in mutate_with:
        peptide[pos] = i
        mutated_peptides.append(''.join(i.upper() for i in peptide))
    return mutated_peptides


def histlist_entry(compound='',
                   peptide='',
                   library='',
                   score=0,
                   target='trka-d5',
                   name='NA',
                   sub_region='nterm'):
    holder = {
        'compound': compound,
        'score': score,
        'sub_region': sub_region,
        'target': target,
        'library': library,
        'peptide': peptide,
        'name': name,
    }
    return holder


# PSQL Searches:
def psql_execute(query, con, params=None, autocommit=False):
    if autocommit:
        con.autocommit = autocommit
    with con.connect() as c:
        c.execute(query, params)
        c.close()


def psql_search(query, con, params=None):
    with con.connect() as c:
        res = c.execute(query, params)
        c.close()
    return [*res]


def random_peptide(n):
    return ''.join(choice([*three_letter_code.keys()]) for _ in range(n))


def fasta_to_smiles(seq):
    return MolToSmiles(MolFromFASTA(seq))
