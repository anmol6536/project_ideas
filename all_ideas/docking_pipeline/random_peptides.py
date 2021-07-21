from FPSim2 import FPSim2Engine
from argparse import ArgumentParser
from util import _create_connection
from util import psql_search, random_peptide, fasta_to_smiles


# todo: add parameter to check for the correct table name not just five_pep_aa
def create_random_peptides(fpe, sql_file, n=10, peptides_tested=None, engine='sqlite'):
    """
    Args:
        fpe: Fingerprint similarity engine
        sql_file: SQLite database created in the same format as create database.py recommends
        n: Number of peptides to search
        peptides_tested: list of peptides not consider in search
        engine: Database type, currently only SQLite supported

    Returns:
        list of random dissimilar peptides
    """
    new_peptides = []
    connection = _create_connection(sql_file, db=engine)
    while len(new_peptides) < n:
        peptide = random_peptide(5)
        if peptides_tested:
            while peptide in peptides_tested:
                peptide = random_peptide(5)

        new_peptides.append(peptide)
        query = fasta_to_smiles(new_peptides[-1])
        results = fpe.on_disk_similarity(query, 0)[-2:-1]
        if engine == 'postgres':
            wildcard = '%s'
        elif engine == 'sqlite':
            wildcard = '?'
        else:
            return
        param_results = ",".join(wildcard for _ in range(0, len(results)))
        params = [int(i[0]) for i in results]
        sql_query = f'''SELECT peptide
                        FROM five_aa_pep
                        WHERE peptide_id IN ({param_results})'''
        if peptides_tested:
            param_peptides_tested = ",".join(wildcard for _ in range(0, len(peptides_tested)))
            params.extend(peptides_tested)
            sql_query = f'{sql_query} AND peptide NOT IN ({param_peptides_tested})'

        if engine == 'postgres':
            dissimilar_peptide = psql_search(query=sql_query,
                                             params=params,
                                             con=connection)
        elif engine == 'sqlite':
            dissimilar_peptide = connection.execute(sql_query, params)
        else:
            dissimilar_peptide = []

        new_peptides.extend([i[0] for i in dissimilar_peptide])
    connection.close()
    return new_peptides


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-s', help='FPSim2 similarity File', type=str)
    parser.add_argument('-d', help='SQLite Database with correct format', type=str)
    args = parser.parse_args()

    fingerprint_engine = FPSim2Engine(args.s, in_memory_fps=False)
    res = create_random_peptides(fingerprint_engine, sql_file=args.d)
    print(res)
