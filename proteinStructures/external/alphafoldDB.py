from Bio import SeqIO
import requests
from pathlib import Path
from tqdm import tqdm
from typing import Any
import logging
import json
import asyncio
import aiohttp
import time

# Typing hints
Result = dict[str, Any]
Results = list[Result]

def read_dataset(filename : str) -> list[SeqIO.SeqRecord]:
    """Read fasta file and return list of SeqRecords"""
    filename = Path(filename)
    if not filename.exists():
        raise ValueError(f"File {filename} does not exist")
    records = list(SeqIO.parse(filename, "fasta"))
    return records

async def get_uniprot_acc(uniprot_acc : str, session,url_base=f"https://alphafold.ebi.ac.uk/api/prediction/") -> Tuple[str,Results]:
    """Query AlphaFoldDB for information on uniprot_acc
    
    Args:
    - uniprot_acc: str, uniprot accession number

    Returns:
    - uniprot_acc: str, uniprot accession number
    - Results : List[dict], list of results from AlphaFoldDB
    """
    url = url_base + uniprot_acc
    try:
        async with session.get(url=url,timeout=5000) as response:
            resp = await response.read()
            print("Successfully got url {} with resp of length {}.".format(url, len(resp)))
            return (uniprot_acc,json.loads(resp.decode('utf-8')))
    except Exception as e:
        print("Unable to get url {} due to {}: {}.".format(url, e.__class__, str(e)))
        return (uniprot_acc,[])

async def main_async_queries(uniprot_accs)-> List[Tuple[str,Results]]:
    """Run async queries for all uniprot_accs
    Args:
    - uniprot_accs: List[str], list of uniprot accession numbers
    """
    async with aiohttp.ClientSession(headers={"Accept": "application/json"}) as session:
        ret = await asyncio.gather(*[get_uniprot_acc(uniprot_acc, session) for uniprot_acc in uniprot_accs])
    return ret

def run_alphafold_queries(filename : str,limit:int|None=None) -> List[Tuple[str,Results]]:
    """Query AlphaFoldDB for all uniprot accessions in given fasta file
    
    Args:
    - filename: str, path to fasta file

    Optional Args:
    - limit : int | None : limit the number of queries to run, for testing purposes
    Returns:
    - 
    """
    records = read_dataset(filename)
    uniprot_accs = [record.id for record in records]
    if limit is not None:
        uniprot_accs = uniprot_accs[:limit]
    t0 = time.time()
    results = asyncio.run(main_async_queries(uniprot_accs))
    t1 = time.time()
    print(f"Time elapsed: {t1-t0}")
    return results

if __name__ == "__main__":
    ret = run_alphafold_queries("/home/dzeiberg/proteinStructures/data/swissprot_dataset_fullheader.fasta",limit=None)
    rdict = dict(ret)
    with open("/home/dzeiberg/proteinStructures/data/alphafoldDB_results.json", "w") as f:
        json.dump(rdict, f, indent=4)