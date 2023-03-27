# CRISPR-TE

This is the code base for the CRISPR-TE project

## Content
```
├── PyExtensions
    Contains the Python C Extension code for the modified Trie data structure
├── src
    Contains code for building NGG Trie, CRISPR-TE sgRNA database, and combination calculation
    ├── addtag
    MorenoMateo sgRNA efficiency estimation
    ├── azimuth
    Azimuth sgRNA efficiency estimation
└── website-dev
    Contains codes for the CRISPR-TE website
    ├── api
    Contains the code for the CRISPR-TE API server
```

## PyExtensions
**Installation**
```
cd ./PyExtensions/Trie
python3 setup.py build
```
Then the build library for the Trie data structure is available at `./PyExtensions/Trie/build/lib.linux-x86_64-[PYTHON_VERSION]`

## Building NGG Trie Tree
```python
from src.trie2db import trie2db_1, trie2db_2

PATH_TO_FASTA = "Homo_sapiens.GRCh38.97.dna.primary_assembly.fa"
PATH_TO_ANNOTATION = "hg38_fullAnnotation.bed"
PATH_TO_DATABASE_FILE = "Homo_sapiens.GRCh38.97.dna.primary_assembly.crisprte.db"
PATH_TO_TRIE_FILE = "Homo_sapiens.GRCh38.97.dna.primary_assembly.crisprte.trie.data"

t = trie2db_1(
    PATH_TO_FASTA,
    PATH_TO_DATABASE_FILE,
    PATH_TO_TRIE_FILE
)
trie2db_2(
    t,
    PATH_TO_ANNOTATION,
    PATH_TO_DATABASE_FILE
)
```
Then the `PATH_TO_DATABASE_FILE` contains all information about the NGG sgRNA information and annotations, which can be further accessed using the code from `./src/dbutils.py`.
To build the website serving the database, see `sqlite2postgres.py` to convert the sqlite3 database to PostgresSQL and see `api` to run an API server.

