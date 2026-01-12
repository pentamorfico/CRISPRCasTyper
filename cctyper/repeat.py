import os
import sys
import re
from typing import List

from cctyper.resources import resolve_database_path

class RepeatTyper(object):

    def __init__(self, args) -> None:
       
        self.input: str = args.input
        self.db: str = args.db
        self.threads: int = 1
        self.kmer: int = args.kmer
        self.repeats: List[str] = []

        # Check databases
        self.check_db()
       
        # Read input
        self.read_input()

    def check_db(self) -> None:

        try:
            db_path = resolve_database_path(self.db)
        except RuntimeError as err:
            print(err)
            sys.exit()

        self.db = db_path
        self.xgb = os.path.join(self.db, "xgb_repeats.model")
        self.typedict = os.path.join(self.db, "type_dict.tab")

    def read_input(self) -> None:
        
        # Load input:
        with open(self.input, 'r') as f:
            self.repeats = [ll.rstrip() for ll in f]

        # Check input
        def is_dna(s: str) -> bool:
            match = re.match("^[ACTGactg]*$", s)
            return match is not None

        for rep in self.repeats:
            if not is_dna(rep):
                print('Error - Non-DNA letters found in sequence:')
                print(rep)
                sys.exit()
