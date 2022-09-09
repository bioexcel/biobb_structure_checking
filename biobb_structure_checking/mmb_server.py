"""
    Class to support interaction with MMB PDB API to include support to biounits.
"""

import os
import sys

from urllib.request import urlretrieve
from urllib.request import urlcleanup

from Bio.PDB.PDBList import PDBList

ALT_SERVERS = {
    'mmb': 'http://mmb.irbbarcelona.org/api/pdb',
    'bsc': 'http://mdb-login.bsc.es/api/pdb'
}

class MMBPDBList(PDBList):
    """
    | mmb_server MMBPDBList
    | Replacement class to support access to biounits at MMB PDB Server.
    | Modified from original BioPython code.

    Args:
        pdb_code (str) : PDB code to retrieve
        obsolete (bool) : (False) Retrieve also obsolete entries
        pdir (str) : (None) Path to local cache
        file_format (str) : (None) Format to download (default: mmCIF)
        overwrite (bool): (False) Do not download existing files
        biounit (bool) : (False) Download a Biological assembly
    """

    def retrieve_pdb_file(
            self,
            pdb_code,
            obsolete=False,
            pdir=None,
            file_format=None,
            overwrite=False,
            biounit=False):
        """
            Replacement for Bio.PDB.PDBList.retrieve_pdb_file to support
            MMB PDB API. Defaults to Biopython super() if standard server is used.
        """
        if self.pdb_server.lower() not in ALT_SERVERS:
            return super().retrieve_pdb_file(
                pdb_code, obsolete, pdir, file_format, overwrite
            )

        self._verbose = True

        code = pdb_code.lower()

        if file_format not in ('pdb', 'cif', 'mmCif', 'xml'):
            print(f'Error: MMB/BSC Server: File format {file_format} not supported')
            sys.exit(1)

        if file_format == 'mmCif':
            file_format = 'cif'

        if not biounit:
            url = f'{ALT_SERVERS[self.pdb_server]}/{code}.{file_format}'
        else:
            file_format = 'pdb'
            url = f'{ALT_SERVERS[self.pdb_server]}/{code}_bn{biounit}.pdb'
        #Where does the final PDB file get saved?
        if pdir is None:
            path = self.local_pdb if not obsolete else self.obsolete_pdb
            if not self.flat_tree:  # Put in PDB-style directory tree
                path = os.path.join(path, code[1:3])
        else:  # Put in specified directory
            path = pdir
        if not os.access(path, os.F_OK):
            os.makedirs(path)
        if biounit:
            final = {
                'pdb': '%s_%s.pdb',
                'mmCif': '%s_%s.cif',
                'cif': '%s_%s.cif',
                'xml': '%s_%s.xml'
            }
            final_file = os.path.join(path, final[file_format] % (code, biounit))
        else:
            final = {
                'pdb': '%s.pdb',
                'mmCif': '%s.cif',
                'cif': '%s.cif',
                'xml': '%s.xml'
            }
            final_file = os.path.join(path, final[file_format] % code)

        # Skip download if the file already exists
        if not overwrite:
            if os.path.exists(final_file):
                if self._verbose:
                    print(f"Structure exists: '{final_file}' ")
                return final_file

        # Retrieve the file
        if self._verbose:
            if biounit:
                print(f"Downloading PDB structure '{pdb_code}.{biounit}' from {self.pdb_server} ...")
            else:
                print(f"Downloading PDB structure '{pdb_code}' from {self.pdb_server} ...")
        try:
            urlcleanup()
            urlretrieve(url, final_file)
        except IOError:
            print(f"Desired structure doesn't exist at {self.pdb_server} server")
        return final_file
