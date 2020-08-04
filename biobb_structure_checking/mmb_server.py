"""
    Class to support interaction with MMB PDB API to include support to biounits.
"""

import os
import sys

from Bio.PDB.PDBList import PDBList

#from Bio._py3k import urlcleanup as _urlcleanup
#from Bio._py3k import urlretrieve as _urlretrieve
from urllib.request import urlretrieve
from urllib.request import urlcleanup

URL_PREFIX = 'http://mmb.irbbarcelona.org/api/pdb'

class MMBPDBList(PDBList):
    """ Replacement class to support access to biounits at MMB PDB Server.

        Modified from original BioPython code.
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
        if self.pdb_server != 'mmb':
            return super().retrieve_pdb_file(
                pdb_code, obsolete, pdir, file_format, overwrite
            )

        self._verbose = True

        code = pdb_code.lower()

        if file_format in ('pdb', 'mmCif', 'xml'):
            if file_format == 'mmCif':
                file_format = 'cif'
            if not biounit:
                url = (URL_PREFIX + '/%s.%s' % (code, file_format))
            else:
                file_format = 'pdb'
                url = (URL_PREFIX + '/%s_bn%s.pdb' % (code, biounit))
        else:
            print('Error: MMB Server: File format', file_format, 'not supported')
            sys.exit(1)
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
                    print("Structure exists: '%s' " % final_file)
                return final_file

        # Retrieve the file
        if self._verbose:
            if biounit:
                print("Downloading PDB structure '%s.%s'..." % (pdb_code, biounit))
            else:
                print("Downloading PDB structure '%s'..." % pdb_code)
        try:
            urlcleanup()
            urlretrieve(url, final_file)
        except IOError:
            print("Desired structure doesn't exists")
        return final_file
