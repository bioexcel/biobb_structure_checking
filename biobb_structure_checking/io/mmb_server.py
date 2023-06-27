"""
    Class to support interaction with  biounits/assemblies.
"""

import os
import sys
import gzip

from urllib.request import urlretrieve
from urllib.request import urlcleanup

from Bio.PDB.PDBList import PDBList
from biobb_structure_checking.constants import PDB_DOWNLOAD_PREFIX, ALT_SERVERS

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
            biounit=False,  # Kept for compatibility
            nocache=False
    ):
        """
            Replacement for Bio.PDB.PDBList.retrieve_pdb_file to support
            MMB PDB API. Defaults to Biopython super() if standard server
            is used.
        """
        if nocache:
            pdir = '/tmp'
            self.flat_tree = False

        if biounit:
            return self.retrieve_assembly_file(
                pdb_code,
                biounit,
                pdir,
                file_format,
                overwrite
            )

        if self.pdb_server.lower() not in ALT_SERVERS:
            return super().retrieve_pdb_file(
                pdb_code,
                obsolete,
                pdir,
                file_format,
                overwrite
            )
        self._verbose = True

        code = pdb_code.lower()

        if file_format not in ('pdb', 'cif', 'mmCif', 'xml'):
            print(f'Error: MMB/BSC Server: File format {file_format} not supported')
            sys.exit(1)

        if file_format == 'mmCif':
            file_format = 'cif'

        # if not biounit:
        #     url = f'{ALT_SERVERS[self.pdb_server]}/{code}.{file_format}'
        # else:
        #     file_format = 'pdb'
        #     url = f'{ALT_SERVERS[self.pdb_server]}/{code}_bn{biounit}.pdb'
        # Where does the final PDB file get saved?

        url = f'{ALT_SERVERS[self.pdb_server]}/{code}.{file_format}'

        if pdir is None:
            path = self.local_pdb if not obsolete else self.obsolete_pdb
            if not self.flat_tree:  # Put in PDB-style directory tree
                path = os.path.join(path, code[1:3])
        else:  # Put in specified directory
            path = pdir
        if not os.access(path, os.F_OK):
            os.makedirs(path)
        # if biounit:
        #     final = {
        #         'pdb': '%s_%s.pdb',
        #         'mmCif': '%s_%s.cif',
        #         'cif': '%s_%s.cif',
        #         'xml': '%s_%s.xml'
        #     }
        #     final_file = os.path.join(path, final[file_format] % (code, biounit))
        # else:
        #     final = {
        #         'pdb': '%s.pdb',
        #         'mmCif': '%s.cif',
        #         'cif': '%s.cif',
        #         'xml': '%s.xml'
        #     }
        #     final_file = os.path.join(path, final[file_format] % code)
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
            # if biounit:
            #     print(
            #         f"Downloading PDB structure '{pdb_code}.{biounit}' "
            #         f"from {self.pdb_server} ..."
            #     )
            # else:
            #     print(f"Downloading PDB structure '{pdb_code}' from {self.pdb_server} ...")
            print(f"Downloading structure '{pdb_code}' from {self.pdb_server} as {file_format} ...")
        try:
            urlcleanup()
            urlretrieve(url, final_file)
        except IOError:
            print(f"Desired structure doesn't exist at {self.pdb_server} server")
        return final_file

    def retrieve_assembly_file(
            self,
            pdb_code,
            assembly_num,
            pdir=None,
            file_format=None,
            overwrite=False,
            nocache=False
    ):
        """Fetch one or more assembly structures associated with a PDB entry.
        Unless noted below, parameters are described in ``retrieve_pdb_file``.
        :type  assembly_num: int
        :param assembly_num: assembly number to download.
        :rtype : str
        :return: file name of the downloaded assembly file.
        """
        if nocache:
            pdir = '/tmp'
            self.flat_tree = False

#       retrieve_assembly_file only available on biopython >= 1.80, added here to ensure
#           if self.pdb_server.lower() not in ALT_SERVERS:
#                return super().retrieve_assembly_file(
#                    pdb_code, assembly_num, pdir, file_format, overwrite
#                )

        pdb_code = pdb_code.lower()

        assembly_num = int(assembly_num)
        archive = {
            "pdb": f"{pdb_code}.pdb{assembly_num}.gz",
            "mmcif": f"{pdb_code}-assembly{assembly_num}.cif.gz",
        }
        if file_format == 'cif':
            file_format = 'mmcif'

        file_format = self._print_default_format_warning(file_format)
        file_format = file_format.lower()  # we should standardize this.
        if file_format not in archive:
            raise (
                f"Specified file_format '{file_format}' is not supported. "
                " Use one of ('mmcif', 'pdb')."
            )

        # Get the compressed assembly structure name
        archive_fn = archive[file_format]

        # Provisional before fixing local servers
        self.pdb_server = PDB_DOWNLOAD_PREFIX

        if file_format == "mmcif":
            url = self.pdb_server + f"/pub/pdb/data/assemblies/mmCIF/all/{archive_fn}"
        elif file_format == "pdb":
            url = self.pdb_server + f"/pub/pdb/data/biounit/PDB/all/{archive_fn}"
        else:  # better safe than sorry
            raise ValueError(f"file_format '{file_format}' not supported")

        # Where will the file be saved?
        if pdir is None:
            path = self.local_pdb
            if not self.flat_tree:  # Put in PDB-style directory tree
                path = os.path.join(path, pdb_code[1:3])
        else:  # Put in specified directory
            path = pdir
        if not os.access(path, os.F_OK):
            os.makedirs(path)

        assembly_gz_file = os.path.join(path, archive_fn)
        assembly_final_file = os.path.join(path, archive_fn[:-3])  # no .gz

        # Skip download if the file already exists
        if not overwrite:
            if os.path.exists(assembly_final_file):
                if self._verbose:
                    print(f"Structure exists: '{assembly_final_file}' ")
                return assembly_final_file

        # Otherwise,retrieve the file(s)
        if self._verbose:
            print(
                f"Downloading assembly ({assembly_num}) for entry "
                f"'{pdb_code}'..."
            )
        try:
            urlcleanup()
            urlretrieve(url, assembly_gz_file)
        except OSError as err:
            print(f"Download failed! Maybe the desired assembly does not exist: {err}")
        else:
            with gzip.open(assembly_gz_file, "rb") as gz:
                with open(assembly_final_file, "wb") as out:
                    out.writelines(gz)
            os.remove(assembly_gz_file)
        return assembly_final_file