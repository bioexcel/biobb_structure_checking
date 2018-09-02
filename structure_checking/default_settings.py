"""
Class to initialize global settings like data file path.
"""
class DefaultSettings():
    def __init__(self, base_dir):
        self.base_dir_path = base_dir
#Default locations
        self.help_dir_path = self.base_dir_path + "/helpdocs"
        self.data_dir_path = self.base_dir_path + "/dat"
        self.res_library_path = self.data_dir_path + "/all_amino03.in"
        self.data_library_path = self.data_dir_path + "/data_lib.json"

