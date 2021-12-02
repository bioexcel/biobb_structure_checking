"""
Module to support json-formatted summary of structure checking, and fixes
"""
import json

class JSONWriter:
    """
    | biobb_structure_checking JSONWriter
    | Classe to manage building and output of JSON summary
    """
    def __init__(self):
        self.data = {}

    def set(self, key, val):
        """ JSONWrite.set
        set single key/value
        Args:
            key (str) : Label
            val (mixed) : Value
        """
        self.data[key] = val

    def append_data(self, array_element, val):
        """ JSONWriter.append_data
        Append item to array
        
        Args:
            array_element(str) : Array to modify
            val (mixed) : Added value
        """
        if array_element not in self.data:
            self.data[array_element] = []
        self.data[array_element].append(val)

    def __str__(self):
        """ JSONWriter.__str__
        String representation
        """
        return json.JSONEncoder(sort_keys=True, indent=1).encode(self.data)

    def save(self, file):
        """ JSONWrite.save
        Save json data
        
        Args:
            file (str) : Path to output file
        """
        with open(file, "w+") as jsout:
            jsout.write(self.__str__())
