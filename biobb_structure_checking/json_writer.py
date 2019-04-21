"""
Module to support json-formatted summary of checking, and fixes
"""
import json

class JSONWriter:
    """
    Classe to manage building and output of JSON summary
    """
    def __init__(self):
        self.data = {}

    def set(self, key, val):
        """set single key/value"""
        self.data[key] = val

    def append_data(self, array_element, val):
        """append item to array"""
        if array_element not in self.data:
            self.data[array_element] = []
        self.data[array_element].append(val)

    def __str__(self):
        return json.JSONEncoder(sort_keys=True, indent=1).encode(self.data)

    def save(self, file):
        """save json data"""
        with open(file, "w+") as jsout:
            jsout.write(self.__str__())
