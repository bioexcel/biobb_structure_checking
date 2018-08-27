
import json

class JSONWriter:
    def __init__(self):
        self.data = {}

    def set(self,k,v):
        self.data[k]=v

    def appendData(self, ar, v):
        if ar not in self.data:
            self.data[ar]=[]
        self.data[ar].append(v)

    def __str__(self):
        return json.JSONEncoder(sort_keys=True, indent=1).encode(self.data)

    def save(self,file):
        jsout = open(file,"w+")
        jsout.write(self.__str__())
        jsout.close()


