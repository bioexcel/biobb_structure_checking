# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.
class GraphmlWriter():
    def __init__(self):
        self.header =  """<?xml version="1.0" encoding="UTF-8"?>
    <graphml xmlns="http://graphml.graphdrawing.org/xmlns"
        xmlns:java="http://www.yworks.com/xml/yfiles-common/1.0/java"
        xmlns:sys="http://www.yworks.com/xml/yfiles-common/markup/primitives/2.0"
        xmlns:x="http://www.yworks.com/xml/yfiles-common/markup/2.0"
        xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
        xmlns:y="http://www.yworks.com/xml/graphml"
        xmlns:yed="http://www.yworks.com/xml/yed/3"
        xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns http://www.yworks.com/xml/schema/graphml/1.1/ygraphml.xsd">
        <key for="node" id="d6" yfiles.type="nodegraphics"/>
        <key for="edge" id="d10" yfiles.type="edgegraphics"/>
        <graph id="G" edgedefault="directed">
"""
        self.footer = """
        </graph>
    </graphml>
"""
        self.body=[]

    def addResidue(self,i,r):
        chcolors = ["#FF0000","#00FF00","#00FF00","#FFFF00","#FF00FF","#00FFFF"]
        self.body.append("\
            <node id=\"" + r.resid(1) + "\"><data key=\"d6\"><y:ShapeNode><y:Fill color=\""\
            + chcolors[i] + "\" transparent=\"false\"/><y:NodeLabel>"+\
            r.resid(1) + "</y:NodeLabel><y:Shape type=\"ellipse\"/></y:ShapeNode></data></node>"
            )

    def addBond(self,pref,r1,r2):
        self.body.append("\
           <edge id=\"" + pref + r1.resid(1)+ r2.resid(1) + "\" source=\"" + r1.resid(1) + "\" target=\"" + r2.resid(1) + "\"></edge>"
            )
    def __str__(self):
        return self.header + '\n'.join(self.body) +"\n"+ self.footer

    def save(self,file):
        gmlout = open(file,"w+")
        gmlout.write(str(self))
        gmlout.close()


