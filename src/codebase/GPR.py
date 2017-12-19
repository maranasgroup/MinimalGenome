from pyparsing import *
import types
import pdb

class Node(object):
    def __init__(self, data):
        self.data = data
        self.children = []

    def add_child(self, obj):
        self.children.append(obj)
    
    def other_name(self, level=0):
        print '\t' * level + repr(self.data)
        for child in self.children:
            child.other_name(level+1)

# traverse GPR list and create trees
def traverse_GPR(GPR_list,root):
    if "AND" not in GPR_list and "OR" not in GPR_list:
        if type(GPR_list) is list: 
            root.add_child(Node(GPR_list[0]))
        else:
            root.add_child(Node(GPR_list))
        return

    if GPR_list != "AND" and GPR_list != "OR":    
        if "AND" in GPR_list:
            logic = "AND"
        if "OR" in GPR_list:
            logic = "OR"
        new_node = Node(logic)        
        root.add_child(new_node)
        root = new_node

    for child in GPR_list:
        # only go deeper is a list
        if child == "AND" or child == "OR":
            continue
        traverse_GPR(child,root)
    
def parathesis_to_list(GPR):
    #pyparse to break GPR into nested list
    enclosed = Forward()
    nestedParens = nestedExpr('(', ')', content=enclosed) 
    enclosed << (Word(alphanums+'.') | ',' | nestedParens)
    GPR_list = enclosed.parseString(GPR).asList()
    return GPR_list    

if __name__ == '__main__':
    # test examples   
    GPR = "((g1 OR g2) AND (g2 OR (g1 AND (g3 OR g4)) OR (g1 AND g4) OR (g3 AND g4)))"
    # GPR = "(((g1) OR (g2)) AND ((g2) OR ((g1) AND ((g3) OR (g4))) OR ((g1) AND (g4)) OR ((g3) AND (g4))))"
    # GPR = "((g1) OR (g2))"
    # GPR = "(g1 OR g2)"
    # GPR = "((g1 and g2) or (g2 and g3))"

    GPR = GPR.replace("or","OR").replace("and","AND")
    GPR_list = parathesis_to_list(GPR)
    print GPR_list

    # convert to trees
    GPR_tree = Node("R1")
    traverse_GPR(GPR_list[0],GPR_tree)
    GPR_tree.other_name()

    
