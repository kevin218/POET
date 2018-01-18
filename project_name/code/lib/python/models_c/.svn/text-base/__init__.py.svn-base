import os
import re
import sys

temp = str(__import__("models_c"))
dir  = re.findall('.+\'(.+)__init',temp)[0]

sys.path.append(dir+'ext_func/')
sys.path.append(dir+'py_func/')

ext = os.listdir(dir+'ext_func/')

for mod in ext:
    if mod.endswith('.so'):
        try:
            exec('from '+mod.partition('.')[0]+' import *')
        except:
            print('Warning: Could not import ' + mod)

pys = os.listdir(dir+'py_func/')
for mod in pys:
    if mod.endswith('.py'):
        try:
            exec('from '+mod.partition('.')[0]+' import *')
            #exec('reload(' + mod.partition('.')[0] + ')')
        except:
            print('Warning: Could not import ' + mod)

#sys.path.remove(dir+'ext_func/')
#sys.path.remove(dir+'ext_func/')
