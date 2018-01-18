#Written by Eric Tollerud (erik.tollerud@gmail.com)
#Modification History
#   2010-11-13  kevin   Modified to work with esp01 pipeline

def unpickle_old_pyfits(fn):
    """
    This function unpickles everything in the specified filename, and
    correctly adapts pyfits files pre 2.3 to the new module structure.
    """
    
    import cPickle
    from cPickle import Unpickler
    import imp,sys

    def fg(modname,classname):        
        if 'NP_pyfits' in modname:
            modname = modname.replace('NP_pyfits','core')
            
        if '.' in modname:
            mod = __import__(modname,fromlist=[0])
        else:
            mod = __import__(modname)
            
        return getattr(mod,classname)
    
    objs = []
    if isinstance(fn, str):
        f = open(fn)
    else:
        f = fn
    
    try:
        u = Unpickler(f)
        u.find_global = fg
        
        while True:
            objs.append(u.load())
            
    except EOFError:
        pass
    finally:
        if isinstance(fn, str):
            f.close()
    if len(objs) == 1:
        return objs[0]
    else:
        return objs
    
def pickle_fn(objlist,fn):
    import cPickle
    
    print objlist,fn
    f = open(fn,'w')
    try:
        for o in objlist:
            cPickle.dump(o,f,protocol=cPickle.HIGHEST_PROTOCOL)
    finally:
        f.close()
    
if __name__ == '__main__':
    import sys
    
    if len(sys.argv)<2:
        print 'Need to provide files to convert'
        sys.exit(1)
        
    fns = sys.argv[1:]
    
    for fn in fns:
        olist = unpickle_old_pyfits(fn)
        pickle_fn(olist,fn+'.new')
    
    
