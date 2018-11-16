import sys

def init(printout, event=[]):
    #If using standard output to screen
    if printout == sys.stdout:
        return sys.stdout
    #If using a file object
    elif type(printout).__name__ == 'file':
        if printout.closed == True:
            return open(printout, 'a')
        else:
            return printout
    #If using a string filename
    elif type(printout).__name__ == 'str':
        if type(event).__name__ == 'instance':
            return open(event.modeldir + "/" + printout, 'a')
        if (type(event).__name__ == 'list') and (len(event) > 0):
            return open(event[0].modeldir + "/" + printout, 'a')
        else:
            return open(printout, 'a')
    else:
        return -1

def close(printout):
    #If using standard output to screen
    if printout == sys.stdout:
        return
    #If using a file object
    elif type(printout).__name__ == 'file':
        printout.close()
        return
    else:
        print("File was not closed.")
        return
