import numpy as np
import string
import os

#	model name	model headers 			array dimensions
def init_models():
    models = [#ECLIPSES / TRANSITS
            ["mandelecl",    "midpt     \twidth     \tdepth     \tt12       \tt34       \tflux\n", np.zeros((4, 6))],
            ["mandelecl2",   "midpt2    \twidth2    \tdepth2    \tt122      \tt342      \tflux2\n", np.zeros((4, 6))],
            ["mandelecl3",   "midpt3    \twidth3    \tdepth3    \tt123      \tt343      \tflux3\n", np.zeros((4, 6))],
            ["mandeltr",     "trmidpt   \ttrrprs    \ttrcosi    \ttrars     \ttrflux    \ttrper\n", np.zeros((4, 6))],
            ["mandeltr2",    "trmidpt   \ttrrprs2   \ttrcosi2   \ttrars2    \ttrflux2   \ttrper2\n", np.zeros((4, 6))],
            ["trquad",       "trqmid    \ttrqrprs   \ttrqcosi   \ttrqars    \ttrqf      \ttrqp      \ttrqc1     \ttrqc2\n", np.zeros((4, 8))],
            ["trquad2",      "trq2mid   \ttrq2rprs  \ttrq2cosi  \ttrq2ars   \ttrq2f     \ttr2qp     \ttrq2c1    \ttrq2c2\n", np.zeros((4, 8))],
            ["trnlldsp",     "trspmid   \trprs      \tcosi      \tars       \ttrspf     \ttrspp     \ttrspc1    \ttrspc2    \ttrspc3    \ttrspc4\n", np.zeros((4, 10))],
            ["trnlldsp2",    "trspmid2  \trprs2     \tcosi2     \tars2      \ttrspf2    \ttrspp2    \ttrspc12   \ttrspc22   \ttrspc32   \ttrspc42\n", np.zeros((4, 10))],
            ["mandelgeom",   "midpt\twidth\tRp/Rs\tb\tflux\n", np.zeros((4, 5))],
	    ["mandelorbit",  "e\tomega\ti\tperiod\trplanet\trstar\tmstar\tecldepth\tflux\n", np.zeros((4, 9))],
            ["batman_trquad","t0        \trprs      \tperiod      \tars       \tcosi       \tecc       \tomega     \tu1     \tu2\n", np.zeros((4,9))],
            ["batman_ecl",   "eclmidpt  \tfpfs      \trprs      \tperiod      \tars       \tcosi       \tecc       \tomega\n", np.zeros((4,8))],
            #RAMPS
            ["constant",      "c\n", np.zeros((4, 1))],
            ["hook",         "hgoal     \thr0       \thr1       \thpm\n", np.zeros((4, 4))],
            ["hook2",        "h2goal    \th2r0      \th2r1      \th2pm      \th2per\n", np.zeros((4, 5))],
            ["heqramp",      "heqt0     \theqr0     \theqr1     \theqr2     \theqr3     \theqpm     \theqper\n", np.zeros((4, 7))],
            ["seramp",       "segoal    \tser0      \tser1      \tsepm\n", np.zeros((4, 4))],
            ["selramp",      "selgoal   \tselr0     \tselr1     \tselr2     \tselt0     \tselpm\n", np.zeros((4, 6))],
            ["seqramp",      "seqgoal   \tseqr0     \tseqr1     \tseqr2     \tseqr3     \tseqt0     \tseqpm\n", np.zeros((4, 7))],
            ["se2ramp",      "se2goal   \tse2r0     \tse2r1     \tse2pm0    \tse2r4     \tse2r5     \tse2pm1\n", np.zeros((4, 7))],
            ["llramp",       "llt0      \tllr6      \tllr2      \tllc       \tllt1\n", np.zeros((4, 5))], 
            ["lqramp",       "lqt0      \tlqr6      \tlqr7      \tlqr2      \tlqc       \tlqt1\n", np.zeros((4, 6))], 
            ["logramp",      "logt0     \tlogr9     \tlogr8     \tlogr7     \tlogr6     \tlogc\n", np.zeros((4, 6))],
            ["log4qramp",    "l4qt0     \tl4qr9     \tl4qr8     \tl4qr7     \tl4qr6     \tl4qr3     \tl4qr2     \tl4qc      \tl4qt1\n", np.zeros((4, 9))],
            ["linramp",      "linr2     \tlinc      \tlint0\n", np.zeros((4, 3))],
            ["quadramp",     "qrr3      \tqrr2      \tqrc       \tqrt0\n", np.zeros((4, 4))],
            ["risingexp",    "goal      \tm         \tt0\n", np.zeros((4, 3))],
            ["relramp",      "goal      \tm         \tt0        \ta         \tb         \tt1\n", np.zeros((4, 6))],
            ["reqramp",      "goal      \tm         \tt0        \ta         \tb         \tc         \tt1\n", np.zeros((4, 7))],
            ["re2ramp",      "goal      \ta         \tm1        \tt1        \tb         \tm2        \tt2\n", np.zeros((4, 7))],
            ["fallingexp",   "goal      \tm         \tt0        \n", np.zeros((4, 3))],
            ["felramp",      "goal      \tm         \tt0        \ta         \tt1\n", np.zeros((4, 5))],
            #SINUSOIDAL
            ["sindecay",     "x0        \ta         \tb         \tc         \td\n", np.zeros((4, 5))],
            ["sincos",       "a         \tp1        \tt1        \tb         \tp2        \tt2        \tc\n", np.zeros((4, 7))],
            ["sincos2",      "c1a       \tc1o       \tc2a       \tc2o       \ts1a       \ts1o       \ts2a       \ts2o       \tp         \tc         \tmidpt     \tt14       \tt12\n", np.zeros((4,13))],
            ["cosine8",      "c8a1      \tc8p1      \tc8t1      \tc8a2      \tc8p2      \tc8t2      \tc8a3      \tc8p3      \tc8t3      \tc8a4      \tc8p4      \tc8t4      \tc8a5      \tc8p5      \tc8t5      \tc8a6      \tc8p6      \tc8t6      \tc8a7      \tc8p7      \tc8t7      \tc8a8      \tc8p8      \tc8t8      \tc8c\n", np.zeros((4, 25))],
            #GAUSSIAN PROCESS
            ["gp_exp2",      "amp       \tscale     \tnsamp\n", np.zeros((4, 3))],
            #INSTRUMENTAL
            ["rotation",     "rota      \trotb      \troto\n", np.zeros((4, 3))],
            ["humidity",     "rha       \trhb\n", np.zeros((4, 2))],
            #["evenodd",      "eoa\n", np.zeros((4, 1))],
            #INTRAPIXEL / INTERPOLATION
            ["linip", "y0\tx0\ty1\tx1\ty2\tx2\ty3\tx3\ty4\tx4\ty5\tx5\ty6\tx6\ty7\tx7\ty8\tx8\n", np.zeros((4, 18))],
            ["quadip", "a\tb\tc\td\te\tf\n", np.zeros((4, 6))],
            ["quadip4", "a\tb\tc\td\te\tf\n", np.zeros((4, 24))],
            ["cubicip", "a\tb\tc\td\te\tf\tg\th\ti\tj\n", np.zeros((4, 10))],
            ["sexticip", "y6\tx6\ty5\tx5\ty4\tx4\ty3\tx3\ty2\tx2\ty1\tx1\tc\n", np.zeros((4, 13))],
            ["sexticipc", "y6\tx6\ty5\tx5\ty4\tx4\ty3\tx3\ty2x\tx2y\ty2\tx2\txy\ty1\tx1\tc\n", np.zeros((4, 16))],["posfluxlinip", "p0\tp1\tp2\tp3\tp4\tp5\tp6\tp7\tp8\ty0\tx0\ty1\tx1\ty2\tx2\ty3\tx3\ty4\tx4\ty5\tx5\ty6\tx6\ty7\tx7\ty8\tx8\n", np.zeros((4, 27))],
            ["cubicgw", "x1\tx2\tx3\ty1\ty2\ty3\tc\n", np.zeros((4, 7))],
            ["ballardip", "sigmay\tsigmax\tnbins\n", np.zeros((4, 3))],
            ["medianip", "rad\n", np.zeros((4, 1))],
            ["nnint", "minpts\n", np.zeros((4, 1))],
            ["bilinint", "minpts\n", np.zeros((4, 1))],
            ["pnnint", "minpts\n", np.zeros((4, 1))],
            ["pbilinint", "minpts\n", np.zeros((4, 1))],
            ["ipspline", "yknots\txknots\n", np.zeros((4, 2))],
            ["unispline", "usnk\tusk\n", np.zeros((4, 2))],
            #OTHERS
            ["ortho",        "none\n", np.zeros((4, 1))],
            ["rednoise",     "whinoise  \trednoise  \tgamma\n", np.zeros((4, 3))],
            ["posflux", "p0\tp1\tp2\tp3\tp4\tp5\tp6\tp7\tp8\n", np.zeros((4, 9))],
            ["posflux2", "p0\tp1\tp2\tp3\tp4\tp5\tp6\tp7\tp8\n", np.zeros((4, 9))],
            ["vsll", "vsx0\tvsa\tvsb\tvsc\tvsx1\n", np.zeros((4, 5))],
            ["vsspline", "k0\tk1\tk2\tk3\tk4\tk5\tk6\tk7\tk8\tk9\tk10\tk11\n", np.zeros((4, 12))],
            ["flatfield3", "f0\tf1\tf2\tf3\tf4\tf5\tf6\tf7\tf8\tflux\n", np.zeros((4, 10))],
            ["expramp", "goal\tm\ta\n", np.zeros((4, 3))],
            ["not0risingexp", "goal\tm\ta\tb\tc\n", np.zeros((4,5))]]  # FINDME ccampo added 2010-08-20
    return models

def init_comment():
    comment = '''****COMMENT SPACE****
This file lists the models and their parameters.
The parameters are presented in a table headed
by the name of each model. Column headers are 
unique to each model, but rows follow a 
standard format:
	Row 1: Initial parameters
	Row 2: Lower bounds
	Row 3: Upper bounds
	Row 4: Step size
NOTE1: To set one parameter equal to another, set its stepsize to the
negative value of the location of the paramter you wish it to be equal to.
CAUTION: The first location is 1, NOT 0, since you can't have -0.
    eg. To set t12 = t34, set stepsize[3] = -5
NOTE2: Zero stepsize results in a fixed value.
****MODELS BELOW****
'''
    return comment
#Above is the default comment for the output text file.  The actual comment field within the file can be edited.
'''
To add a new model to this array, simply folow the folowing format:
	["mymodelid", "myheading1\tmyheading2\tmyheading3\tmyendheading\n", np.zeros((4, mynumber_of_parameters))]
This creates an array of zeros which should look like:
modelid
myheading1	myheading2	myheading3	myendheading
0		0		0		0
0		0		0		0
0		0		0		0
0		0		0		0
Call the function write_array("myparams.txt", models) to add this change to the parameter file.
Open that file and edit the parameters. 
modelid
heading1	heading2	heading3	endheading
0.1234		0.456		7.89		0.0
....
....
....

Understand the purpose of each row:
	Row 1: Initial parameters
	Row 2: Lower bounds
	Row 3: Upper bounds
	Row 4: Step size

To access any particular models in this list, say the first one, use this format:
models[0] returns the first model on the list.
models[0][0] returns the model title of the first model.
models[0][1] returns the column headers of the first model.
models[0][2] returns the parameters of this model.
models[0][2][0] returns the first row of this model's parameters.
models[0][2][0][0] returns the first value of the first row of this models's parameters.

 '''

def read(filename, selected_models, event):
    #global models
    if os.path.isfile(filename) == False:
        print(filename, " does not exist.  Creating new file...")
        defaultfile = filename[:filename.rfind('/')+1] + event.eventname[:2] + '-params.txt'
        models = read_parameters(defaultfile) #Read from default/backup file
        write_parameters(filename, models)
    models = read_parameters(filename)
    models_returned = [0]*len(selected_models)
    for i in models:
        k = 0
        for j in selected_models:
            if i[0] == j:
                models_returned[k] = i
            k+=1
    return models_returned #Equal to models_used in the next function.

def write(filename, models_used):
    #global models
    if os.path.isfile(filename) == False:
        print(filename, " does not exist.  Creating new file...")
        models = read_parameters('params2.txt') #Read from default/backup file
        write_parameters(filename, models)
    else:
        models = read_parameters(filename)
    
    print('Updating the following models:')
    for i in models_used:
        k = 0
        print(i[0])
        for j in models:
            if i[0] == j[0]: #Test if model names are equal
                models[k] = i
            k+=1
    write_parameters(filename, models)
    return

def read_parameters(file_name):
    models = init_models()
    f = open(file_name, "r")
    #global comment
    comment = ""
    while True: #Read everything in the comment space.
        line = str(f.readline())
        if line == "****MODELS BELOW****\n":
            comment += line 
            break
        else:
            comment += line
    while True:
        modelname = f.readline()
        # Test for EOF
        if len(modelname) == 0:
            break
        modelname = modelname.rstrip('\n')
        # Look for matching model names
        for i in range(0, len(models)):
            if models[i][0] == modelname:
                # Read variable names
                foo2 = f.readline()
                for j in range(0, 4):
                    k=0
                    for element in f.readline().split('\t'):
                        models[i][2][j][k]=element				
                        k+=1
    '''
    offset = False
    for i in range(0, len(models)):	
        if offset == False:
            modelname=str(f.readline())
        if modelname.endswith('\n'):
            modelname=modelname[:len(modelname)-1]
        #print(models[i][0], modelname)
        if models[i][0] == modelname:
            #Check
            #models[i][0] = modelname
            #print(models[i][1])
            models[i][1] = str(f.readline())
            #print(models[i][1])
            #foo2 = str(f.readline())
            for j in range(0, 4):
                k=0
                for element in f.readline().split('\t'):
                    models[i][2][j][k]=element				
                    k+=1
                offset = False
        else:
            offset = True
    '''
    f.close()
    print('Parameters have been successfully read from ' + file_name)
    return models

def write_parameters(file_name, models):
	f = open(file_name, "w")
	f.write(init_comment())
	for i in range (0, len(models)):
		f.write(models[i][0])
		if models[i][0].endswith('\n') == False:
			f.write('\n')
		f.write(models[i][1])
		for j in range(0, 4):
			for k in range(0, len(models[i][2][j])):
				f.write(str('%7.4e' % models[i][2][j][k]))
				if k < len(models[i][2][j]) - 1:
					f.write('\t')
			f.write('\n')
	f.close()
	print('Parameters have been successfully written to ' + file_name)
	return

