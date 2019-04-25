import numpy as np
import string

#	model name	model headers 			array dimensions
models = [["mandelecl", "midpt\twidth\tdepth\tt12\tt34\tflux\n", np.zeros((4, 6))],
		["risingexp", "goal\tm\tt0\n", np.zeros((4, 3))],
		["fallingexp", "goal\tm\tt0\n", np.zeros((4, 3))],
		["quadramp", "a\tb\tc\tx0\t\n", np.zeros((4, 4))],
		["linramp", "a\tb\tx0\t\n", np.zeros((4, 3))],
		["logramp", "t0\ta\tb\tc\td\te\n", np.zeros((4, 6))],
		["llramp", "x0\ta\tb\tc\n", np.zeros((4, 4))], 
		["sindecay", "x0\ta\tb\tc\td\n", np.zeros((4, 5))],
		["quadip", "a\tb\tc\td\te\tf\n", np.zeros((4, 6))],
		["quadip4", "a\tb\tc\td\te\tf\n", np.zeros((4, 24))],
		["cubicip", "a\tb\tc\td\te\tf\tg\th\ti\tj\n", np.zeros((4, 10))]]
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
eg. To set t34 = t12, set stepsize[4] = -3
NOTE2: Zero stepsize results in a fixed value.
****MODELS BELOW****
'''
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
	global models
	try:
		f = open(filename)
	except:
		print(filename, " does not exist.  Creating new file...")
		defaultfile = filename[:filename.rfind('/')+1] + event.eventname[:2] + '-params.txt'
		models = read_parameters(defaultfile) #Read from default/backup file
		write_parameters(filename)
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
	global models
	try:	
		f = open(filename)
	except:
		print(filename, " does not exist.  Creating new file...")
		models = read_parameters('params2.txt') #Read from default/backup file
		write_parameters(filename)
	
	print('Updating the following models:')
	for i in models_used:
		k = 0
		print(i[0])
		for j in models:
			if i[0] == j[0]: #Test if model names are equal
				models[k] = i
			k+=1
	write_parameters(filename)
	return

def read_parameters(file_name):
	f = open(file_name, "r")
	global comment
	comment = ""
	while True: #Read everything in the comment space.
		line = str(f.readline())
		if line == "****MODELS BELOW****\n":
			comment += line 
			break
		else:
			comment += line
	offset = False
	for i in range(0, len(models)):	
		if offset == False:
			modelname=str(f.readline())
		if modelname.endswith('\n'):
			modelname=modelname[:len(modelname)-1]
		#print(models[i][0], modelname)
		if models[i][0] == modelname:
			#Check
			models[i][0] = modelname
			models[i][1] = str(f.readline())
			#foo = str(f.readline())
			for j in range(0, 4):
				k=0
				for element in f.readline().split('\t'):
					models[i][2][j][k]=element				
					k+=1
			offset = False
		else:
			models[i] = models[i]
			if offset == False:
				offset = True
	f.close()
	print('Parameters have been successfully read from ' + file_name)
	return models

def write_parameters(file_name):
	f = open(file_name, "w")
	f.write(comment)
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

