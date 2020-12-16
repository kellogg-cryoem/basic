import pandas as pd
import random
import argparse
import os

parser = argparse.ArgumentParser(description='This is a script used to split up particle sets from star files. Requires star file input from csparc2star.py')
parser.add_argument('-star_file',default='test.star',type=str,help='Provide the path to the star file to be split up')
parser.add_argument('-count',default=1,type=int,help='Provide the number of particle subsets to be created')
parser.add_argument('-output_header',default='TESTSTAR-SPLIT_',type=str,help='Provide the prefix to the output star files. Number ID and .star suffix will be added by script')
parser.add_argument('--random_seed',type=int,default=10,help='Provide integer value to serve as random seed.')
parser.add_argument('--header_temp',type = str,default = 'header.tmp', help = 'Name of temporary file to store header information.')
parser.add_argument('--particle_temp', type=str,default = 'particles.tmp', help = 'Name of temporary file to store particle information.')
parser.add_argument('--particle_count', default = 0, type = int, help = 'This can be used to override the count option and create particle splits of a specific size. Remainder will be put into another file.')
parser.add_argument('--outputonly',default = 0,type = int, help = 'Use this to specify how many particle sets should be outputted if you do not want all sets to be outputted.')
args = parser.parse_args()


random.seed(args.random_seed)

# This is used to collect the information from the file before the particle content, so things like star formatting and what not.
def GetHeader(File):
	Header = False
	Record = False
	headf = open(args.header_temp,'w')
	particlef = open(args.particle_temp,'w')
	for line in open(File):
		if(line[0] == '_'):
			Header = True
		if(Header == True and line[0] != '_'):
			Record = True
			headf.close()
			Header = False

		if(Record == True and line.strip() == ''):
			Record = False
			particlef.close()
			break

		if(Header == False and Record == False):
			headf.write(line)
		if(Header == True and Record == False):
			headf.write(line)
		if(Header == False and Record == True):
			particlef.write(line)

	if(Record == True):
		Record = False
		particlef.close()

# Used if particle_count isn't specified. Determines how many particles will show up in each file.
def CalculateParticles(File):
	particles = sum(1 for line in open(File))
	split = int(particles/args.count)
	remainder = int(particles%args.count)
	print('This is the amount of particles to expect within each split star file. Variance will occur based on remainder: ' + str(split))
	return split,remainder

# Used if particle_count is specified. Determines how many files will be needed.
def CalcFiles(File):
	particles = sum(1 for line in open(File))
	remainder = int(particles%args.particle_count)
	count = int(particles/args.particle_count)
	print('Due to use of the outputonly option, this is the revised amount of particles to expect within each split star file: ' + str(args.particle_count))
	return count,remainder

# Setting up header temp file
GetHeader(args.star_file)

# Particles is always calculated, however this can get modified is particle_count or outputonly is specified
ParticlePerFile,ParticleRemainder = CalculateParticles(args.particle_temp)
if(args.particle_count > 0):
	ChangeFileCount,ParticleSplitRemainder = CalcFiles(args.particle_temp)
	if(args.outputonly) > 0:
		ChangeFileCount = args.outputonly # No point in creating every file needed if only a subsection is needed for output
	args.count = ChangeFileCount

# Convenient way to separate lines of the particle temp file and to also be able to assign values to each
ParticleTable = pd.read_csv(args.particle_temp,header = None)

# NumberDict is used to record how many particles have already been assigned to a file
# PrefixList is a way of tracking which file prefixes to get rid of because they already reached the particle amount
# OGPref is just a copy of PrefixList that can be used later
NumberDict = {}
PrefixList = []
OGPref = []
print(args.count)
for i in range(0,args.count):
	NumberDict[i] = 0
	PrefixList.append(i)
	OGPref.append(i)

# Each line of particle table is iterated through and assigned to a random number
if(args.particle_count == 0):
	for i in range(0,len(ParticleTable)):
		CurrPref = PrefixList[random.randint(0,len(PrefixList)-1)]
		ParticleTable.at[i,'Prefix'] = CurrPref
		NumberDict[CurrPref] = NumberDict.get(CurrPref) + 1
		if(NumberDict.get(CurrPref) == ParticlePerFile):
			PrefixList.remove(CurrPref)
		if(len(PrefixList) == 0):
			break

	for x in range(i+1,len(ParticleTable)):
		RemList = OGPref
		CurrPref = RemList[random.randint(0,len(RemList)-1)]
		ParticleTable.at[i,'Prefix'] = CurrPref
		RemList.remove(CurrPref)
else: # If particle_count is specified, a random prefix is drawn along with a random particle, and they are assigned to one another. This prevents remainder bias.
	PossibleVals = list(range(0,len(ParticleTable)))
	Count = 0
	while(len(PrefixList) > 0):
		CurrPref = PrefixList[random.randint(0,len(PrefixList)-1)]
		CurrLine = PossibleVals[random.randint(0,len(PossibleVals))]
		ParticleTable.at[CurrLine,'Prefix'] = CurrPref
		NumberDict[CurrPref] = NumberDict.get(CurrPref) + 1
		PossibleVals.remove(CurrLine)
		if(NumberDict.get(CurrPref) == args.particle_count):
			PrefixList.remove(CurrPref)
	for i in range(0,len(PossibleVals)):  # -1 is assigned to all remainder particles so that they can be collected later
		ParticleTable.at[PossibleVals[i],'Prefix'] = -1
		if(len(PossibleVals) > 0):
			OGPref.append(-1)


if(args.outputonly == 0): # Basically, prints out every output star file possible, including the remainder
	for i in range(0,len(OGPref)):
		CurrStar = args.output_header + str(OGPref[i]) + '.star'
		if(OGPref[i] == -1):
			CurrStar = args.output_header + 'REMAINDER' + '.star'
		ParticleGrab = ParticleTable.loc[ParticleTable['Prefix'] == OGPref[i]].reset_index(drop = True)
		f = open(CurrStar,'w')
		for line in open(args.header_temp):
			f.write(line)
		for x in range(0,len(ParticleGrab)):
			f.write(ParticleGrab.at[x,0])
			f.write('\n')
		f.close()
else: # Used for when outputonly is specified
	if(-1 in OGPref):
		OGPref.remove(-1) # We don't want the remainder to be a possible output if you're using outputonly
	for i in range(0,args.outputonly):
		CurrPrefix = OGPref[random.randint(0,len(OGPref)-1)]
		CurrStar = args.output_header + str(OGPref[i]) + '.star'
		ParticleGrab = ParticleTable.loc[ParticleTable['Prefix'] == OGPref[i]].reset_index(drop = True)
		f = open(CurrStar,'w')
		for line in open(args.header_temp):
			f.write(line)
		for x in range(0,len(ParticleGrab)):
			f.write(ParticleGrab.at[x,0])
			f.write('\n')
		f.close()

# Removing temporary files created, header content file and particle content file
os.remove(args.header_temp)
os.remove(args.particle_temp)