import pandas as pd
import random
import argparse
import os

parser = argparse.ArgumentParser(description='This is a script used to split up particle sets from star files. Requires star file input from csparc2star.py')
parser.add_argument('-star_file',default='test.pdb',type=str,help='Provide the path to the star file to be split up')
parser.add_argument('-count',default=1,type=int,help='Provide the number of particle subsets to be created')
parser.add_argument('-output_header',default='TESTSTAR-SPLIT_',type=str,help='Provide the prefix to the output star files. Number ID and .star suffix will be added by script')
parser.add_argument('--random_seed',type=int,default=10,help='Provide integer value to serve as random seed.')
parser.add_argument('--header_temp',type = str,default = 'header.tmp', help = 'Name of temporary file to store header information.')
parser.add_argument('--particle_temp', type=str,default = 'particles.tmp', help = 'Name of temporary file to store particle information. ')
args = parser.parse_args()


StartRec = False
random.seed(10)

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

def CalculateParticles(File):
	particles = sum(1 for line in open(File))
	split = int(particles/args.count)
	remainder = int(particles%args.count)
	print('This is the amount of particles to expect within each split star file. Variance will occur based on remainder: ' + str(split))
	return split,remainder


GetHeader(args.star_file)
ParticlePerFile,ParticleRemainder = CalculateParticles(args.particle_temp)

ParticleTable = pd.read_csv(args.particle_temp,header = None)

NumberDict = {}
PrefixList = []
OGPref = []
for i in range(0,args.count):
	NumberDict[i] = 0
	PrefixList.append(i)
	OGPref.append(i)

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

for i in range(0,len(OGPref)):
	CurrStar = args.output_header + str(OGPref[i]) + '.star'
	ParticleGrab = ParticleTable.loc[ParticleTable['Prefix'] == OGPref[i]].reset_index(drop = True)
	f = open(CurrStar,'w')
	for line in open(args.header_temp):
		f.write(line)
	for x in range(0,len(ParticleGrab)):
		f.write(ParticleGrab.at[x,0])
		f.write('\n')
	f.close()

os.remove(args.header_temp)
os.remove(args.particle_temp)