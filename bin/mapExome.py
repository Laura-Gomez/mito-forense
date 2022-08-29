#!/usr/bin/python

"""
Written by Ernesto Picardi - e.picardi@biologia.uniba.it
"""

import getopt, sys, os

def usage():
	print """Map FASTQ onto mtDNA
Options:
	-a		Input Fastq 
	-b 		Input Fastq for pair-end (optional)
	-g		GSNAP excecutable [/usr/local/bin/gsnap]
	-D		GSNAP database location [/usr/local/share]
	-M		GSNAP database for mtDNA [mtdna]
	-H		GSNAP database for complete human genome [hg18]
	-t		GSNAP threads [8]
	-o		Out folder
	"""

try:
	opts, args = getopt.getopt(sys.argv[1:], "ha:b:g:D:M:H:t:o:")
except getopt.GetoptError, err:
	print str(err) 
	usage()
	sys.exit()

if len(opts)==0:
	usage()
	sys.exit()
fastq1=None
fastq2=None
gsnapexe='/usr/local/bin/gsnap'
gsnapdb='/references/gmapindex/hg38.p13/'
mtdb='chrM'
humandb='hg38'
mqual=30
thread=3
folder=os.path.join(os.getcwd(),'OUTfolder2')
for o,a in opts:
	if o == "-h":
		usage()
		sys.exit()
	elif o == "-a": fastq1 = a
	elif o == "-b": fastq2 = a
	elif o == "-g": gsnapexe = a
	elif o == "-D": gsnapdb = a
	elif o == "-M": mtdb = a
	elif o == "-H": humandb = a
	#elif o == "-q": mqual = int(a)
	elif o == "-t": thread = int(a)	
	elif o == "-o": folder = a		
	else:
		assert False, "unhandled option"

def rev(seq):
	d={'A':'T','T':'A','C':'G','G':'C','N':'N'}
	s=''.join([d[x] for x in seq])
	return s[::-1]

if not os.path.exists(folder): os.mkdir(folder)
        
if fastq2!=None:
	map1cmd='%s --gunzip -D %s -d %s -A sam --no-sam-header --nofails --pairmax-dna=500 --query-unk-mismatch=1 -n 1 -Q -O -t %i %s %s > %s 2> %s' %(gsnapexe,gsnapdb,mtdb,thread,fastq1,fastq2,os.path.join(folder,'outmt.sam'),os.path.join(folder,'logmt.txt'))
else:
	map1cmd='%s --gunzip -D %s -d %s -A sam --no-sam-header --nofails --pairmax-dna=500 --query-unk-mismatch=1 -n 1 -Q -O -t %i %s > %s 2> %s' %(gsnapexe,gsnapdb,mtdb,thread,fastq1,os.path.join(folder,'outmt.sam'),os.path.join(folder,'logmt.txt'))
print 'Mapping onto mtDNA...'
#print (map1cmd)
#exit() 
os.system(map1cmd)
#

print 'Extracting FASTQ from SAM...'
mtoutsam=os.path.join(folder,'outmt.sam')
dics={}
f=open(mtoutsam)
for i in f:
	if i.strip()=='': continue
	l=(i.strip()).split('\t')
	if l[2]=='*': continue
	if dics.has_key(l[0]): dics[l[0]].append(l)
	else: dics[l[0]]=[l]
f.close()
single,pair1,pair2=[],[],[]

for i in dics:
	ll=dics[i]
	if len(ll)==1:
		strand,seq,qual=int(ll[0][1]) & 16,ll[0][9],ll[0][10]
		if strand==16: seq,qual=rev(seq),qual[::-1]
		entry='\n'.join(['@'+ll[0][0],seq,'+',qual])+'\n'
		single.append(entry)
	else:
		strand,seq,qual=int(ll[0][1]) & 16,ll[0][9],ll[0][10]
		if strand==16: seq,qual=rev(seq),qual[::-1]
		entry='\n'.join(['@'+ll[0][0],seq,'+',qual])+'\n'
		pair1.append(entry)
		strand,seq,qual=int(ll[1][1]) & 16,ll[1][9],ll[1][10]
		if strand==16: seq,qual=rev(seq),qual[::-1]
		entry='\n'.join(['@'+ll[1][0],seq,'+',qual])+'\n'		
		pair2.append(entry)

sig,pai=0,0
if len(single)!=0:
	mtoutfastq=os.path.join(folder,'outmt.fastq')
	out=open(mtoutfastq,'w')
	out.writelines(single)
	out.close()
	sig=1
if len(pair1)!=0:
	mtoutfastq1=os.path.join(folder,'outmt1.fastq')
	out=open(mtoutfastq1,'w')
	out.writelines(pair1)
	out.close()	
	mtoutfastq2=os.path.join(folder,'outmt2.fastq')
	out=open(mtoutfastq2,'w')
	out.writelines(pair2)
	out.close()
	pai=1


if sig:
	print 'Mapping onto complete human genome...single reads'
	map2cmd='%s -D %s -d %s -A sam --no-sam-header --nofails --query-unk-mismatch=1 -O -t %i %s > %s 2> %s' %(gsnapexe,gsnapdb,humandb,thread,mtoutfastq,os.path.join(folder,'outhumanS.sam'),os.path.join(folder,'loghumanS.txt'))
	os.system(map2cmd)
if pai:
	print 'Mapping onto complete human genome...pair reads'
	map3cmd='%s -D %s -d %s -A sam --no-sam-header --nofails --query-unk-mismatch=1 -O -t %i %s %s > %s 2> %s' %(gsnapexe,gsnapdb,humandb,thread,mtoutfastq1,mtoutfastq2,os.path.join(folder,'outhumanP.sam'),os.path.join(folder,'loghumanP.txt'))
	os.system(map3cmd)	

print 'Reading Results...'
if sig:
	hgoutsam=os.path.join(folder,'outhumanS.sam')
	dicsingle={}
	f=open(hgoutsam)
	for i in f:
		if i.strip()=='': continue
		l=(i.strip()).split('\t')
		if l[2]=='*': continue
		if dicsingle.has_key(l[0]):
			dicsingle[l[0]].append(l)
		else:
			dicsingle[l[0]]=[l]
	f.close()
if pai:
	hgoutsam2=os.path.join(folder,'outhumanP.sam')
	dicpair={}
	f=open(hgoutsam2)
	for i in f:
		if i.strip()=='': continue
		l=(i.strip()).split('\t')
		if l[2]=='*': continue
		if dicpair.has_key(l[0]):
			dicpair[l[0]].append(l)
		else:
			dicpair[l[0]]=[l]
	f.close()
	
print 'Filtering reads...'
good=[]
for i in dics:
	ll=dics[i]
	if len(ll)==1:
		if dicsingle.has_key(i):
			r=dicsingle[i]
			#print ll
			#print r
			if len(r)==1:
				if r[0][2]==ll[0][2] and ll[0][3]==r[0][3]: good.append('\t'.join(ll[0])+'\n')
		else: good.append('\t'.join(ll[0])+'\n')
	else:
		if dicpair.has_key(i):
			r=dicpair[i]
			if len(r) == 2:
				if r[0][2]==ll[0][2] and ll[0][3]==r[0][3] and r[1][2]==ll[1][2] and ll[1][3]==r[1][3]:
					good.append('\t'.join(ll[0])+'\n')
					good.append('\t'.join(ll[1])+'\n')
		else:
			good.append('\t'.join(ll[0])+'\n')
			good.append('\t'.join(ll[1])+'\n')

"""
hgoutsam=os.path.join(folder,'outhuman.sam')
goodhits={}
f=open(hgoutsam)
for i in f:
	if i.strip()=='': continue
	l=(i.strip()).split('\t')
	nn=l[0].split()[0]
	if l[2]=='chrM': mt=1
	else: mt=2
	if goodhits.has_key(nn): goodhits[nn].append(mt)
	else: goodhits[nn]=[mt]
f.close()
"""
finalsam=os.path.join(folder,'OUT.sam')
out=open(finalsam,'w')
out.writelines(good)
out.close()

"""
f=open(mtoutsam)
xx=1
for i in f:
	if i.strip()=='': continue
	l=(i.strip()).split('\t')
	nn=l[0].split()[0]+'.'+str(xx)
	if not goodhits.has_key(nn): out.write(i)
	else:
		if sum(goodhits[nn])==1: out.write(i)
	xx+=1
f.close()
out.close()
"""
print 'Outfile saved on %s.' %(finalsam)
print 'Done.'

