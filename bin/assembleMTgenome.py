#!/usr/bin/python

"""
Written by Ernesto Picardi - e.picardi@biologia.uniba.it
"""

import getopt, sys, os, re

mt_track="""track db="hg18" type="bed" name="mtGenes" description="Annotation" visibility="3" itemRgb="On"
chrM 0 578 D-Loop 0 - 1 578 165,42,42
chrM 16024 16571 D-Loop 0 - 16024 16571 165,42,42
chrM 649 1603 RNR1 0 - 649 1603 255,140,0
chrM 1672 3230 RNR2 0 - 1672 3230 255,140,0
chrM 3307 4263 ND1 0 + 3307 4263 255,0,0
chrM 4470 5512 ND2 0 + 4470 5512 255,0,0
chrM 5904 7446 COX1 0 + 5904 7446 255,0,0
chrM 7586 8270 COX2 0 + 7586 8270 255,0,0
chrM 8527 9208 ATP6 0 + 8527 9208 255,0,0
chrM 8366 8573 ATP8 0 + 8366 8573 255,0,0
chrM 9207 9991 COX3 0 + 9207 9991 255,0,0
chrM 10059 10405 ND3 0 + 10059 10405 255,0,0
chrM 10470 10767 ND4L 0 + 10470 10767 255,0,0
chrM 10760 12138 ND4 0 + 10760 12138 255,0,0
chrM 12337 14149 ND5 0 + 12337 14149 255,0,0
chrM 14149 14674 ND6 0 - 14149 14674 0,0,255
chrM 14747 15888 CYTB 0 + 14747 15888 255,0,0
"""

mtt={}

def usage():
	print """Assembling MT-DNA from SAM/BAM/Pileup files
Version 1.0 - Ernesto Picardi - 2011-2012
Options:
	-f		Reference MT-DNA in fasta 
	-i 		Input File [.pileup .sam or .bam]
	-a		Human Genome in fasta [SAM or BAM only]
	-q		min per base quality score [25]
	-c		min. confidence level [0.75]
	-d		min. coverage depth [5]
	-g		min. gap lentgh [10]
	-o		output base name [mtDNAassembly]
	-s		Samtools executable [/usr/local/bin/samtools]
	-F		generate fasta output [no]
	-C		generate coverage file [no]
	-U		generate UCSC track file [no]
	-P		print out basic statistics [no]
	-N		normalize bed graph [no]
	-A		add a value to name field of UCSC track [None]
	-D		add a value to description field of UCSC track [None]
	
	"""

try:
	opts, args = getopt.getopt(sys.argv[1:], "hf:i:q:c:d:o:g:a:s:FCUPNA:D:")
except getopt.GetoptError, err:
	print str(err) 
	usage()
	sys.exit()
mtdnafile='/references/hg38.p13/GCF_000001405.39_GRCh38.p13_genomic.MT.fna'
inputfile=''
hgenome='/home/amellalli/analisis/003-sam-conver/GCF_000001405.39_GRCh38.p13_genomic.fasta.fai'
mqual=25
clev=0.75
cov=5
glen=10
basename='mtDNAassembly'
sexe='samtools'
crf=0
crc=0
cru=0
pout=0
normb=0
addv=''
addd=''

print(opts, args)

for o,a in opts:
	if o == "-h":
		usage()
		sys.exit()
	elif o == "-f": mtdnafile = a
	elif o == "-i": inputfile = a
	elif o == "-q": mqual = int(a)
	elif o == "-c": clev = float(a)
	elif o == "-d": cov = int(a)
	elif o == "-g": glen = int(a)
	elif o == "-a": hgenome = a
	elif o == "-o": basename = a
	elif o == "-s": sexe = a
	elif o == "-F": crf = 1
	elif o == "-C": crc = 1
	elif o == "-U": cru = 1	
	elif o == "-P": pout = 1		
	elif o == "-N": normb = 1
	elif o == "-A": addv = a
	elif o == "-D": addd = a	
	else:
		assert False, "Unhandled option."


if not os.path.exists(mtdnafile):
	usage()
	sys.exit('File %s does not exist.' %(mtdnafile))
if not os.path.exists(inputfile):
	usage()
	sys.exit('File %s does not exist.' %(inputfile))	
ext=inputfile.split('.')[-1]
basext=inputfile.replace('.'+ext,'')
if ext not in ['sam','bam','pileup']:
	usage()
	sys.exit('Input file name must contain: .sam, .bam or .pileup.')
samfile=basext+'.sam'
bamfile=basext+'.bam'
pileupfile=basext+'.pileup'
if ext in ['sam','bam'] and not os.path.exists(hgenome):
	sys.exit('Human genome file does not exist.')
if ext in ['sam','bam'] and not os.path.exists(hgenome+'.fai'):
	sys.exit('Human genome indices do not exist. Run samtools faidx fist.')

r=re.compile("#+")
r1=re.compile("""\^.{1}""")
rr=re.compile("[\+\-]{1}[0-9]+")

def normS(s,ref):
	c=re.finditer(rr,s)
	sl=list(s)
	cc=[(x.start(),x.end()) for x in c]
	for i in cc:
		n=int(''.join(sl[i[0]+1:i[1]]))
		sl[i[0]:i[1]+n]=['#' for xx in range(len(sl[i[0]:i[1]+n]))]
	ns=''.join(sl)
	ns=ns.replace('#','')
	ss=''
	for i in ns:
		if i in '.,ACGTNacgtN<>*': ss+=i
	return (ss.replace('.',ref)).replace(',',ref)
	
def nuc(seq):
	d={'A':0,'C':0,'G':0,'T':0,'N':0}
	for i in seq:
		if d.has_key(i): d[i]+=1
		else: d['N']+=1
	return d

dn={'A':'T','T':'A','C':'G','G':'C'}
def comp(s):
	ss=''
	for i in s:
		if dn.has_key(i): ss+=dn[i]
		else: ss+='N'
	return ss

def ff(v,l):
	for i in l:
		x=0
		for j in i:
			if j in v: x+=1
		if x==len(v): return i
	return 0

dIUPAC={'AG':'R','CT':'Y','GC':'S','AT':'W','GT':'K','AC':'M','CGT':'B','AGT':'D','ACT':'H','ACG':'V'}
def getIUPAC(f):
	vv=''.join([i[1] for i in f if i[0]>0])
	k=ff(vv,dIUPAC.keys())
	if k!=0: return dIUPAC[k]
	else: return '#'

def freq(d):
	f=[]
	for i in d:
		try: v=float(d[i])/sum(d.values())
		except: v=0.0
		f.append((v,i))
	f.sort()
	f.reverse()
	maxv=[f[0]]
	for i in f[1:]:
		if i[0]==maxv[0][0]: maxv.append(i)
	if len(maxv)==1:
		if maxv[0][0]>=clev: return maxv[0][1]
		else: return getIUPAC(f)
	elif len(maxv)>1: return getIUPAC(f)

if mtdnafile==None:
	usage()
	sys.exit('Please insert a valid mtDNA file in fasta format.')
if inputfile==None:
	usage()
	sys.exit('Please insert a valid pileup file.')


#if ext=='sam':
	#print 'Converting SAM to BAM...'
	#cmd='%s view -bt %s.fai %s > %s' %(sexe,hgenome,samfile,bamfile)
	#os.system(cmd)
	#ext='bam'
#if ext=='bam':
	#print 'Sorting and indexing BAM...'
	#cmd1='%s sort %s.bam %s-sorted' %(sexe,basext,basext)
	#cmd2='%s index %s-sorted.bam' %(sexe,basext)
	#os.system(cmd1)
	#os.system(cmd2)
	#print 'Creating pileup...'
	#cmd3='%s mpileup -f %s %s-sorted.bam > %s.pileup' %(sexe,hgenome,basext,basext)
	#os.system(cmd3)
	
mtdna={}
x=1
print 'Reading mtDNA sequence...'
f=open(mtdnafile)
for i in f:
	if i.strip()=='': continue
	if i.startswith('>'): continue
	for j in i.strip():
		mtdna[x]=(j.upper(),['#',(0,0,0,0),0,0.0])
		x+=1
f.close()

print 'Reading pileup file...'
f=open(pileupfile)
for i in f:
	if i.strip()=='': continue
	l=(i.strip()).split('\t')
	if l[0]!='chrM': continue
	pos=int(l[1])
	ref,seq,qual=l[2],normS(re.sub(r1,"",l[4]),l[2]),l[5]
	s,q='',0
	for j in range(len(seq)):
		if seq[j] not in '<>*' and ord(qual[j])-33 >= mqual:
			s+=seq[j].upper()
			q+=(ord(qual[j])-33)
	try: mq=float(q)/len(s)
	except: mq=0.0
	dnuc=nuc(s)
	mfreq=freq(dnuc)
	lnuc=(dnuc['A'],dnuc['C'],dnuc['G'],dnuc['T'])
	cnuc='#'
	if len(s) >= cov: cnuc=mfreq
	#print pos,cnuc,s,dnuc
	mtdna[pos][1][0]=cnuc
	mtdna[pos][1][1]=lnuc
	mtdna[pos][1][2]=len(s)
	mtdna[pos][1][3]=mq
f.close()



print 'Assembling...'
fastafile=basename+'-genome.fasta'
tablefile=basename+'-table.txt'
statfile=basename+'-statistics.txt'
coveragefile=basename+'-coverage.txt'
contigfile=basename+'-Contigs.fasta'
trackfile=basename+'-UCSCtrack.bed'
track=['browser position chrM\nbrowser hide all\n']
track.append(mt_track)
track.append('track db="hg18" type="bedGraph" name="Reads%s" description="Coverage%s" visibility="full" color=0,128,0\n' %(addv,addd))
aseq=''
f=open(tablefile,'w')
f.write('Position\tRefNuc\tConsNuc\tCov\tMeanQ\tBaseCount(A,C,G,T)\n')
assb,totb=0,0
cop=0
maxCval=1
for i in range(len(mtdna)):
	#print i+1, mtdna[i+1]
	line=[str(i+1),mtdna[i+1][0],mtdna[i+1][1][0],str(mtdna[i+1][1][2]),"%.2f" %(mtdna[i+1][1][3]),str(mtdna[i+1][1][1])]
	f.write('\t'.join(line)+'\n')
	aseq+=mtdna[i+1][1][0]
	totb+=1
	if mtdna[i+1][1][0] !='#':
		assb+=1
		cop+=mtdna[i+1][1][2]
		track.append('chrM %i %i %i\n' %(i,i+1,mtdna[i+1][1][2]))
		if mtdna[i+1][1][2] > maxCval: maxCval=mtdna[i+1][1][2]

try:passb=(float(assb)/totb)*100
except: passb=0.0
try:covmt=(float(cop)/assb)
except: covmt=0.0
fseq=aseq.replace('#','N')
rseq=comp(fseq)
bcomp=nuc(fseq)
bcomp2=nuc(rseq)
try:pa=float(bcomp['A'])/sum(bcomp.values())
except: pa=0.0
try:pc=float(bcomp['C'])/sum(bcomp.values())
except: pc=0.0
try:pg=float(bcomp['G'])/sum(bcomp.values())
except: pg=0.0
try:pt=float(bcomp['T'])/sum(bcomp.values())
except: pt=0.0
try: pgc=float(bcomp['G']+bcomp['C'])/sum(bcomp.values())
except: pgc=0.0
try:gcskl=float(bcomp['C']-bcomp['G'])/(bcomp['C']+bcomp['G'])
except: gcskl=0.0
try:gcskh=float(bcomp2['C']-bcomp2['G'])/(bcomp2['C']+bcomp2['G'])
except: gcskh=0.0
gaps=[]
for i in re.finditer(r,aseq):
	cc=(i.start()+1,i.end())
	if (cc[1]-cc[0])+1 >= glen: gaps.append(cc)
"""
f=open(statfile,'w')
f.write('Fraction of assembled bases: %.4f\n' %(passb))
f.write('Number of Gaps: %i\n' %(len(gaps)))
f.write('Base composition:\n')
f.write('Fraction of As: %.2f\n' %(pa))
f.write('Fraction of Cs: %.2f\n' %(pc))
f.write('Fraction of Gs: %.2f\n' %(pg))
f.write('Fraction of Ts: %.2f\n' %(pt))
f.write('G+C content: %.2f\n' %(pgc))
#f.write('GCskew L strand: %.2f\n' %(gcskl))
#f.write('GCskew H strand: %.2f\n' %(gcskh))
f.close()
"""
contigs=[]
if len(gaps)!=0:
	for i in range(len(gaps)-1):
		cc=(gaps[i][1]+1,gaps[i+1][0]-1)
		contigs.append((cc,fseq[cc[0]-1:cc[1]]))
	if gaps[0][0]!=1:
		cc=(1,gaps[0][0]-1)
		contigs.insert(0,(cc,fseq[cc[0]-1:cc[1]]))
	if gaps[-1][1]!=len(aseq):
		cc=(gaps[-1][1]+1,len(aseq))
		contigs.append((cc,fseq[cc[0]-1:cc[1]]))
	contigs.sort()
else:
	cc=(1,len(aseq))
	contigs=[(cc,fseq[cc[0]-1:cc[1]+1])]

# print out option
if pout:
	print 'Basic statistics:'
	print 'Assembled bases: %.2f' %(passb)+'%'
	print 'Mean coverage depth: %.2f' %(covmt)
	print 'Number of Contigs: %i' %(len(contigs))
	print 'Base composition [A,C,G,T]: %.2f,%.2f,%.2f,%.2f' %(pa,pc,pg,pt)
#

track.append('track db="hg18" type="bed" name="Assembly%s" description="Contigs%s" color=255,0,0 visibility="1"\n' %(addv,addd))
x=1
for i in contigs:
	track.append('chrM %i %i Contig.%i 0\n' %(i[0][0]-1,i[0][1],x))
	x+=1
if len(gaps)!=0:
	track.append('track db="hg18" type="bed" name="GAPS%s" description="Gaps%s" color=0,0,0 visibility="1"\n' %(addv,addd))
	x=1
	for i in gaps:
		track.append('chrM %i %i Gap.%i 0\n' %(i[0]-1,i[1],x))
		x+=1

if cru:
	f=open(trackfile,'w')
	for i in track:
		ll=i
		if normb:
			line=(ll.strip()).split(' ')
			if line[0].startswith('chrM') and len(line)==4:
				try:v=float(line[3])/maxCval
				except: v=0.0
				line[3]='%.3f' %(v)
				line=' '.join(line)+'\n'
				ll=line
		f.write(ll)
	f.close()

dass={}

if crf: f=open(contigfile,'w')
x=1
for i in contigs:
	if crf: f.write('>Contig.%i|%i-%i\n' %(x,i[0][0],i[0][1]))
	dass[i[0]]=[0,0,0,0,0]
	for j in range(0,len(i[1]),60):
		if crf: f.write(i[1][j:j+60]+'\n')
	x+=1
if crf: f.close()

#
dann={(1,578):['D-Loop1',0,0],(16025,16571):['D-Loop2',0,0],(650,1603):['RNR1',0,0],(1673,3230):['RNR2',0,0],(3308,4263):['ND1',0,0],(4470,5512):['ND2',0,0],(5904,7446):['COX1',0,0],(7586,8270):['COX2',0,0],(8527,9208):['ATP6',0,0],(8366,8573):['ATP8',0,0],(9207,9991):['COX3',0,0],(10059,10405):['ND3',0,0],(10470,10767):['ND4L',0,0],(10760,12138):['ND4',0,0],(12337,14149):['ND5',0,0],(14149,14674):['ND6',0,0],(14747,15888):['CYTB',0,0]}
#

for i in range(len(mtdna)):
	for j in dass:
		if j[0]<=i+1<=j[1]:
			if mtdna[i+1][1][2]>0:
				dass[j][1]+=1
				dass[j][2]+=mtdna[i+1][1][2]
			if mtdna[i+1][1][0] !='#':
				dass[j][3]+=1
				dass[j][4]+=mtdna[i+1][1][2]
				dass[j][0]+=1
	for j in dann:
		if j[0]<=i+1<=j[1]:
			if mtdna[i+1][1][0] !='#':
				dann[j][1]+=1
				dann[j][2]+=mtdna[i+1][1][2]			


if crc: f=open(coveragefile,'w')
x=1
if crc:f.write('Coverage of Assembled mtDNA. Coverage %.4f - Per base depth %.3f.\n' %(passb,covmt))
if crc:f.write('Name\tType\tStart\tEnd\tLength\tCoverage\tMeanDepth\n')
for i in contigs:
	vv=dass[i[0]]
	#cv1=float(vv[1])/vv[0]
	#cvd1=float(vv[2])/vv[0]
	#cv2=float(vv[3])/vv[0]
	cvd2=float(vv[4])/vv[0]
	contiglen=(i[0][1]-i[0][0])+1
	fcov=(float(vv[0])/contiglen)*100
	if crc:f.write('Contig.%i\tContig\t%i\t%i\t%i\t%.3f\t%.3f\n' %(x,i[0][0],i[0][1],contiglen,fcov,cvd2))
	x+=1
x=1
for i in gaps:
	gaplen=(i[1]-i[0])+1
	if crc:f.write('Gap.%i\tGap\t%i\t%i\t%i\t0.000\t0.000\n' %(x,i[0],i[1],gaplen))
	x+=1
for i in dann:
	vv=dann[i]
	flen=(i[1]-i[0])+1
	try: fcov=(float(vv[1])/flen)*100
	except: fcov=0.0
	try:cvd=float(vv[2])/vv[1]
	except: cvd=0.0
	if crc:f.write('%s\tAnnotation\t%i\t%i\t%i\t%.3f\t%.3f\n' %(vv[0],i[0],i[1],flen,fcov,cvd))

if crc:f.close()

if crf:f=open(fastafile,'w')
if crf:f.write('>mtDNAassembly\n')
for i in range(0,len(fseq),60):
	if crf:f.write(fseq[i:i+60]+'\n')
f.close()



