import subprocess, sys, re

inputfile = open(sys.argv[1])
discarded = open('discarded.tsv', 'w')
voted = open('voted.tsv', 'w')
result = open('results.tsv', 'w')
forwardaa = 0
reverseendtt = 1
errors = []
readinput = ""
reads = []
counter = 0

for line in inputfile.readlines():
	
	if line.startswith('Chr') or len(line) < 20:
		continue
	bedline = line.split('\t')
	
	if len(bedline) < 4:
		continue
	if line.startswith("Positive") or line.startswith("Sample") or line.startswith("Table"):
		print line.strip()
		continue
	counter+=1
	if counter%100 == 0:
		print str(counter) +" "
	
	sample = str(bedline[4].split('.')[0])
	
	cramfile =  subprocess.check_output('ls crams/*' +sample +'_*.*am', shell=True).strip()
	start = int(bedline[1])
	end = int(bedline[2])
	middle = start + (end-start)
	newstart = middle - 200
	newend = middle + 200
	chr = bedline[0]
	query = chr +":" +str(newstart) +"-" +str(newend)
	readinput = subprocess.check_output('samtools view ' +cramfile +' ' +query, shell=True)
	reads = readinput.split('\n')
	forwardReads = []
	forwardaa = []
	reversett = []
	aaarepeats = []
	tttrepeats = []
	
	for readline in reads:
		read = readline.split('\t')
		if len(read) < 2:
			continue		
		if read[5] == '*':
			continue
		dir = int(read[1])%pow(2,5) < pow(2,4)
		readstart = int(read[3])-1
		seq = str(read[9])		
		
		if 'S' in read[5]:
			match = re.findall(r'(\d+)S', read[5])
			if len(match) == 1:
				if read[5].endswith('S'):
					readend += 0		
				else:
					readstart -= int(match[0])
			else:
				first = True
				for nro in match:
					if first:
						readstart -= int(nro)
					else:
						readend += 0
					first = False

		readend = readstart + len(seq)
		if 'D' in read[5]:
			match = re.findall(r'(\d+)D', read[5])
			for nro in match:
				readend += int(nro)
		if seq.startswith('AAA'):
			if dir == True:
				forwardaa.append({ 'start': readstart, 'name': read[0] })
		else:
			if 'AAA' in seq and 'S' not in read[5]:
				for m in re.finditer('(?=AAA)', seq):					
					if readstart + m.start() not in aaarepeats:
						aaarepeats.append(readstart + m.start())
					if readstart + m.start() +1 not in aaarepeats:
						aaarepeats.append(readstart + m.start()+1)
					if readstart + m.start() +2 not in aaarepeats:
						aaarepeats.append(readstart + m.start()+2)
					if readstart + m.start() +3 not in aaarepeats:
						aaarepeats.append(readstart + m.start()+3)
					if readstart + m.start() +4 not in aaarepeats:
						aaarepeats.append(readstart + m.start()+4)		
		if seq.endswith('TTT'):
			if dir == False:
				reversett.append({ 'end': readend, 'name': read[0] })
		else:
			if 'TTT' in seq and 'S' not in read[5]:
				for m in re.finditer('(?=TTT)', seq):
					if readstart + m.start() not in tttrepeats:
						tttrepeats.append(readstart + m.start())
					if readstart + m.start() +1 not in tttrepeats:
						tttrepeats.append(readstart + m.start()+1)
					if readstart + m.start() +2 not in tttrepeats:
						tttrepeats.append(readstart + m.start()+2)
					if readstart + m.start() +3 not in tttrepeats:
						tttrepeats.append(readstart + m.start()+3)
					if readstart + m.start() +4 not in tttrepeats:
						tttrepeats.append(readstart + m.start()+4)								
	
	maxindex = 0
	maxvalue = 0
	maxname = ""
	realaa = []
	realtt = []
	forwardstring = ''
	reversestring = ''
	
	if len(forwardaa) > 1:		
		readstring = ''
		for read in forwardaa:
			if read['start'] not in aaarepeats and abs(read['start'] - middle) < 100:
				realaa.append(read)
				readstring += str(read['start']) +':' +read['name'] +','
		if len(realaa) > 1:
			forwardstring = chr +":" +str(realaa[0]['start']+1) +"\tforward-AAA\t" +str(len(realaa)) +"\t" +readstring +"\n"
	
	if len(reversett) > 1:					
		readstring = ''
		for read in reversett:
			if read['end'] not in tttrepeats and abs(read['end'] - middle) < 100:
				realtt.append(read)
				readstring += str(read['end']) +':' +read['name'] +','
		if len(realtt) > 1:
			reversestring = chr +":" +str(realtt[0]['end']) +"\treverse-TTT\t" +str(len(realtt)) +"\t" +readstring +"\n"	
			
	if len(realaa) > 1 and len(realtt) > 1 and len(realtt) == len(realaa):
		discarded.write(line)
		continue
	
	if len(realaa) > len(realtt) and len(realaa) > 1:
		result.write(line.strip() +"\t" +forwardstring)
		if len(realtt) > 1:
			voted.write(line.strip() +"\t" +reversestring)
	if len(realtt) > len(realaa) and len(realtt) > 1:
		result.write(line.strip() +"\t" +reversestring)
		if len(realaa) > 1:
			voted.write(line.strip() +"\t" +forwardstring)

discarded.close()
voted.close()
result.close()