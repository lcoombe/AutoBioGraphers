import sys
import re

sequences = []

if len(sys.argv) != 3:
	print "Usage: %s <interaction txt> <protein sequences to filter>" % (sys.argv[0])
	sys.exit(1)

interaction_fileName = sys.argv[1]
sequences_fileName = sys.argv[2]

interaction = open(interaction_fileName, 'r')
for line in interaction:
	line = line.strip()
	split_line = re.split(r'\s+', line)
	if split_line[0] not in sequences:
		sequences.append(split_line[0])
	if split_line[1] not in sequences:
		sequences.append(split_line[1])

interaction.close()

seq_header = ""
filter_seq = False
seq = ""

sequences_file = open(sequences_fileName, 'r')
for line in sequences_file:
	line = line.strip()
	header_match = re.match(r'>(\S+)', line)
	if header_match != None:
		if filter_seq is True:
			print ">" + seq_header
			print seq
			seq = ""
			filter_seq = False
		seq_header = header_match.group(1)
		if seq_header in sequences:
			filter_seq = True
	else:
		if filter_seq is True:
			seq += line 

if filter_seq is True and seq != "":
	print ">" + seq_header
	print seq

sequences_file.close()