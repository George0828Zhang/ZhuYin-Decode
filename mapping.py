#!/usr/bin/env python3
import sys
import codecs

table = {}
single = []

if __name__ == "__main__":
	with codecs.open(sys.argv[1], "r", "Big5-hkscs", "replace") as filobj:
		for line in filobj:
			# print(line)
			wrd = line[0]
			single.append(wrd)
			# wrd seg/seg/seg
			# seg = punc-?-?
			for seg in line[1:].strip().split('/'):
				# print(seg)
				punc = seg[0]
				if punc in table:
					if wrd not in table[punc]:
						table[punc].append(wrd)
				else:
					table[punc] = [wrd]
			# print(st[0])
	of = open(sys.argv[2], 'w', encoding="Big5-hkscs")

	for punc in table:
		of.write(punc+"\t")
		for wd in table[punc]:
			of.write(wd+" ")
		of.write("\r\n")
	for wd in single:
		of.write(wd+"\t"+wd+"\r\n")	
	of.close()

