#!/usr/bin/env python3

import codecs

table = {}
single = []

if __name__ == "__main__":
	with codecs.open("Big5-ZhuYin.map", "r", "Big5-hkscs", "replace") as filobj:
		for line in filobj:
			# print(line)
			wrd = line[0]
			single.append(wrd)
			for (i,seg) in enumerate(line.split()):
				if line[i-1] in [wrd, '/', ' ', None]:
					punc = seg[0]
					if punc in table:
						if wrd not in table[punc]:
							table[punc].append(wrd)
					else:
						table[punc] = [wrd]
			# print(st[0])
	of = open('ZhuYin-Big5.map', 'w', encoding="Big5-hkscs")

	for punc in table:
		of.write(punc+"\t")
		for wd in table[punc]:
			of.write(wd+" ")
		of.write("\r\n")
	for wd in single:
		of.write(wd+"\t"+wd+"\r\n")	
	of.close()

