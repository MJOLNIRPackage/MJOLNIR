import sys

with open('docs/Tutorials/Keep.txt') as f:
	keepFiles = [x.strip() for x in f.readlines()]

if not len(sys.argv)==1:
	folder = sys.argv[1]

	import os
	content = os.listdir(folder)
	if folder[-1]!='/':
		folder+='/'
	if not len(content)==0:
		for f in content:
			if not f in keepFiles:
				os.remove(folder+f)
#			else:
#				print('Keeping: '+f)




