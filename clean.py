import sys

if not len(sys.argv)==1:
	folder = sys.argv[1]

	import os
	content = os.listdir(folder)
	if folder[-1]!='/':
		folder+='/'
	if not len(content)==0:
		for f in content:
			os.remove(folder+f)




