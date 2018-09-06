#This small script is for turning a username into a 5-digit seed.
#should work for both python 2.7 and 3.6
import hashlib
import sys
tohash = " ".join(sys.argv[1:])
if len(tohash) < 1:
    print("usage: python3 myseed.py kjetilmg")
    sys.exit(1)
tohash = tohash.encode('latin-1').lower().strip()
hexstr = hashlib.sha224(tohash).hexdigest()
dec = str(int(hexstr, 16))
print(int(dec[-5:]))
