import glob, re, sys
from stats_C import read_arp, freqs
d, tag = sys.argv[1], sys.argv[2]
print("simid\tsfs_low\tsfs_high")
fs = sorted(glob.glob(d+"/out%s_*.arp"%tag),
            key=lambda x:int(re.search(r'_(\d+)\.arp$',x).group(1)))
for f in fs:
    sid = re.search(r'_(\d+)\.arp$', f).group(1)
    try:
        pops = read_arp(f)
        nl = len(next(iter(pops.values()))[0])
        allg = [i for pop in pops.values() for i in pop]
        lo=hi=n=0
        for L in range(nl):
            p,_ = freqs(allg,L)
            if p is None: continue
            m=min(p,1-p); n+=1
            if m<0.1: lo+=1
            elif m>=0.4: hi+=1
        print("%s\t%.5f\t%.5f" % (sid, lo/n, hi/n))
    except Exception:
        print("%s\tNA\tNA" % sid)
