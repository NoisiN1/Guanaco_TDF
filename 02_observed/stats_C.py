#!/usr/bin/env python3
"""
Calcula los mismos estadisticos sobre .arp simulados y sobre el VCF empirico,
usando identicas formulas para que sean comparables.

Estadisticos:
  HT        heterocigosidad total (pool completo)
  FST_glob  Nei: (HT - HS)/HT sobre las 8 localidades
  FST_IM    Nei entre dos pools: isla vs continente
  H_isl     heterocigosidad esperada media, pool isla
  H_main    heterocigosidad esperada media, pool continente
  FIS_isl   1 - Ho/He en el pool isla

Uso:
  python3 stats_C.py arp  outscenario1_C/out1_1.arp
  python3 stats_C.py vcf  ../dataset_C_final.vcf
  python3 stats_C.py arpdir outscenario1_C 1 > sim1_multi_C.txt
"""
import sys, os, re, glob

MAINLAND = {"sg"}          # unica localidad continental


# ---------------------------------------------------------------- lectura ARP
def read_arp(path):
    """Devuelve {localidad: [(a1,a2), ...] por locus} -> dict de listas de
    genotipos; cada genotipo es una lista de pares (alelo1, alelo2) por locus."""
    pops, cur, hap = {}, None, []
    with open(path) as f:
        for line in f:
            s = line.rstrip("\n")
            m = re.search(r'SampleName="([^"]+)"', s)
            if m:
                cur = m.group(1); pops[cur] = []; hap = []
                continue
            if cur is None or "SampleData" in s or s.strip() in ("{", "}", ""):
                continue
            parts = s.split()
            if not parts:
                continue
            # linea con ID (ej. 1_1) o continuacion del segundo haplotipo
            if re.match(r'^\d+_\d+$', parts[0]):
                alleles = parts[1:]
            else:
                alleles = parts
            if all(a.strip().isdigit() for a in alleles):
                hap.append([int(a.strip()) for a in alleles])
                if len(hap) == 2:
                    pops[cur].append(list(zip(hap[0], hap[1])))
                    hap = []
    return {k: v for k, v in pops.items() if v}


# ---------------------------------------------------------------- lectura VCF
def read_vcf(path):
    samples, pops = [], {}
    with open(path) as f:
        for line in f:
            if line.startswith("##"):
                continue
            p = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                samples = p[9:]
                locs = [s.split("_")[0] for s in samples]
                for l in set(locs):
                    pops[l] = [[] for _ in range(locs.count(l))]
                idx = {}
                for i, l in enumerate(locs):
                    idx[i] = (l, sum(1 for j in range(i) if locs[j] == l))
                continue
            for i, g in enumerate(p[9:]):
                gt = g.split(":")[0].replace("|", "/")
                if gt in ("0/0", "0/1", "1/0", "1/1"):
                    a = tuple(int(x) for x in gt.split("/"))
                else:
                    a = None
                l, k = idx[i]
                pops[l][k].append(a)
    return pops


# ------------------------------------------------------------- estadisticos
def freqs(genos, locus):
    """frecuencia del alelo 1 y Ho en un locus, para una lista de individuos"""
    n = het = alt = 0
    for ind in genos:
        g = ind[locus]
        if g is None:
            continue
        if g[0] > 1 or g[1] > 1:
            return None, None
        n += 2; alt += g[0] + g[1]
        if g[0] != g[1]:
            het += 1
    if n == 0:
        return None, None
    return alt / n, het / (n / 2)


def he(p):
    return 2 * p * (1 - p)


def pool_stats(poplist, nloci):
    """He medio y Ho medio de un pool (lista de localidades)"""
    allg = [ind for pop in poplist for ind in pop]
    sHe = sHo = k = 0
    for L in range(nloci):
        p, ho = freqs(allg, L)
        if p is None:
            continue
        sHe += he(p); sHo += ho; k += 1
    return (sHe / k, sHo / k) if k else (float("nan"), float("nan"))


def nei_fst(poplist, nloci):
    """Nei: (HT - HS)/HT promediado sobre loci"""
    sHT = sHS = k = 0
    for L in range(nloci):
        ps, ws, allg = [], [], []
        for pop in poplist:
            p, _ = freqs(pop, L)
            if p is None:
                continue
            ps.append(p); ws.append(len(pop)); allg += pop
        if len(ps) < 2:
            continue
        pt, _ = freqs(allg, L)
        HS = sum(he(p) * w for p, w in zip(ps, ws)) / sum(ws)
        HT = he(pt)
        sHS += HS; sHT += HT; k += 1
    if not k or sHT == 0:
        return float("nan"), float("nan")
    return sHT / k, (sHT - sHS) / sHT


def compute(pops):
    nloci = len(next(iter(pops.values()))[0])
    isl  = [v for k, v in pops.items() if k not in MAINLAND]
    main = [v for k, v in pops.items() if k in MAINLAND]
    HT, FSTg = nei_fst(list(pops.values()), nloci)
    pool_i = [ind for pop in isl for ind in pop]
    pool_m = [ind for pop in main for ind in pop]
    _, FSTim = nei_fst([pool_i, pool_m], nloci) if pool_m else (0, float("nan"))
    Hi, Hoi = pool_stats(isl, nloci)
    Hm, _   = pool_stats(main, nloci) if main else (float("nan"), 0)
    FISi = 1 - Hoi / Hi if Hi else float("nan")
    return HT, FSTg, FSTim, Hi, Hm, FISi


# ------------------------------------------------------------------- main
def main():
    mode = sys.argv[1]
    if mode == "arp":
        r = compute(read_arp(sys.argv[2]))
    elif mode == "vcf":
        r = compute(read_vcf(sys.argv[2]))
    elif mode == "arpdir":
        d, tag = sys.argv[2], sys.argv[3]
        print("simid\tHT\tFST_glob\tFST_IM\tH_isl\tH_main\tFIS_isl")
        for f in sorted(glob.glob(os.path.join(d, f"out{tag}_*.arp")),
                        key=lambda x: int(re.search(r'_(\d+)\.arp$', x).group(1))):
            sid = re.search(r'_(\d+)\.arp$', f).group(1)
            try:
                v = compute(read_arp(f))
                print(sid + "\t" + "\t".join(f"{x:.6f}" for x in v))
            except Exception as e:
                print(f"{sid}\tNA\tNA\tNA\tNA\tNA\tNA", file=sys.stderr)
        return
    else:
        sys.exit("modo: arp | vcf | arpdir")
    print("HT\tFST_glob\tFST_IM\tH_isl\tH_main\tFIS_isl")
    print("\t".join(f"{x:.6f}" for x in r))


if __name__ == "__main__":
    main()
