"""Merge best verified per-pattern programs from a solutions dir. Usage: merge.py input soldir output"""
import sys, os, glob
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'verify.py')).read().split("def main")[0])
pats = parse_input(open(sys.argv[1]).read())
best = {}
for f in glob.glob(os.path.join(sys.argv[2], '*.txt')):
    i, n = map(int, os.path.basename(f).split('_')[:2])
    if i in best and best[i][0] <= n: continue
    prog = [l.split() for l in open(f).read().split('\n') if l.strip()]
    if simulate(pats[i], prog) >= pats[i][3]: best[i] = (n, prog)
    else: print('invalid', f)
missing = [i for i in range(len(pats)) if i not in best]
print('missing', missing)
with open(sys.argv[3], 'w') as o:
    o.write('\n\n'.join('\n'.join(' '.join(l) for l in best.get(i, (0, [['-'] * 6]))[1]) for i in range(len(pats))) + '\n')
print('total', sum(best[i][0] for i in best))
