import commands, sys

ress = {}

Ns = (20, 22, 24, 26)

for N in Ns:
    print N
    res = dict(map(lambda x: x.split(':'), commands.getoutput("./selection %i" % 2**N).split('\n')))
    for k in res:
        ress[k] = ress.get(k, []) + [res[k].rjust(4)]

print 'Method', '&', ' & '.join('N = $2^%i$' % n for n in Ns), r'\\'

for k in sorted(ress):
    print k.rjust(16), '&', ' & '.join(ress[k]), r'\\'

