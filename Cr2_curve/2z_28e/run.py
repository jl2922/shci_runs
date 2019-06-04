import os
import re

dirs = [d for d in os.listdir('.') if re.match('^r\d\.\d{2}$', d)]

print(dirs)

for d in dirs:
    if d in ['r2.25', 'r2.50', 'r3.00']:
        continue

    r = re.search('\d.*', d).group(0)
    print(r)
    os.system('rm -rf %s/*' % (d, ))
    fcidump = 'FCIDUMP_%s' % (r, )
    cmd = 'cp FCIDUMPS/%s %s && cd %s && ln -s %s FCIDUMP' % (fcidump, d, d, fcidump)
    print(cmd)
    os.system(cmd)
    cmd = 'cp config.json %s && cd %s && sbatch ../run_lm.sh' % (d, d)
    print(cmd)
    os.system(cmd)
