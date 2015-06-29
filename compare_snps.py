from optparse import OptionParser

parser = OptionParser()
parser.add_option("--ref", type=str)
parser.add_option("--test", type=str)
parser.add_option("--out", type=str)
options,args = parser.parse_args()

ref = open(options.ref, 'r')
test = open(options.test, 'r')
out = open(options.out+'_acc.txt', 'w')

def sampleDict(vcf):
    samples = {}
    for line in vcf:
        if line.startswith('#'):
            if line.startswith('#CHROM'):
                header = line.strip().split()
                for i in range(9, len(header)):
                    samples[i-9] = header[i]
                return samples

def snpDict(line, sampledict, gprobs=False):
    ref = line[3]
    alt = line[4]
    snps = {}
    for i in range(9, len(line)):
        info = line[i].strip().split(':')
        snp = info[0]
        if gprobs:
            gps = info[2].split(',')
        if snp == '0/0' or snp == '0|0':
            gt = ref+'/'+ref
            if gprobs:
                gp = gps[0]
        elif snp == '1/0' or snp == '0/1' or snp == '1|0' or snp == '0|1':
            gt = ref+'/'+alt
            if gprobs:
                gp = gps[1]
        elif snp == '1/1' or snp == '1|1':
            gt = alt+'/'+alt
            if gprobs:
                gp = gps[2]
        elif snp == './.':
            gt = './.'
            gp = 0
        else:
            gt = './.'
            gp = 0
        if gprobs:
            snps[sampledict[i-9]] = str(gt)+':'+str(gp)
        else:
            snps[sampledict[i-9]] = gt
    return snps

def compareSnps(ref, test):
    match = 0
    missing = 0
    mismatch = 0
    na = 0
    rp = []
    wp = []
    for key in test:
        try:
            tgt = test[key].split(':')
            if tgt[0] == ref[key]:
                match += 1
#                rp.append(tgt[1])
            elif test[key][:3] == './.' or ref[key] == './.':
                missing += 1
            elif tgt[0][0] == ref[key][2] or tgt[0][2] == ref[key][0]:
                match += 1
            else:
                mismatch += 1
#                wp.append(tgt[1])
        except KeyError:
            na += 1
    return match,mismatch,missing,na,rp,wp

ref_samples = sampleDict(ref)
test_samples = sampleDict(test)

tmatch = 0
tmismatch = 0
tmissing = 0
keep_snps = 0

for line in test:
        if line.startswith('#'):
            continue
        else:
            tl = line.strip().split()
            found = False
            while not found:
                rline = ref.readline()
                if rline.startswith('#'):
                    continue
                else:
                    rl = rline.strip().split()
                    if not rl == []:
                        if tl[0] == rl[0] and tl[1] == rl[1]:
                            found = True
            if not found:
                continue
            elif found:
                match,mismatch,missing,na,rp,wp = compareSnps(snpDict(rl, ref_samples), snpDict(tl, test_samples))
                tmatch += match
                tmismatch += mismatch
                if (match+mismatch > 0):
                    new = tl[2]+'\t'+str(match/(match+mismatch))+'\t'+str(missing)+'\n'
                    keep_snps += 1
                    out.write(new)

print('Matching:', tmatch)
print('Mismatched:', tmismatch)
print('% Matching:', tmatch/(tmatch+tmismatch))
print('Missing:', tmissing)
print('Keep:', keep_snps)

