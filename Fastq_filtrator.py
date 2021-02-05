import os, sys, re

args = sys.argv

minlen = False
keepfil = False
gcfilt = False

if 'C:\\' in args[-1]:
    fqfile = args[-1]
else:
    fqfile = os.getcwd() + '\\' + sys.argv[-1]

if '--min-length' in args:
    minlen = True
    if args[args.index('--min-length') + 1].isdigit():
        min_read_len = int(args[args.index('--min-length') + 1])

if '--output_base_name' in args:
    outname = args[args.index('--output_base_name') + 1] + '.fastq'
else:
    for i in range(-1, -len(args[-1]), -1):
        if args[-1][i] == '.':
            outname = args[-1][:i]
            break
    else:
        outname = args[-1]

if '--keep_filtered' in args:
    keepfil = True
    outname_p = outname + '__passed.fastq'
    outname_f = outname + '__failed.fastq'


if '--gc_bounds' in args:
    gcfilt = True
    gccount = [int(args[args.index('--gc_bounds') + 1])]
    if args[args.index('--gc_bounds') + 2].isdigit():
        gccount += [int(args[args.index('--gc_bounds') + 2])]

def write_outfq(name, read):
    with open(os.getcwd() + '\\' + name, "a") as fqout:
        fqout.write(read)

with open(fqfile, 'r') as fqsample:
    i = 0
    Read = []
    failed = False
    for line in fqsample:
        Read += [line]
        i += 1
        if i == 2:
            if minlen:
                if len(Read[1]) < min_read_len:
                    if keepfil:
                        failed = True
            if gcfilt:
                gc_perc = (Read[1].count("G") + Read[1].count("C")) / len(Read[1]) * 100
                if len(gccount) == 1:
                    if gc_perc < gccount[0]:
                        failed = True
                if len(gccount) == 2:
                    if gc_perc < gccount[0] or gc_perc > gccount[1]:
                        failed = True
        if i == 4:
            if failed and keepfil:
                for j in range(4):
                    write_outfq(outname_f, Read[j])
            elif keepfil:
                for j in range(4):
                    write_outfq(outname_p, Read[j])
            elif not failed:
                for j in range(4):
                    write_outfq(outname, Read[j])
            i = 0
            Read = []
            failed = False
