import os, sys

args = sys.argv

minlen = False
keepfil = False
gccount = False

if 'C:\\' in args[-1]:
    fqfile = args[-1]
else:
    fqfile = os.getcwd() + '\\' + sys.argv[-1]

if '--min-length' in args:
    minlen = True
    min_read_len = args[args.index('--min-length') + 1]

if '--output_base_name' in args:
    outname = args[args.index('--output_base_name') + 1]
else:
    for i in range(-1, -len(args[-1]), -1):
        if args[-1][i] == '.':
            break
    outname = args[-1][:-6]

if '--keep_filtered' in args:
    keepfil = True
    outname_p = outname + '__passed.fastq'
    outname_f = outname + '__failed.fastq'
else:
    outname = outname

if '--gc_bounds' in args:
    gccount = [int(args[args.index('--gc_bounds') + 1])]
    if args[args.index('--gc_bounds') + 2].isdigit():
        gccount += [int(args[args.index('--gc_bounds') + 2])]

def write_outfq(name, read):
    with open(os.getcwd() + '\\' + name, "a") as fqout:
        fqout.write(read)

with open(fqfile, 'r') as fqsample:
    i = 0
    Read = []
    k = 0
    for line in fqsample:
        Read += [line]
        i += 1
        if i == 2:
            if minlen:
                #if
                pass
            pass
        if i == 4:

                print(Read)
                for i in range(4):
                    write_outfq(outname, Read[i])
            i = 0
            Read = []
        k += 1
        if k > 4:
            break

