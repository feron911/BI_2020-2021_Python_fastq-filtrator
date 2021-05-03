import os, sys

args = sys.argv
minlen = False
keepfil = False
gcfilt = False
names_of_files = {}

def file_name_parse(args):
    if 'C:\\' in args[-1]:
        fqfile = args[-1]
    else:
        fqfile = os.getcwd() + '\\' + args[-1]
    return fqfile

if '--help' in args:
    print("This script filter input fastQ file by minimal length \n "
          "(option --min_length [number]) and by GC content of read \n"
          "(option --gc_bounds [lower] [upper])")
    exit()

def min_length_parse(args):
    if '--min_length' in args:
        global minlen
        minlen = True
        if args[args.index('--min_length') + 1].isdigit():
            min_read_len = int(args[args.index('--min_length') + 1])
            return min_read_len
        else:
            print("Please check --min-length argument. Probably you forgot to add threshold.")
            exit()

def outname_parse(args):
    if '--output_base_name' in args:
        outname = args[args.index('--output_base_name') + 1]
    else:
        for i in range(-1, -len(args[-1]), -1):
            if args[-1][i] == '.':
                outname = args[-1][:i]
                break
        else:
            outname = args[-1]
    global names_of_files
    names_of_files["outname"] = outname
    return outname

def keep_filtered_parse(args):
    outname = outname_parse(args)
    if '--keep_filtered' in args:
        global keepfil
        global names_of_files
        keepfil = True
        outname_p = outname + '__passed.fastq'
        outname_f = outname + '__failed.fastq'
        names_of_files["outname_p"] = outname_p
        names_of_files["outname_f"] = outname_f
        return [outname_p, outname_f]
    else:
        return [outname, outname]

def gc_boundry_parse(args):
    if '--gc_bounds' in args:
        global gcfilt
        gcfilt = True
        if args[args.index('--gc_bounds') + 1].isdigit():
            gccount = [int(args[args.index('--gc_bounds') + 1])]
        else:
            print("Please check --gc_bounds argument. Probably you forgot to add threshold.")
            exit()
        if args[args.index('--gc_bounds') + 2].isdigit():
            gccount += [int(args[args.index('--gc_bounds') + 2])]
        return gccount

def write_outfq(name, read):
    with open(os.getcwd() + '\\' + name, "a") as fqout:
        fqout.write(read + "\n")

def filter_function(args):
    with open(file_name_parse(args), 'r') as fqsample:
        i = 0
        Read = []
        failed = False
        file_names = keep_filtered_parse(args)
        gccount = gc_boundry_parse(args)
        minlen_count = min_length_parse(args)
        file_names.append(names_of_files["outname"] + '.fastq')
        for name in file_names:
            if os.path.isfile(os.getcwd() + '\\' + name):
                raise ValueError ("Please, restart function with --output_base_name [base_name_for_new_file] option")
        for line in fqsample:
            Read += [line.strip()]
            i += 1
            if i == 2:
                if minlen:
                    if len(Read[1]) < minlen_count:
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
                        write_outfq(names_of_files["outname_f"], Read[j])
                elif keepfil:
                    for j in range(4):
                        write_outfq(names_of_files["outname_p"], Read[j])
                elif not failed:
                    for j in range(4):
                        write_outfq(names_of_files["outname"] + '.fastq', Read[j])
                i = 0
                Read = []
                failed = False
    return

if __name__ == "__main__":
    filter_function(args)
    file_names.append(names_of_files["outname"])