from operator import itemgetter
import sys, os, math, operator
from optparse import OptionParser , IndentedHelpFormatter
from collections import OrderedDict




def  process_file(idxDir,options,outdir):
    
    back_tags = {}
    for fname in os.listdir(idxDir):
        if not (fname.endswith(".tab") or fname.endswith(".idx")):
            continue
        print "Calculating signal and background for "+fname
        idxData = {}
        total_reads = 0
        in0 = open(os.path.join(idxDir,fname),"rt")
        for line in in0:
            if line.startswith("#") or line.startswith("chrom"):
                continue
            cols = line.rstrip().split("\t")
            idxData[cols[0]+":"+cols[1]] = float(cols[2]) + float(cols[3])
            total_reads = total_reads + float(cols[2]) + float(cols[3])
        in0.close()
        
        # Find signal in GTFs
        
        total_signal = find_signal_in_GTF(options.ref,options.up,options.down,idxData)
        #total_background = total_reads - total_signal
        #back_tags[fname] = total_background
        back_tags[fname] = total_signal
    
    #max_back = max(back_tags.values())
    normalize_tags(back_tags,idxDir,outdir)

def normalize_tags(back_tags,idxDir,outdir):
    max_back = max(back_tags.values())
    for k,v in back_tags.items():
        print "Normalizing file = "+k
        in0 = open(os.path.join(idxDir,k),"rt")
        cf = round(max_back/v,2)
        outfile = os.path.join(outdir,os.path.splitext(k)[0]+"_"+str(cf)+".tab")
        out = open(outfile,"w")
        for line in in0:
            if line.startswith("#") or line.startswith("chrom"):
                continue
            cols = line.rstrip().split("\t")
            fwd = float(cols[2])*cf
            rev = float(cols[3])*cf
            tot = fwd+rev
            out.write(cols[0]+"\t"+cols[1]+"\t"+str(fwd)+"\t"+str(rev)+"\t"+str(tot)+"\n")
        out.close()
        in0.close()
            
    
    
    

def find_signal_in_GTF(ref,up,down,idxData):
    in1 = open(ref,"rt")
    total_sum = 0
    for line in in1:
        if line.startswith("#"):
            continue
        cols = line.rstrip().split("\t")
        if cols[6] == "+":
            start = int(cols[3]) - up
            end   = int(cols[3]) + down
        elif cols[6] == "-":
            start = int(cols[4]) - down
            end = int(cols[4]) + up
        total_sum = total_sum + return_sum(cols[0],start,end,idxData)
        
    return(total_sum)


              

def return_sum(chrom,start,end,idxData):
    summation = 0
    for j in range(start,end+1):
        key = chrom+":"+str(j)
        if key in idxData:
            summation = summation + idxData[key]
    
    return(summation)





usage = '''
input_paths may be:
- a directory to run on all files in them

example usages:
python Background_norm_for_GTFs_and_Rpb3.py [OPTIONS] /path/to/idx/directory/

# The total reads in the Input should be less than the total reads in
the treatment.
'''.lstrip()


 
# We must override the help formatter to force it to obey our newlines in our custom description
class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description


def run():
    parser = OptionParser(usage='%prog [options] input_paths', description=usage, formatter=CustomHelpFormatter())
    parser.add_option('-r', action='store', type='string', dest='ref',
                      help='The reference file in gff format.')
    parser.add_option('-u', action='store', type='int', dest='up',default=100,
                      help='Upstream distance from Ref., Default=100')
    parser.add_option('-d', action='store', type='int', dest='down',default=100,
                      help='Downstream distance from Ref. default=100')
    #parser.add_option('-f', action='store', type='int', dest='flag',default = 1,
    #                  help='1 => GTFs. 0 => Rpb3 only!. Default = 1')
    
    (options, args) = parser.parse_args()
    
    outdir = os.path.join(args[0],"norm_tab")
    if not os.path.exists(outdir): os.makedirs(outdir)
    
    #tmpdir = os.path.join(args[0],"tmpdir")
    #if not os.path.exists(tmpdir): os.makedirs(tmpdir)
    
    process_file(args[0],options,outdir)
    
    
if __name__ == "__main__":
    run() 