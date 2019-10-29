from optparse import OptionParser , IndentedHelpFormatter
import sys, os, re, operator, math


def process_file(idxDir,options,outdir):
    
    input0 = open(options.chrom_len,"rt")
    max_size = {}
    for line in input0:
        cols = line.rstrip().split("\t")
        max_size[cols[0]] = int(cols[1])
    input0.close()
    
    
    ## Reading the shift distance file
    #input1 = open(options.config,"rt")
    #shift_dist = {}
    #for line in input1:
    #    if line.startswith("fname"):
    #        continue
    #    cols = line.rstrip().split("\t")
    #    if options.flag == 0:
    #        # Get every thing in the list except last item
    #        tab_fname = "_".join(cols[0].split("_")[:-1])+".tab"
    #        if float(cols[1]) < 0:
    #            shift_dist[tab_fname] = 0
    #        else:
    #            shift_dist[tab_fname] = abs(int(math.ceil((float(cols[1])/2))))
    #        
    #    else:
    #        shift_dist[cols[0]] = int(cols[1])
    #input1.close()
    
    for fname in os.listdir(idxDir):
        fwd = {}
        rev = {}
        index = {}
        if not (fname.endswith(".tab") or fname.endswith("idx")):
            continue
        #dist = shift_dist[fname]
        dist = 6
        outfile = os.path.join(outdir,os.path.splitext(fname)[0]+"_shift"+str(dist)+".tab")
        out = open(outfile,"w")
        tmp = os.path.join(outdir,"tmp.tab")
        in1 = open(os.path.join(idxDir,fname),"rt")
        for line in in1:
            if line.startswith("#") or line.startswith("chrom"):
                continue
            cols = line.rstrip().split("\t")
            chrom = cols[0]
            if float(cols[2]) > 0:
                new_fwd_pos = int(cols[1]) + dist
                if new_fwd_pos <= max_size[chrom]:
                    fwd[chrom+"\t"+str(new_fwd_pos)] = float(cols[2])
                    index[chrom+"\t"+str(new_fwd_pos)] = 1
                else:
                    fwd[chrom+"\t"+cols[1]] = float(cols[2])
                    index[chrom+"\t"+cols[1]] = 1
                
            if float(cols[3]) > 0:
                new_rev_pos = int(cols[1]) - dist  
                if new_rev_pos > 0:
                    rev[chrom+"\t"+str(new_rev_pos)] = float(cols[3])
                    index[chrom+"\t"+str(new_rev_pos)] = 1
                else:
                    rev[chrom+"\t"+cols[1]] = float(cols[3])
                    index[chrom+"\t"+cols[1]] = 1
        in1.close()
        
        for k,v in index.items():
            if k in fwd and k in rev:
                out.write(k+"\t"+str(fwd[k])+"\t"+str(rev[k])+"\n")
            elif k in fwd and k not in rev:
                out.write(k+"\t"+str(fwd[k])+"\t0.0\n")
            elif k not in fwd and k in rev:
                out.write(k+"\t0.0\t"+str(rev[k])+"\n")
        out.close()
        os.system("sort -k1,1 -k2,2n "+outfile+" >"+tmp)
        os.system("mv  "+tmp+" "+outfile)












##Make sure if your input is shift.config file then it only contains two columns. Column 1 = > Filename [ WITHOUT PATHS APPENDED]
##column 2 = > shift distance.


usage = '''
input_paths may be:
- a directory to run on all files in them

example usages:
python shift_tags.py  [OPTIONS] <path_to_tab_folder>

BY DEFAULT SHIFT DISTANCE WILL BE 6.
'''.lstrip()


 
# We must override the help formatter to force it to obey our newlines in our custom description
class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description


def run():
    parser = OptionParser(usage='%prog [options] input_paths', description=usage, formatter=CustomHelpFormatter())
    #parser.add_option('-i', action='store', type='string', dest='config',
    #                  help='Input file containing shift distance.')
    #parser.add_option('-f', action='store', type='int', dest='flag',default=0,
    #                  help='0 => the shift file == "statistics.txt" from cw-peak-pair output, 1 => "shift.config" file')
    parser.add_option('-g', action='store', type='string', dest='chrom_len',
                      help='File containing chromosome length.')
    
    (options, args) = parser.parse_args()
    
    if not args:
        parser.print_help()
        sys.exit(1)
        
    outdir = os.path.join(os.path.dirname(args[0]),"shifted_tab")
    if not os.path.exists(outdir): os.makedirs(outdir)
    
     
    if not os.path.exists(args[0]):
        parser.error('Path %s does not exist.' % args[0])
    
    process_file(args[0],options,outdir)
    
    
    
if __name__ == "__main__":
    run() 