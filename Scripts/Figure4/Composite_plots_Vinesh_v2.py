import sys, os, re
from operator import add
from optparse import OptionParser , IndentedHelpFormatter
from pylab import *
import numpy as np
#from matplotlib.font_manager import FontProperties
#import matplotlib.ticker as plticker
from scipy import stats

list1 = {}
## color schema Dark-RED-ROY-G-BIV Black, 
colors = ["#FF3333","#FF9933","#FFFF00","#4C9900","#0000FF","#4B0082","#9F00FF","#000000"]
#colors = ["#000000","#FF0000","#FF6600","#00FF00","#0000FF","#551a8b","#663300", "#660066", "#009999"]
#colors = ["#551a8b","#00FF00","#000000","#FF0000","#0000FF", "#660066", "#009999"]
#colors = ["#000000","#FF6600","#00FF00","#0000FF","#551a8b","#663300", "#660066", "#009999"]


    

def process_onestrand_files(infile,options,output_folder,ax,count):
    print "processing "+infile
    X = []
    label = os.path.basename(infile).split("_")[0]
    # Process the only CDT file 
    in_sense = open(infile,"rt")
    noL = 0
    for line in in_sense:
        if line.startswith("Uniqe") or line.startswith("ID") or line.startswith("gene"):
            tmp = line.rstrip().split("\t")[2:]
            X = [int(x) for x in tmp]
            xmin = min(X)
            xmax = max(X)
            Y = [0]*len(X)
            continue
        
        #tmplist = line.rstrip().split("\t")[602:1402]
        noL = noL + 1
        tmplist = line.rstrip().split("\t")[2:]
        newList = [float(x) for x in tmplist]
        Y = map(add,Y,newList)
        
    list1[label] = smoothListGaussian(Y,options.window)
    ##### REMOVE THE HASH WHEN YOU WANT TO DIVIDE BY NO OF GENES. ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
    #Y = [float(x)/noL for x in Y]
    plot_graph(X,Y,0,xmin,xmax,options,ax,label,count,noL)
    return(max(movingaverage(Y,options.window)))
    

def plot_graph(X,Y1,Y2,xmin,xmax,options,ax,label,count,noL):
    
    #X = smoothListGaussian(X,options.window)
    #Y1 = smoothListGaussian(Y1,options.window)
    #Y1 = [float(x)/noL for x in Y1]
    X = movingaverage(X,options.window)
    Y1 = movingaverage(Y1,options.window)
    # Allows you to plot y axis from 0-1
    if options.norm == 1:
        Y1 = [float(x)/max(Y1) for x in Y1]
        ax.plot(X, Y1, color=colors[count],label=label,lw=3.0)
        ## Change here for x axis limits and y axis limits. Nitika ########################################################################
        ax.set_xlim(-500,500)
        ax.set_ylim(0,1)
    
    else:
        
        ax.plot(X, Y1, color=colors[count],label=label,lw=3.0)
        ## Change here for x axis limits and y axis limits. Nitika ########################################################################
        ax.set_xlim(-500,500)
        ax.set_ylim(0,3500) ################################################################################################################################################
        
        
    

def movingaverage(interval, window_size):
    window= np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'valid')


def smoothListGaussian(list,degree):  

     window=degree*2-1  
     weight=np.array([1.0]*window)  
     weightGauss=[]  
     for i in range(window):  
         i=i-degree+1  
         frac=i/float(window)  
         gauss=1/(np.exp((4*(frac))**2))  
         weightGauss.append(gauss)  
     weight=np.array(weightGauss)*weight  
     smoothed=[0.0]*(len(list)-window)  
     for i in range(len(smoothed)):  
         smoothed[i]=sum(np.array(list[i:i+window])*weight)/sum(weight)
     #print smoothed
     new_smooth = [int(i) for i in smoothed]
     return new_smooth  


usage = '''
input_paths may be:
- a directory to run on all files in them

example usages:
python composite_for_many_factors_one_plot.py /usr/local/folder_containing_CDT_files
'''.lstrip()


 
# We must override the help formatter to force it to obey our newlines in our custom description
class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description


def run():
    parser = OptionParser(usage='%prog [options] input_paths', description=usage, formatter=CustomHelpFormatter())
    parser.add_option('-w', action='store', type='int', dest='window', default = 5,
                      help='Window size of moving average., Default=5')
    parser.add_option('-o', action='store', type='int', dest='norm',default=0,
                      help='1=> Divde by max normalization, Default => Plot absolute occupancy.')
    #parser.add_option('-p', action='store', type='string', dest='ppDir',
    #                  help='Directory containing peak-pairs.')
    #parser.add_option('-g', action='store', type='string', dest='genome_file',
    #                  help='File with chromosome lengths.')
    #parser.add_option('-i', action='store', type='int', dest='ilength', default = 40,
    #                  help='interval length., Default=40')

    (options, args) = parser.parse_args()
    
    if not args:
        parser.print_help()
        sys.exit(1)
        
    antisense_files = []
    sense_files = []
    
    output_folder = os.path.join(args[0],'_composite/') 
    if not os.path.exists(output_folder): os.makedirs(output_folder)
    outfile = os.path.join(output_folder,"composite_plot_all_factors.svg")
    
    
    if not os.path.exists(args[0]):
        parser.error('Path %s does not exist.' % args[0])
    if os.path.isdir(args[0]):
        for fname in os.listdir(args[0]):
            if fname.endswith(".txt") or fname.endswith(".cdt"):
                # intellignetly join paths without worrying about '/'
                fpath = os.path.join(args[0], fname)
                sense_files.append(fpath)
                    
        # No of subplots to plot
        nof = len(sense_files)
        
         # Declaring plotting parameters
        f,ax = subplots(1,1,sharex='all')
        #f.subplots_adjust(hspace=0)
        count = -1        
        
    ymax_list = []
    for f in sense_files:
        count = count + 1
        ymax = process_onestrand_files(f,options,output_folder,ax,count)
        ymax_list.append(ymax)
    
    if options.norm != 1:    
        ax.set_ylim(0,max(ymax_list))
    #lfp = FontProperties()
    #lfp.set_size(12)
    ax.legend()
    #ax.grid()
    savefig(outfile)
       
    
    
if __name__ == "__main__":
    run() 
