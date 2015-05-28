from optparse import OptionParser
import ast
import math
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.path as path
import time
import operator
import random

random.seed(2)

options = OptionParser(usage='%prog [input_motif] [list_proteomes] [output_file] ', description="Specify input file"\
                       "containing the motif matrix. Specify text file containing "\
                       "a list of proteomes to be examined. Specify output file for the outliers to be "\
                       "saved in")


def transpose(m):
    mT=[[]for x in m[0]]
    for i in range(0,len(m[0])):
        for j in range(0,len(m)):
            mT[i].append(m[j][i])
    return mT

def scramble_profile(prof,k):
    values=prof.values()
    aas=prof.keys()
    vT=transpose(values)
    random.shuffle(vT)
    v=transpose(vT)
    scrambled_prof={aas[i]:v[i] for i in range(0,k)}
    return scrambled_prof

#probability of a kmer occurring, given a profile p
def probability(kmer,p):
    prob=1
    for i in range(0,len(kmer)):
        try:
            prob=prob*p[kmer[i]][i]
        except KeyError:
            prob=0
    return prob

#returns most probable kmer in a given sequence, according
#to profile p
def most_probable(seq,k,p):
    k_most_probable=""
    maxp=-1
    for i in range(0,len(seq)-k):
        kmer=seq[i:i+k]
        probability=1
        for j in range(0,k):
            try:
                probability=probability*p[kmer[j]][j]
            except KeyError:pass
        if probability>maxp:
            maxp=probability
            k_most_probable=kmer
    return k_most_probable,maxp

#opens and reads a fasta file (filename), returns dictionary mapping each
#title (>...) to its corresponding protein sequence.
def read_sequences(filename):
    f=open(filename)
    lines=f.readlines()
    f.close()
    titles={}
    for i in lines:
        if i[0]==">":
            title=i[:-1]
            titles[title]=""
        else:
            titles[title]+=i.replace("X","A")[:-1]
    return titles


def histogram(scores,high,low): #histogram of scores
    fig,ax=plt.subplots()
    low=int(math.log(low[1],10))
    high=int(math.log(high[1],10))
    data = scores.values()
    pl.hist(data, bins=np.logspace(low, high, 50))
    pl.gca().set_xscale("log")
    pl.show()

def boxplot(scores, titles, outfile):    #boxplot of logged scores, whiskers at 3IQR
    logged_data=[[] for i in scores]
    for i in range(0,len(scores)):
        for j in scores[i]:
            if j>0:
                logged_data[i].append(math.log(j,10))
    r=plt.boxplot(logged_data, whis=3)

    pl.xticks(range(1,len(scores)+1),titles)
    #plt.show()
    pl.savefig(outfile)
    #print the number of outliers
    #print str(len(r["fliers"][0].get_data()[1]))+" outliers"


def get_scores(proteome,prof,k):
    x1=time.time()
    titles=read_sequences(proteome)
    proteins=titles.values()
    ptitles=titles.keys()
    maxp=0
    scores={}
    for i in ptitles:
        kmer,prob=most_probable(titles[i],k,prof)
        scores[i]=prob

    print proteome +" done"
    print str(time.time()-x1)+" s elapsed"
    return scores, titles

#get the index of scores where threshold c is reached    
def get_outliers_cutoff(scores,c):
    nn=0
    for i in scores:
        if math.log(i[1],10)>=c:nn+=1
        else:break
    return nn

#get the score value at the 3IQR upper whisker
def get_3iqr_cutoff(sortedscores):
    justscores=[math.log(i[1],10) for i in sortedscores if i[1]>0]
    scores_q1=len(justscores)/4
    scores_q3=scores_q1*3
    iqr=justscores[scores_q1]-justscores[scores_q3]
    scores_whisker=justscores[scores_q1]+3*iqr
    return scores_whisker
    
def main():
    opts, args = options.parse_args()
    
    #require 3 arguments
    if len(args) < 3:
        options.print_help()
        return
    
    motif_file=args[0]
    proteomes=args[1]
    outliers_file=args[2]
    
    f=open(motif_file)
    prof=f.readlines()
    f.close()
    use_iqr=True
    if prof[0][:3]!="IQR":
        use_iqr=False
        cutoff_value=float(prof[0])
        log_cutoff=math.log(cutoff_value,10)
    profile=ast.literal_eval(prof[1])

    k=len(profile['A'])
    
    #Add ambiguous amino acids:
    #B: Asx (N+D)
    profile['B']=[profile['N'][i]+profile['D'][i] for i in range(k)]
    #Z: Glx (Q+E)
    profile['Z']=[profile['Q'][i]+profile['E'][i] for i in range(k)]
    #J: Xle (I+L)
    profile['J']=[profile['I'][i]+profile['L'][i] for i in range(k)]
    #X: any (1)
    profile['X']=[1 for i in range(k)]
    
    f2=open(proteomes)
    organisms=f2.readlines()
    f2.close()

    org_names=[i[:-7] for i in organisms]
    org_scores=[]
    
    f3=open(outliers_file,'w')
    for i in organisms:
        f3.write(i)
        scores,titles=get_scores(i[:-1],profile,k)
        org_scores.append(scores.values())
        
        sorted_scores=sorted(scores.iteritems(),key=operator.itemgetter(1), reverse=True)
        
        if use_iqr==True:
            log_cutoff=get_3iqr_cutoff(sorted_scores)
            
        outlier_index=get_outliers_cutoff(sorted_scores,log_cutoff)
        for i in sorted_scores[:outlier_index]:
            f3.write(str(i[1])+"\t"
                     +str(most_probable(titles[i[0]],k,profile)[0])+"\t"
                     +str(i[0])+"\n")
        
    f3.close()

    #boxplot(org_scores,org_names,image_file)


if __name__=="__main__":
    main()
    
