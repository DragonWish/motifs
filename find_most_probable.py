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

options = OptionParser(usage='%prog input output ', description="Specify input file"\
                       "containing the motif matrix. Specify output file for the boxplot to be "\
                       "saved as")


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

##print top 25 sequences in the sorted dictionary
##for i in sorted_scores[:25]:
##    prob_motif=most_probable(i[0],k,hpt_profile)[0]
##    print str(i[1])+"\t"+prob_motif+"\t"+str(titles.keys()[titles.values().index(i[0])])

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
    print str(len(r["fliers"][0].get_data()[1]))+" outliers"


def get_scores(proteome,prof,k):
    x1=time.time()
    titles=read_sequences(proteome)
    proteins=titles.values()
    maxp=0
    scores={}
    for i in proteins:
        kmer,prob=most_probable(i,k,prof)
        scores[i]=prob

    #sorted_scores=sorted(scores.iteritems(),key=operator.itemgetter(1), reverse=True)

    print str(time.time()-x1)+"  elapsed"
    return scores, titles

def main():
    opts, args = options.parse_args()
    
    #require 2 arguments
    if len(args) < 2:
        options.print_help()
        return
    
    motif_file=args[0]
    image_file=args[1]
    
    f=open(motif_file)
    prof=f.read()
    f.close()
    profile=ast.literal_eval(prof)

    k=len(profile['A'])
    
    organism_proteome="1773.fasta"

    scores, titles = get_scores(organism_proteome,profile,k)
    #t=read_sequences(organism_proteome)
    #first20=[probability(i[1:21],hpt_profile) for i in t.values()]
    #first20_e=[probability(i[1:21],hpt_profile_e) for i in t.values()]
    
##    sorted_scores=sorted(scores.iteritems(),key=operator.itemgetter(1), reverse=True)
##    for i in sorted_scores[:25]:
##        print str(i[1])+"\t"+str(most_probable(i[0],k,hpt_profile)[0])+"\t"+str(titles.keys()[titles.values().index(i[0])])
##    boxplot(scores)
    boxplot([scores.values()],["org"],image_file)


if __name__=="__main__":
    main()
    
