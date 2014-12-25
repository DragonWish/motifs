import matplotlib.pyplot as plt
import math
from optparse import OptionParser
import random
import time

#Motif finding
#Program takes .txt or .fasta file of unaligned sequences as input, finds a representative motif
#of the peptide sequences provided. 
#Default is Gibbs

options = OptionParser(usage='%prog input outputfile k motif_file ', description="Specify a protein sequence file and an output filename \n"\
                "usage: FIND_MOTIF.py inputfile outputfile k motif_file \n"\
                "This will take the unaligned sequences in inputfile (could be .txt or .fasta) "\
                "and generate a motif profile that is saved as a python dictionary in the designated"\
                "motif_file. k is the specified length of the motif, and outputfile contains "\
                "the most likely motif selections from each sequence in inputfile. ")

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

#Transpose a matrix
def transpose(m): 
    for j in m:
        if len(j)!=10: print j, len(j)
    mT=[['' for j in m] for x in m[0]]
    for i in range(0,len(m[0])):
        for j in range(0,len(m)):
            mT[i][j]=m[j][i]
            
    return mT

#Generate a profile matrix (dimensions 20xk) given a collection of sequences 
#that are k-mers. The profile is generated with pseudocounts.
def profile(motifs): 
    m_=[[j for j in i] for i in motifs]
    m=transpose(m_)
    l=len(m)
    P={'A':[0 for i in range(0,l)],
       'L':[0 for i in range(0,l)],
       'V':[0 for i in range(0,l)],
       'I':[0 for i in range(0,l)],
       'F':[0 for i in range(0,l)],
       'M':[0 for i in range(0,l)],
       'S':[0 for i in range(0,l)],
       'Y':[0 for i in range(0,l)],
       'T':[0 for i in range(0,l)],
       'H':[0 for i in range(0,l)],
       'R':[0 for i in range(0,l)],
       'K':[0 for i in range(0,l)],
       'D':[0 for i in range(0,l)],
       'E':[0 for i in range(0,l)],
       'N':[0 for i in range(0,l)],
       'Q':[0 for i in range(0,l)],
       'C':[0 for i in range(0,l)],
       'G':[0 for i in range(0,l)],
       'W':[0 for i in range(0,l)],
       'P':[0 for i in range(0,l)]}
    aas=P.keys()
    h=[0 for i in range(0,l)]
    for i in range(0,l):
        h_=0
        for j in aas:
            x=(m[i].count(j)+1)/float(len(motifs))
            P[j][i]=x
            try:
                h_-=x*math.log(x,2)
            except ValueError:pass
        h[i]=h_
    return P

def score_motif(motifs,prof):
    x=0
    for i in range(0,len(motifs[0])):
        p=[prof[k][i]  for k in prof.keys()]
        col_score=len(motifs)*(1-max(p))
        x+=col_score
    return x

def score_m(motifs,t):
    x=0
    for i in range(0,len(motifs[0])):
        jj=[tots[i] for tots in motifs]
        col_score=t-max([jj.count('A'),jj.count('T'),
                         jj.count('G'),jj.count('C'),
                         jj.count('L'),jj.count('P'),
                         jj.count('V'),jj.count('W'),
                         jj.count('F'),jj.count('Q'),
                         jj.count('M'),jj.count('N'),
                         jj.count('Y'),jj.count('E'),
                         jj.count('S'),jj.count('D'),
                         jj.count('I'),jj.count('R'),
                         jj.count('H'),jj.count('K')])
        x+=col_score
    return x

def greedy_motif_search(dna,k,t):
    best_motifs=[i[:k] for i in dna]
    best_score=t*k
    for i in range(0,len(dna[0])-k):
        motif=[dna[0][i:i+k]]
        for j in range(1,t):
            p=profile(motif)
            motif.append(most_probable(dna[j],k,p))
        s=score_motif(motif,profile(motif))
        if s<best_score:
            best_score=s
            best_motifs=motif
    return best_motifs

def motif(p, k, dna):
    M=["" for i in dna]
    for i in range(0,len(dna)):
        M[i]=most_probable(dna[i],k,p)
    return M

def random_motif_select(dna,k,t):
    ms=["" for i in dna]
    for i in range(0,t):
        x=int(random.random()*(len(dna[i])-k))
        ms[i]=dna[i][x:x+k]
    return ms

def random_motif_search(dna,k,t):
    motifs=random_motif_select(dna,k,t)
    best_motifs=motifs
    best_score=50000
    for i in range(0,100):
        p=profile(motifs)
        motifs=motif(p,k,dna)
        s=score_m(motifs,t)
        if s<best_score:
            best_score=s
            best_motifs=motifs
    return best_motifs,best_score

def probability(kmer,p):            #calculates the probability of 
    prob=1                          #occurrence of a k-mer, given a
    for i in range(0,len(kmer)):    #motif profile p
        prob=prob*p[kmer[i]][i]
    return prob

def loaded_die(x):                      #biased probability selection
    xx=[float(sum(x[:i]))/sum(x) for i in range(1,len(x)+1)]
    p=random.random()
    i=len(x)-1
    for j in xx:
        if j>=p:
            i=xx.index(j)
            break
    return i

def gibbs_gen(pep,p,k):                     #randomly generate a kmer with a probability 
    probs=[0 for i in range(0,len(pep)-k)]  #bias according to profile
    for i in range(0,len(pep)-k):
        km=pep[i:i+k]
        probs[i]=probability(km,p)
    j=loaded_die(probs)
    kmer=pep[j:j+k]
    return kmer

def gibbs_motif_search(peptides,k,t,N):     #generates a set of k-mers from a sequence of
    motifs=random_motif_select(peptides,k,t)#peptides. 
    best_motifs=motifs
    best_score=5000
    for j in range(0,N):
        i=int(random.random()*t)
        motifs_=motifs[:i]+motifs[i+1:]
        p=profile(motifs_)
        motifs[i]=gibbs_gen(peptides[i],p,k)
        s=score_m(motifs,t)
        if s<best_score:
            best_score=s
            best_motifs=motifs
    return best_motifs,best_score

def most_probable(peptide,k,p):         #returns the most probable k-mer in a sequence
    nucs={"A":0,"C":1,"G":2,"T":3}      #given profile p
    k_most_probable=""
    maxp=-1
    for i in range(0,len(dna)-k):
        kmer=peptide[i:i+k]
        probability=1
        for j in range(0,k):
            probability=probability*p[kmer[j]][j]
        if probability>maxp:
            maxp=probability
            k_most_probable=kmer
    return k_most_probable

def main_greedy():        
    f=open("greedy_motif.txt")
    l=f.readlines()
    f.close()
    L=l[0].split()
    k=int(L[0])
    t=int(L[1])
    dna=l[1:]
    best=greedy_motif_search(dna,k,t)
    print best
def main_randomized():
    f=open("randomizes_motif.txt")
    l=f.readlines()
    f.close()
    L=l[0].split()
    k=int(L[0])
    t=int(L[1])
    dna=[i[:-1] for i in l[1:]]
    best_motif=[]
    best_score=40000
    for i in range(0,1000):
        motif,score=random_motif_search(dna,k,t)
        if score<best_score:
            best_score=score
            best_motif=motif
    return best_motif, best_score

def main_gibbs(seqfile,k):
    #for fasta formatted sequences:
    if seqfile[-5:]=="fasta":
        fasta_seqs=read_sequences(seqfile)
        proteins=fasta_seqs.values()
    #for text formatted files:
    elif seqfile[-3:]=="txt":
        f=open(seqfile)
        lines=f.readlines()
        f.close()
        proteins=[i[:-1] for i in lines]
    t=len(proteins)
    best_motif=[]
    best_score=40000
    for i in range(0,20):
        print "cycle no:"+str(i) +" "+ str(time.time()-x1) + " s"
        motif,score=gibbs_motif_search(proteins,k,t,2000)
        if score<best_score:
            best_score=score
            best_motif=motif
    return best_motif


x1=time.time()

def main():
    opts, args = options.parse_args()
    #require 4 arguments
    if len(args) < 4:
        options.print_help()
        return
    seqfile=args[0]
    outfile=args[1]
    k=int(args[2])
    motif_file=args[3]
    
    #write motif kmers in outfile
    best_gibbs=main_gibbs(seqfile,k)
    f=open(outfile,'w')
    for i in best_gibbs:f.write(i+"\n")
    f.close()
    
    #write motif profile matrix in motif_file
    f2=open(motif_file,'w')
    f2.write(str(profile(best_gibbs)))
    f2.close()
    
    print time.time()-x1, "s elapsed"

if __name__ == '__main__':
    main()
