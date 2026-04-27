from Bio import SeqIO
from Bio.Seq import Seq

def contaNucleotideos(sequencia):
    a=sequencia.count("A")
    t=sequencia.count("T")
    c=sequencia.count("C")
    g=sequencia.count("G")
    n=sequencia.count("N")

    return a,t,c,g,n