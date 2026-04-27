from contaNucleotideos import contaNucleotideos
import math
from Bio import SeqIO
from Bio.Seq import Seq

def processaDados(dados):
    resultados=[]
    
    for entrada in dados:
        header = entrada.description
        n = contaNucleotideos(entrada.seq)

        # Calcula o Conteúdo GC
        gc = (n[2]+n[3])*100/math.fsum(n)

        # Calcula temperatura de anelamento
        tm = 64.9 + 0.41*gc - (500/len(entrada.seq))

        # Salva os cabeçalhos e as contagens como uma tupla em um vetor
        tupla = (header, n, gc, tm)
        resultados.append(tupla)
        
    return resultados