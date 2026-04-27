import math
from matplotlib import pyplot as plt
from Bio import SeqIO
from Bio.Seq import Seq
from processaDados import processaDados
from plotGrafico import plotGrafico
from montaCSV import *

zaire = list(SeqIO.parse("Zaire Mpox cds_from_genomic.fna", "fasta"))

# Acessar com dados[i][j]
# i = Entrada desejada
# j = (0: header, 1: contagem de nucleotídeos, 2: GC, 3: Tm)
dados = processaDados(zaire)

plotGrafico(dados, "Resultado_Mpox_Zaire.png")
nucleotideosCSV(dados)
cgCSV(dados)
temperaturaXgcCSV(dados)

# TO-DO
# * Imprima os dados em aquivos de saída no formato csv. Vem Envio item D.
# * Faça um gráfico (pontos/scatter) colocando na abscissa (eixo-x) a  temperatura de anelamento (melting)  e na ordenada o (eixo y) o conteúdo GC (em %) de cada sequência. Cada sequencia é um para ordenado (x,y).