# Montar CSV
def cgCSV(dados):
    arquivo = []
    arquivo.append("Cabeçalho,Razão CG")
    for d in dados:
        # Arredonda e força a exibição com 3 casas decimais
        cg = f"{round((d[2]/100),3):.3f}"
        arquivo.append(f"{d[0]},{cg}")

    with open("Conteudo_CG.csv", "w") as file:
        for line in arquivo:
            file.write(f"{line}\n")


# Contar nucleotídeos
from Bio import SeqIO
from Bio.Seq import Seq

def contaNucleotideos(sequencia):
    a=sequencia.count("A")
    t=sequencia.count("T")
    c=sequencia.count("C")
    g=sequencia.count("G")
    n=sequencia.count("N")

    return a,t,c,g,n

# Ler arquivo Fasta
from Bio import SeqIO
from Bio.Seq import Seq
zaire = list(SeqIO.parse("Zaire Mpox cds_from_genomic.fna", "fasta"))

# Plotar gráfico
from matplotlib import pyplot as plt

def plotGrafico(resultados, nome_saida="grafico_gc_tm.png"):
    # Extração dos dados
    tm_eixo_x = [r[3] for r in resultados]
    gc_eixo_y = [r[2] for r in resultados]

    # Configuração do gráfico
    plt.figure(figsize=(10, 6))
    plt.scatter(tm_eixo_x, gc_eixo_y, alpha=0.6, color='blue', edgecolors='black')
    
    # Rótulos e Títulos
    plt.xlabel("Temperatura de Anelamento (Tm) [°C]")
    plt.ylabel("Conteúdo GC (%)")
    plt.title("Relação entre Conteúdo GC e Temperatura de Melting")
    plt.suptitle("Pedro Ernesto RA: 156.331 - Análise Mpox Zaire", fontsize=10, y=0.95)
    plt.grid(True, linestyle='--', alpha=0.7)

    # Salva o arquivo
    plt.savefig(nome_saida, dpi=300, bbox_inches='tight')