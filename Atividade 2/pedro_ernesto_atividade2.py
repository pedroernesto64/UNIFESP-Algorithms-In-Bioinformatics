# Pedro Ernesto Duarte Pilchowski - RA: 156.331
# Atividade 2 - Algoritmos em Bioinformática

from Bio import SeqIO
from Bio.Seq import Seq


# Escaneia o arquivo e separa seu conteúdo (descrição e genoma) em uma classe
# e salva na variável "brazil"
brazil = next(SeqIO.parse("Mpox Brazil partial genome sequence.fasta", "fasta"))
# Guarda apenas a descrição na variável "headerBrazil"
headerBrazil = (brazil.description)
# Adiciona aspas ao descrição, para que ele seja tratado inteiramente como uma
# string ao ser adicionado ao CSV (já que a presença de vírgulas seria interpretada
# como uma nova coluna sem as aspas)
headerBrazil = '"' + headerBrazil + '"'
# Guarda o genoma em "seqBrazil"
seqBrazil = (brazil.seq)

zaire = next(SeqIO.parse("Mpox Zaire-9 sequence.fasta", "fasta"))
headerZaire = (zaire.description)
headerZaire = '"' + headerZaire + '"'
seqZaire = (zaire.seq)

# Inicia vetores vazios. Cada item do vetor será uma linha no arquivo CSV escrito
arquivo1=[]
arquivo2=[]

# Função de análise do genoma usando BioPython. Recebe o genoma inteiro, e incrementa as
# variáveis correspondentes a cada nucleotídeo à medida que vai encontrando-as na sequência
def contaComBiopython(sequencia):
    a=sequencia.count("A")
    t=sequencia.count("T")
    c=sequencia.count("C")
    g=sequencia.count("G")
    n=sequencia.count("N")

    # Retorna uma tupla com os resultados
    return a,t,c,g,n

# Função de contagem elaborada na Atividade 1
def contaSemBiopython(sequencia):
    a=t=c=g=n=0
    for nucleotideo in sequencia:
        match nucleotideo:
            case 'A':
               a+=1
            case 'T':
                t+=1
            case 'C':
                c+=1
            case 'G':
                g+=1
            case 'N':
                n+=1
    return a,t,c,g,n

# Adiciona o cabeçalho ao vetor. Cada append é um novo item, isto é, uma nova linha.
# Cada vírgula representa uma coluna no arquivo CSV.
arquivo1.append("Resultados utilizando código da Atividade 1,,,,,,")
arquivo1.append("Nome: Pedro Ernesto Duarte Pilchowski,RA: 156.331,,,,,")
arquivo1.append("Vírus,A,T,C,G,N,Total")

# Chama a função de contagem e guarda a tupla retornada na variável "b"
b=contaSemBiopython(seqBrazil)
# Soma os valores da variável "b", acessando os índices da tupla
sumB=b[0]+b[1]+b[2]+b[3]+b[4]
# "lineBrazil" contém a descrição seguida dos valores dos nucleotídeos, separados por vírgulas
lineBrazil=headerBrazil+","+(f"{b[0]},{b[1]},{b[2]},{b[3]},{b[4]},{sumB}")
arquivo1.append(lineBrazil)

z=contaSemBiopython(seqZaire)
sumZ=z[0]+z[1]+z[2]+z[3]+z[4]
lineZaire=headerZaire+","+(f"{z[0]},{z[1]},{z[2]},{z[3]},{z[4]},{sumZ}")
arquivo1.append(lineZaire)

# Abre (ou cria, caso não exista) o arquivo CSV em modo de escrita
with open("Ativ2_Resultado_codigo_Ativ1.csv", "w") as file:
    # Percorre o vetor "arquivo1", e escreve linha-a-linha no arquivo
    for line in arquivo1:
        file.write(f"{line}\n")

########################################################################


# Repete o mesmo processo para os resultados utilizando módulo Seq
arquivo2.append("Resultados utilizando módulo Seq,,,,,,")
arquivo2.append("Nome: Pedro Ernesto Duarte Pilchowski,RA: 156.331,,,,,")
arquivo2.append("Vírus,A,T,C,G,N,Total")

b=contaComBiopython(seqBrazil)
sumB=b[0]+b[1]+b[2]+b[3]+b[4]
lineBrazil=headerBrazil+","+(f"{b[0]},{b[1]},{b[2]},{b[3]},{b[4]},{sumB}")
arquivo2.append(lineBrazil)

z=contaComBiopython(seqZaire)
sumZ=z[0]+z[1]+z[2]+z[3]+z[4]
lineZaire=headerZaire+","+(f"{z[0]},{z[1]},{z[2]},{z[3]},{z[4]},{sumZ}")
arquivo2.append(lineZaire)

with open("Ativ2_Resultado_Modulo_Seq.csv", "w") as file:
    for line in arquivo2:
        file.write(f"{line}\n")