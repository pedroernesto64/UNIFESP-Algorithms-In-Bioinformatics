# Pedro Ernesto Duarte Pilchowski - RA: 156.331
# Atividade 3 - Algoritmos em Bioinformática

from Bio import SeqIO
from Bio.Seq import Seq

# Dicionário chave-valor, com a chave sendo a trinca
# genética e o valor sendo a proteína sintetizada
dicionarioProteina = {
    "UUU" : "F", "CUU" : "L", "AUU" : "I", "GUU" : "V",
    "UUC" : "F", "CUC" : "L", "AUC" : "I", "GUC" : "V",
    "UUA" : "L", "CUA" : "L", "AUA" : "I", "GUA" : "V",
    "UUG" : "L", "CUG" : "L", "AUG" : "M", "GUG" : "V",
    "UCU" : "S", "CCU" : "P", "ACU" : "T", "GCU" : "A",
    "UCC" : "S", "CCC" : "P", "ACC" : "T", "GCC" : "A",
    "UCA" : "S", "CCA" : "P", "ACA" : "T", "GCA" : "A",
    "UCG" : "S", "CCG" : "P", "ACG" : "T", "GCG" : "A",
    "UAU" : "Y", "CAU" : "H", "AAU" : "N", "GAU" : "D",
    "UAC" : "Y", "CAC" : "H", "AAC" : "N", "GAC" : "D",
    "UAA" : "*", "CAA" : "Q", "AAA" : "K", "GAA" : "E",
    "UAG" : "*", "CAG" : "Q", "AAG" : "K", "GAG" : "E",
    "UGU" : "C", "CGU" : "R", "AGU" : "S", "GGU" : "G",
    "UGC" : "C", "CGC" : "R", "AGC" : "S", "GGC" : "G",
    "UGA" : "*", "CGA" : "R", "AGA" : "R", "GGA" : "G",
    "UGG" : "W", "CGG" : "R", "AGG" : "R", "GGG" : "G" 
}

def geraFastaRNA(dna):
    with open("RNA_Zaire_mpox.fasta", "w") as f:
        # Itera sobre cada registro (cabeçalho + sequência)
        for registro in dna:
            # Adiciona nome e RA em cada cabeçalho
            cabecalho = f">Nome: Pedro Ernesto Duarte Pilchowski, RA: 156.331, {registro.description}"
            # Transcreve para RNA usando o método .transcribe()
            rna = str(registro.seq.transcribe())
            # Escreve o cabeçalho no arquivo, e em seguida a
            # sequência de RNA, limitando a 70 bases por linha
            f.write(f"{cabecalho}\n")
            for i in range(0, len(rna), 70):
                f.write(f"{rna[i:i + 70]}\n")

def geraFastaProteinaDicionario(dna):
    with open("Proteina_Zaire_mpox_Dicionario.fasta", "w") as f:
        for registro in dna:
            # Reinicia a sequência de proteínas a cada novo registro
            seqProteina = ""
            cabecalho = f">Nome: Pedro Ernesto Duarte Pilchowski, RA: 156.331, {registro.description}"
            # Obtém o RNA do registro em forma de string
            rna = str(registro.seq.transcribe())
            # Separa o RNA em trincas
            for i in range(0, len(rna), 3):
                trinca = rna[i:i+3]
                # Confere a proteína relativa à trinca pelo dicionário,
                # e guarda o resultado na string seqProteina
                seqProteina = seqProteina + dicionarioProteina[trinca]
            # Escreve o cabeçalho e as proteínas no arquivo,
            # limitando a 70 caracteres por linha
            f.write(f"{cabecalho}\n")
            for i in range(0, len(seqProteina), 70):
                f.write(f"{seqProteina[i:i + 70]}\n")

def geraFastaProteinaBiopython(dna):
    with open("Proteina_Zaire_mpox_Biopython.fasta", "w") as f:
        for registro in dna:
            cabecalho = f">Nome: Pedro Ernesto Duarte Pilchowski, RA: 156.331, {registro.description}"
            # Faz a mesma coisa que a função anterior, porém as proteínas
            # são traduzidas usando o método .translate(), que guarda o
            # resultado em forma de string na variável "proteina"
            proteina = str(registro.seq.translate())
            f.write(f"{cabecalho}\n")
            for i in range(0, len(proteina), 70):
                f.write(f"{proteina[i:i + 70]}\n")

# Abre o arquivo e chama as respectivas funções
zaire = list(SeqIO.parse("Zaire Mpox cds_from_genomic.fna", "fasta"))
geraFastaRNA(zaire)
geraFastaProteinaDicionario(zaire)
geraFastaProteinaBiopython(zaire)