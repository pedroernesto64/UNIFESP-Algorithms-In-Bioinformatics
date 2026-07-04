# Pedro Ernesto Duarte Pilchowski - RA: 156.331
# Atividade 6 - Algoritmos em Bioinformática

import os
from Bio.Blast import NCBIWWW, NCBIXML

def ler_fasta(nome_arquivo):
    """
    Lê um arquivo FASTA e retorna uma lista com todas as sequências encontradas.
    Ignora as linhas de cabeçalho (>) e junta sequências que tiverem quebras de linha.
    """
    sequencias = []
    with open(nome_arquivo, 'r') as f:
        seq_atual = ""
        for linha in f:
            linha = linha.strip()
            if linha.startswith(">"):
                if seq_atual:
                    sequencias.append(seq_atual)
                    seq_atual = ""
            else:
                seq_atual += linha
        if seq_atual: # Adiciona a última sequência lida
            sequencias.append(seq_atual)
    return sequencias

def obter_melhor_sobreposicao(s1, s2):
    """
    Calcula a maior sobreposição onde o sufixo de s1 é igual ao prefixo de s2.
    Retorna o tamanho da sobreposição e a string fundida.
    """
    max_sobreposicao = 0
    string_fundida = s1 + s2
    
    tamanho_minimo = min(len(s1), len(s2))
    
    for i in range(1, tamanho_minimo + 1):
        sufixo = s1[-i:]
        prefixo = s2[:i]
        
        if sufixo == prefixo:
            max_sobreposicao = i
            string_fundida = s1 + s2[i:]
            
    return max_sobreposicao, string_fundida

def montar_contig_guloso(reads):
    """
    Aplica o algoritmo guloso para montar o contig iterativamente.
    Busca sempre o par com a maior sobreposição e os funde.
    """
    reads_atuais = reads.copy()
    
    while len(reads_atuais) > 1:
        max_sobreposicao_global = -1
        melhor_par = (0, 0)
        melhor_fusao = ""
        
        for i in range(len(reads_atuais)):
            for j in range(len(reads_atuais)):
                if i != j:
                    sobreposicao, fusao = obter_melhor_sobreposicao(reads_atuais[i], reads_atuais[j])
                    
                    if sobreposicao > max_sobreposicao_global:
                        max_sobreposicao_global = sobreposicao
                        melhor_par = (i, j)
                        melhor_fusao = fusao
                        
        idx1, idx2 = melhor_par
        if idx1 > idx2:
            reads_atuais.pop(idx1)
            reads_atuais.pop(idx2)
        else:
            reads_atuais.pop(idx2)
            reads_atuais.pop(idx1)
            
        reads_atuais.append(melhor_fusao)
        
    return reads_atuais[0]

def salvar_fasta(sequencia, nome_arquivo, cabecalho):
    """ Salva o contig em um arquivo formato FASTA. """
    with open(nome_arquivo, "w") as f:
        f.write(f"{cabecalho}\n")
        for i in range(0, len(sequencia), 70):
            f.write(sequencia[i:i+70] + "\n")

def fazer_blast_api(sequencia_fasta):
    """
    Submete a sequência ao BLAST via API utilizando Biopython.
    """
    print("\nSubmetendo sequência ao BLAST via API (Pode demorar)")
    
    result_handle = NCBIWWW.qblast("blastn", "nt", sequencia_fasta)
    
    print("Lendo resultados do XML")
    blast_record = NCBIXML.read(result_handle)
    
    print("\n--- TOP 5 ALINHAMENTOS ENCONTRADOS ---")
    
    if len(blast_record.alignments) == 0:
        print("Nenhum alinhamento encontrado para esta sequência.")
        return
    
    for i, alignment in enumerate(blast_record.alignments[:5]):
        for hsp in alignment.hsps:
            print(f"\nAlinhamento {i+1}:")
            print(f"Organismo/Sequência: {alignment.title}")
            print(f"Tamanho: {alignment.length}")
            print(f"E-value: {hsp.expect}")
            print(f"Score (Bits): {hsp.bits}")
            print(f"Identidade: {hsp.identities}/{hsp.align_length} ({(hsp.identities/hsp.align_length)*100:.2f}%)")

if __name__ == "__main__":
    # Lê os reads a partir do arquivo reads4.fasta
    reads = ler_fasta("reads4.fasta")
    
    # Monta o contig
    contig_final = montar_contig_guloso(reads)
    print("\nContig Montado:", contig_final)
    
    # Salva o resultado
    cabecalho = f">Contig montado utilizando Algoritmo Guloso RA 156.331"
    salvar_fasta(contig_final, "contig.fasta", cabecalho)
    
    # Lê a sequência do gene Mpox a partir do arquivo gene_Mpox.fasta para o Blast
    gene_mpox = ler_fasta("gene_Mpox.fasta")[0]
    fazer_blast_api(gene_mpox)