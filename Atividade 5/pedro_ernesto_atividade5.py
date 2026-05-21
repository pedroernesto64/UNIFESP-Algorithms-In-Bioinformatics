# Pedro Ernesto Duarte Pilchowski - RA: 156.331
# Atividade 5 - Algoritmos em Bioinformática

import os
from Bio import SeqIO
from Bio.Seq import Seq

def encontrar_orfs(sequencia_nucleotideos, frame, seq_len, is_reverse=False):
    """
    Traduz a sequência no frame específico, encontra proteínas delimitadas por 
    um códon de iniciação (M) e parada (*), e calcula a localização original.
    """
    # Ajusta o deslocamento (offset) baseado no frame
    # Frame 1/4 -> offset 0 | Frame 2/5 -> offset 1 | Frame 3/6 -> offset 2
    offset = frame - 1 if not is_reverse else (frame - 4)
    
    # Garante que a sequência tenha tamanho múltiplo de 3 para não dar erro de tradução
    tamanho_util = len(sequencia_nucleotideos) - offset
    seq_analise = sequencia_nucleotideos[offset : offset + (tamanho_util - tamanho_util % 3)]
    
    # Traduz para proteína (código de 1 letra)
    proteina_seq = seq_analise.translate()
    
    orfs = []
    in_orf = False
    start_aa_idx = -1
    
    for i, aa in enumerate(proteina_seq):
        if not in_orf and aa == 'M':  # Encontrou códon iniciador (ATG)
            in_orf = True
            start_aa_idx = i
        elif in_orf and aa == '*':    # Encontrou códon de parada
            tamanho_proteina = i - start_aa_idx
            
            # Condição da atividade: considerar somente > 50 resíduos
            if tamanho_proteina > 50:
                proteina_encontrada = proteina_seq[start_aa_idx:i]
                
                # Calcula posições relativas na fita atual (0-based)
                start_nt = offset + start_aa_idx * 3
                end_nt = offset + i * 3 + 2
                
                # Converte para as posições reais da fita molde (forward) 1-based
                if is_reverse:
                    pos_inicio = seq_len - end_nt
                    pos_fim = seq_len - start_nt
                else:
                    pos_inicio = start_nt + 1
                    pos_fim = end_nt + 1
                    
                orfs.append({
                    'proteina': str(proteina_encontrada),
                    'inicio': pos_inicio,
                    'fim': pos_fim
                })
            
            # Reseta para buscar a próxima proteína
            in_orf = False
            
    return orfs

def main():
    ra = "156.331"
    arquivo_entrada = "Mpox Brazil partial genome sequence.fasta"

    if not os.path.exists(arquivo_entrada):
        print(f"Erro: O arquivo '{arquivo_entrada}' não foi encontrado.")
        return

    # Lê o genoma utilizando o Biopython
    record = SeqIO.read(arquivo_entrada, "fasta")
    seq_genoma = record.seq
    seq_genoma_rc = seq_genoma.reverse_complement()
    seq_len = len(seq_genoma)
    
    # Extrai o código do genoma do cabeçalho original 
    partes_id = record.id.split('|')
    codigo_genoma = partes_id[3] if len(partes_id) > 3 else "ON751962.1"

    print("Iniciando a busca por proteínas (ORFs > 50aa)...")
    
    # Loop pelos 6 frames (1, 2, 3 na fita direta | 4, 5, 6 na reversa)
    for frame in range(1, 7):
        if frame <= 3:
            orfs = encontrar_orfs(seq_genoma, frame, seq_len, is_reverse=False)
        else:
            orfs = encontrar_orfs(seq_genoma_rc, frame, seq_len, is_reverse=True)
        
        nome_arquivo_saida = f"gc_frame{frame}_ativ5_{ra}.fasta"
        
        # Grava os resultados
        with open(nome_arquivo_saida, "w") as f_out:
            for idx, orf in enumerate(orfs, start=1):
                # Formata o cabeçalho idêntico ao exigido
                cabecalho = f'> "{codigo_genoma}", Frame {frame}, proteína {idx}, [location={orf["inicio"]}..{orf["fim"]}], {ra}'
                
                # Quebra a sequência de aminoácidos a cada 80 caracteres
                seq_prot = orf['proteina']
                seq_formatada = "\n".join([seq_prot[j:j+80] for j in range(0, len(seq_prot), 80)])
                
                f_out.write(f"{cabecalho}\n{seq_formatada}\n")
                
        print(f"Frame {frame}: {len(orfs)} proteínas encontradas -> Salvo em {nome_arquivo_saida}")

    print("\nProcesso finalizado com sucesso!")

if __name__ == "__main__":
    main()