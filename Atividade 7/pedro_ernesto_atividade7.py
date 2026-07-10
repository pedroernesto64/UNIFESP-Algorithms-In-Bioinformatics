# Pedro Ernesto Duarte Pilchowski - RA: 156.331
# Atividade 7 - Algoritmos em Bioinformática

import numpy as np
import pandas as pd

def smith_waterman_gerar_csv(seq1, seq2, match=5, mismatch=-3, gap=-4):
    cols = len(seq1) + 1
    rows = len(seq2) + 1
    
    # Inicializa a matriz com zeros
    H = np.zeros((rows, cols), dtype=int)
    
    max_score = 0
    max_pos = []
    
    # Preenchimento da matriz de escore
    for i in range(1, rows):
        for j in range(1, cols):
            # Calcula os escores vindos da diagonal, de cima e da esquerda
            diag_score = H[i-1, j-1] + (match if seq2[i-1] == seq1[j-1] else mismatch)
            up_score = H[i-1, j] + gap
            left_score = H[i, j-1] + gap
            
            # Pega o maior valor (mínimo é 0)
            H[i, j] = max(0, diag_score, up_score, left_score)
            
            # Rastreia o(s) ponto(s) de pontuação máxima
            if H[i, j] > max_score:
                max_score = H[i, j]
                max_pos = [(i, j)]
            elif H[i, j] == max_score:
                max_pos.append((i, j))
                
    # SALVAR A MATRIZ EM CSV
    # Criando os cabeçalhos com as letras das sequências
    cols_labels = ['-'] + list(seq1)
    rows_labels = ['-'] + list(seq2)
    
    # Converte para DataFrame e salva
    df_matrix = pd.DataFrame(H, columns=cols_labels, index=rows_labels)
    nome_arquivo_matriz = 'matriz_escore.csv'
    df_matrix.to_csv(nome_arquivo_matriz)
    
    # TRACEBACK E SALVAR ALINHAMENTOS EM CSV
    alignments_data = []
    
    # Faz o traceback para todos os caminhos que atingiram o escore máximo
    for idx, pos in enumerate(max_pos):
        i, j = pos
        align1, align2, symbols = "", "", ""
        
        while H[i, j] != 0:
            current_score = H[i, j]
            diag_score = H[i-1, j-1]
            up_score = H[i-1, j]
            left_score = H[i, j-1]
            
            if current_score == diag_score + (match if seq2[i-1] == seq1[j-1] else mismatch):
                align1 += seq1[j-1]
                align2 += seq2[i-1]
                symbols += "|" if seq1[j-1] == seq2[i-1] else " "
                i -= 1
                j -= 1
            elif current_score == left_score + gap:
                align1 += seq1[j-1]
                align2 += "-"
                symbols += " "
                j -= 1
            elif current_score == up_score + gap:
                align1 += "-"
                align2 += seq2[i-1]
                symbols += " "
                i -= 1
                
        # Armazena os dados do alinhamento atual (invertidos pois o traceback é de trás pra frente)
        alignments_data.append({
            "Opcao_Alinhamento": f"Alinhamento {idx + 1}",
            "Escore_Total": max_score,
            "Sequencia_Topo (Seq1)": align1[::-1],
            "Matches": symbols[::-1],
            "Sequencia_Lateral (Seq2)": align2[::-1]
        })
        
    # Converte os alinhamentos para DataFrame e salva
    df_alignments = pd.DataFrame(alignments_data)
    nome_arquivo_alinhamentos = 'resultados_alinhamento.csv'
    df_alignments.to_csv(nome_arquivo_alinhamentos, index=False)

# EXECUÇÃO DO PROGRAMA
seq_lateral = "GGATCGA"
seq_topo = "GAATTCAGTTA"

smith_waterman_gerar_csv(seq_topo, seq_lateral)