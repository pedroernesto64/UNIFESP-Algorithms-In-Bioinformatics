from matplotlib import pyplot as plt

def plotGrafico(resultados, nome_saida):
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