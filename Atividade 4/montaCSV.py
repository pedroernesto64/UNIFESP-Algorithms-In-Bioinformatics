def nucleotideosCSV(dados):
    arquivo = []
    arquivo.append("Cabeçalho,A,T,C,G")
    for d in dados:
        arquivo.append(f"{d[0]},{d[1][0]},{d[1][1]},{d[1][2]},{d[1][3]}")

    with open("Dados das sequencias.csv", "w") as file:
        for line in arquivo:
            file.write(f"{line}\n")

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

def temperaturaXgcCSV(dados):
    arquivo = []
    arquivo.append("Temperatura de Melting,Conteúdo CG(%)")
    for d in dados:
        melt=valor = f"{round((d[3]),3):.3f}"
        cg=valor = f"{round((d[2]),2):.2f}"
        arquivo.append(f"{melt},{cg}")

    with open("Temperatura_x_CG.csv", "w") as file:
        for line in arquivo:
            file.write(f"{line}\n")