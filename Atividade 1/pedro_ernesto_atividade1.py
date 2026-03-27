# Pedro Ernesto Duarte Pilchowski - RA: 156.331
# Atividade 1 - Algoritmos em Bioinformática

arquivo = "Monkeypox_virus-trechos.fasta"

# Recebe uma linha do genoma e conta quantos nucleotídeos de cada tipo há
def conta_nucleotideos(sequencia):
    a=t=c=g=u=0
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
            case 'U':
                u+=1            
    return a,t,c,g,u

sum_a=sum_t=sum_c=sum_g=sum_u=0
# Prepara o arquivo resultado.txt para escrita
with open("resultado.txt", 'w') as w:
    with open(arquivo, 'r') as f:
        # Abre o arquivo e o percorre linha-a-linha
        for l in f:
            # Verifica se é um cabeçalho
            if l[0] == '>':
                soma=sum_a+sum_t+sum_c+sum_g+sum_u
                # Verifica se há um genoma já processado que ainda não foi impresso
                if (soma!=0):
                    print(f"A: {sum_a}\nT: {sum_t}\nC: {sum_c}\nG: {sum_g}\nU: {sum_u}\nTotal: {soma} nucleotídeos\n", file=w)
                    sum_a=sum_t=sum_c=sum_g=sum_u=0
                # Imprime o cabeçalho
                print(l, end="", file=w)

            # Acumula a contagem dos nucleotídeos
            else:
                a,t,c,g,u = conta_nucleotideos(l)
                sum_a+=a
                sum_t+=t
                sum_c+=c
                sum_g+=g
                sum_u+=u
        
        soma=sum_a+sum_t+sum_c+sum_g+sum_u
        # Imprime o último genoma do arquivo
        print(f"A: {sum_a}\nT: {sum_t}\nC: {sum_c}\nG: {sum_g}\nU: {sum_u}\nTotal: {soma} nucleotídeos", end="", file=w)