import numpy as np
import matplotlib.pyplot as plt

# Configura uma semente fixa para reprodutibilidade
np.random.seed(42)

N = 500
MAX_ITER = 20

def run_cg_experiment(tau):
    print("\n" + "="*50)
    print(f" Executando Gradientes Conjugados para tau = {tau}")
    print("="*50)

    # 1. Construção da Matriz A
    A = np.random.uniform(-1.0, 1.0, (N, N))
    A = (A + A.T) / 2.0
    np.fill_diagonal(A, 1.0)
    
    mask_off_diagonal = ~np.eye(N, dtype=bool)
    A[mask_off_diagonal & (np.abs(A) > tau)] = 0.0

    non_zeros = np.count_nonzero(A)
    print(f"Número de elementos não-zeros na matriz A: {non_zeros}")

    # 2. Inicialização dos vetores
    x = np.zeros(N)
    b = np.random.uniform(-1.0, 1.0, N)

    # 3. Algoritmo do Slide #23
    r = b - np.dot(A, x)
    v = r.copy()

    # Listas para guardar o histórico do gráfico
    iterations_history = [0]
    residue_history = [np.linalg.norm(r)]

    print("\nIteração\t||r_n||")
    print("-" * 36)
    print(f"0\t\t{residue_history[0]:.4e}")

    for it in range(1, MAX_ITER + 1):
        Av = np.dot(A, v)
        r_dot_r = np.dot(r, r)
        v_dot_Av = np.dot(v, Av)

        if abs(v_dot_Av) < 1e-15:
            print("Aviso: v^T * A * v próximo de zero. Abortando iterações.")
            # Mesmo abortando, vamos repetir o último resíduo até o fim 
            # para o gráfico ficar contínuo igual ao do livro
            while len(residue_history) <= MAX_ITER:
                iterations_history.append(len(residue_history))
                residue_history.append(residue_history[-1])
            break

        alpha = r_dot_r / v_dot_Av
        x = x + alpha * v
        r_next = r - alpha * Av

        v_dot_A_rnext = np.dot(Av, r_next)
        beta = -v_dot_A_rnext / v_dot_Av

        v = r_next + beta * v
        r = r_next.copy()

        current_norm = np.linalg.norm(r)
        iterations_history.append(it)
        residue_history.append(current_norm)

        print(f"{it}\t\t{current_norm:.4e}")
        
    return iterations_history, residue_history

# Executa para os 4 casos de tau e armazena os resultados para o gráfico
taus = [0.01, 0.05, 0.1, 0.2]
results = {}

for t in taus:
    iters, residues = run_cg_experiment(t)
    results[t] = (iters, residues)

# --- CONSTRUÇÃO DO GRÁFICO (IDÊNTICO AO LIVRO) ---
plt.figure(figsize=(9, 6))

# Plota a linha de cada tau com marcadores
markers = ['o', 's', '^', 'd']
for i, t in enumerate(taus):
    iters, residues = results[t]
    plt.plot(iters, residues, label=f'$\\tau = {t}$', marker=markers[i], markersize=4, linewidth=1.2)

# Configurações dos Eixos
plt.yscale('log') # Escala Logarítmica no eixo Y (crucial para ver a convergência)
plt.xlim(0, MAX_ITER)
plt.xticks(range(0, MAX_ITER + 1, 2)) # Marcas de 2 em 2 iterações

# Limites do eixo Y baseados no comportamento do livro (de 10^-16 até 10^4)
plt.ylim(1e-16, 1e4) 

# Legendas e Títulos
plt.xlabel('n (Iterações)', fontsize=11)
plt.ylabel('$||r_n||$ (Norma do Resíduo)', fontsize=11)
plt.title('Curvas de Convergência do Gradiente Conjugado ($500 \\times 500$ Sparse Matrix)', fontsize=12, fontweight='bold')
plt.grid(True, which="both", ls="--", linewidth=0.5) # Linhas de grade
plt.legend(fontsize=10, loc='upper right')

# Salva a imagem automaticamente
plt.savefig('convergencia_cg.png', dpi=300, bbox_inches='tight')
print("\n[SUCESSO] Gráfico salvo como 'convergencia_cg.png'!")

# Mostra o gráfico na tela
plt.show()
