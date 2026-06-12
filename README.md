
# 🧮 Método dos Gradientes Conjugados

## 📌 Sobre o Projeto

Este repositório contém a implementação e a análise de convergência do **Método dos Gradientes Conjugados (CG)** aplicado à resolução de sistemas lineares de grande porte ($500 \times 500$).

O foco principal do estudo é avaliar o impacto do condicionamento da matriz dos coeficientes na velocidade e na estabilidade da convergência do método iterativo. Essa avaliação é feita controlando a esparsidade e as propriedades da matriz através da variação de um parâmetro de limiarização $\tau$ (tau).

Projeto prático desenvolvido para a disciplina de **Álgebra Linear Computacional**.

---

## 🎯 Objetivos Principais

* **Implementação Algorítmica:** Desenvolver o solver do Gradiente Conjugado passo a passo para sistemas da forma $Ax = b$, baseado na formulação clássica iterativa.
* **Geração de Matrizes:** Criar matrizes simétricas $500 \times 500$ com diagonais estritamente positivas e aplicar uma máscara de corte (limiar $\tau$) para os elementos fora da diagonal principal.
* **Análise de Convergência:** Investigar a correlação direta entre o parâmetro $\tau \in \{0.01, 0.05, 0.1, 0.2\}$ e o comportamento da norma do resíduo ( $||r_n||$ ).
* **Avaliação de Desempenho Visual:** Gerar um gráfico comparativo em escala logarítmica com a curva de convergência para cada valor de limiarização durante 20 iterações.

---

## 🛠️ Tecnologias Utilizadas

* **Linguagem Principal:** Python 3
* **Computação Numérica:** `NumPy` (para operações vetorizadas eficientes, cálculo de normas de Frobenius e produtos internos)
* **Visualização de Dados:** `Matplotlib` (para construção da curva de convergência em escala logarítmica)

---

## 🚀 Como Executar

**1. Clone o repositório**

```bash
git clone https://github.com/DanielySoares/Gradientes-Conjugados.git
cd Gradientes-Conjugados

```

**2. Instale as dependências requeridas**

```bash
pip install numpy matplotlib

```

**3. Execute a análise**

```bash
python main.py

```

*(Certifique-se de substituir `main.py` pelo nome correto do seu arquivo script).*

---

## 📊 Estrutura e Resultados

O script principal realiza as seguintes etapas operacionais de forma automática:

1. Fixa a *seed* para garantia de reprodutibilidade dos dados.
2. Executa as simulações para os 4 casos de teste de $\tau$ ($0.01, 0.05, 0.1$ e $0.2$).
3. Monitora se o denominador $\langle v, Av \rangle$ se aproxima de zero para evitar instabilidade por divisão zero.
4. Gera e salva localmente a imagem **`convergencia_cg.png`**, retratando de forma explícita que valores menores de $\tau$ (que preservam a dominância diagonal) levam a um declínio suave e rápido do resíduo para próximo de $10^{-10}$, enquanto valores maiores resultam em matrizes mal-condicionadas e consequente não-convergência nos primeiros 20 passos.

---