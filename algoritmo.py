import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')


# Parâmetros do problema
T0 = 90.0  # Temperatura inicial (°C)
Ta = 25.0  # Temperatura ambiente (°C)
C = T0 - Ta  # Constante da solução (T0 - Ta)
k = 0.1    # Constante de resfriamento (s^-1)
t_max = 30.0  # Tempo máximo (s)
t_pontos = [0, 10, 20, 30] 

# Solução analítica: T(t) = (T0 - Ta) * exp(-k * t) + Ta
def solucao_analitica(t):
    return C * np.exp(-k * t) + Ta

# Validar solução analítica
print("Validação da Solução Analítica (T(t) = 65 * exp(-0.1 * t) + 25):")
for t in t_pontos:
    T_ana = solucao_analitica(t)
    print(f"t = {t:<2}: T_ana = {T_ana:.8f}")

# Método de Euler para tempos específicos
def metodo_euler_pontos(t_pontos, h, T0, k, Ta):
    T = [T0]  # Garante que o primeiro ponto é exatamente a condição inicial
    t_prev = t_pontos[0]
    T_prev = T0
    
    for t_target in t_pontos[1:]:
        # Calcula quantos passos são necessários para chegar ao próximo ponto
        steps = int(np.ceil((t_target - t_prev) / h))
        adjusted_h = (t_target - t_prev) / steps 
        
        for _ in range(steps):
            T_next = T_prev + adjusted_h * (-k * (T_prev - Ta))
            t_prev += adjusted_h
            T_prev = T_next
        
        T.append(T_prev)
    
    return np.array(t_pontos), np.array(T)

# Método de Runge-Kutta de 4ª ordem para tempos específicos
def runge_kutta_4_pontos(t_pontos, h, T0, k, Ta):
    T = [T0]  
    t_prev = t_pontos[0]
    T_prev = T0
    
    for t_target in t_pontos[1:]:
        # Calcula quantos passos são necessários para chegar ao próximo ponto
        steps = int(np.ceil((t_target - t_prev) / h))
        adjusted_h = (t_target - t_prev) / steps
        
        for _ in range(steps):
            k1 = -k * (T_prev - Ta)
            k2 = -k * (T_prev + adjusted_h/2 * k1 - Ta)
            k3 = -k * (T_prev + adjusted_h/2 * k2 - Ta)
            k4 = -k * (T_prev + adjusted_h * k3 - Ta)
            T_next = T_prev + (adjusted_h/6) * (k1 + 2*k2 + 2*k3 + k4)
            t_prev += adjusted_h
            T_prev = T_next
        
        T.append(T_prev)
    
    return np.array(t_pontos), np.array(T)

# Método de Euler para gráficos
def metodo_euler(h, t_max, T0, k, Ta):
    t = np.arange(0, t_max + h, h)
    T = np.zeros(len(t))
    T[0] = T0
    for i in range(1, len(t)):
        T[i] = T[i-1] + h * (-k * (T[i-1] - Ta))
    return t, T

# Método de Runge-Kutta de 4ª ordem para gráficos
def runge_kutta_4(h, t_max, T0, k, Ta):
    t = np.arange(0, t_max + h, h)
    T = np.zeros(len(t))
    T[0] = T0
    for i in range(1, len(t)):
        k1 = -k * (T[i-1] - Ta)
        k2 = -k * (T[i-1] + h/2 * k1 - Ta)
        k3 = -k * (T[i-1] + h/2 * k2 - Ta)
        k4 = -k * (T[i-1] + h * k3 - Ta)
        T[i] = T[i-1] + (h/6) * (k1 + 2*k2 + 2*k3 + k4)
    return t, T

# Calcula o desvio relativo percentual
def desvio_relativo(T_num, T_ana):
    return 100 * np.abs((T_num - T_ana) / T_ana)

# Solução analítica para gráficos
t_analitico = np.linspace(0, t_max, 1000)
T_analitico = solucao_analitica(t_analitico)

# Calcula soluções numéricas para gráficos
t_euler_1, T_euler_1 = metodo_euler(1.0, t_max, T0, k, Ta)
t_euler_05, T_euler_05 = metodo_euler(0.5, t_max, T0, k, Ta)
t_euler_01, T_euler_01 = metodo_euler(0.1, t_max, T0, k, Ta)
t_rk4_2, T_rk4_2 = runge_kutta_4(2.0, t_max, T0, k, Ta)
t_rk4_1, T_rk4_1 = runge_kutta_4(1.0, t_max, T0, k, Ta)
t_rk4_05, T_rk4_05 = runge_kutta_4(0.5, t_max, T0, k, Ta)

# Gráfico de temperatura
plt.figure(figsize=(10, 6))
plt.plot(t_analitico, T_analitico, 'k-', label='Solução Analítica', linewidth=2)
plt.plot(t_euler_1, T_euler_1, 'r--', label='Euler h=1.0')
plt.plot(t_euler_05, T_euler_05, 'g--', label='Euler h=0.5')
plt.plot(t_euler_01, T_euler_01, 'b--', label='Euler h=0.1')
plt.plot(t_rk4_2, T_rk4_2, 'c-.', label='RK4 h=2.0')
plt.plot(t_rk4_1, T_rk4_1, 'm-.', label='RK4 h=1.0')
plt.plot(t_rk4_05, T_rk4_05, 'y-.', label='RK4 h=0.5')
plt.xlabel('Tempo (s)')
plt.ylabel('Temperatura (°C)')
plt.title('Resfriamento de Componente Eletrônico')
plt.legend()
plt.grid(True)
plt.savefig('grafico_resfriamento.png')
print("Gráfico de temperatura salvo como 'grafico_resfriamento.png'")
plt.show()
plt.close()

# Gráfico de erro absoluto (escala logarítmica)
plt.figure(figsize=(10, 6))
plt.semilogy(t_euler_1, np.abs(T_euler_1 - solucao_analitica(t_euler_1)), 'r--', label='Erro Euler h=1.0')
plt.semilogy(t_euler_05, np.abs(T_euler_05 - solucao_analitica(t_euler_05)), 'g--', label='Erro Euler h=0.5')
plt.semilogy(t_euler_01, np.abs(T_euler_01 - solucao_analitica(t_euler_01)), 'b--', label='Erro Euler h=0.1')
plt.semilogy(t_rk4_2, np.abs(T_rk4_2 - solucao_analitica(t_rk4_2)), 'c-.', label='Erro RK4 h=2.0')
plt.semilogy(t_rk4_1, np.abs(T_rk4_1 - solucao_analitica(t_rk4_1)), 'm-.', label='Erro RK4 h=1.0')
plt.semilogy(t_rk4_05, np.abs(T_rk4_05 - solucao_analitica(t_rk4_05)), 'y-.', label='Erro RK4 h=0.5')
plt.xlabel('Tempo (s)')
plt.ylabel('Erro Absoluto (°C)')
plt.title('Erro Absoluto das Soluções Numéricas')
plt.legend()
plt.grid(True)
plt.savefig('grafico_erro_absoluto.png')
print("Gráfico de erro absoluto salvo como 'grafico_erro_absoluto.png'")
plt.show()
plt.close()

# Calcula soluções numéricas nos tempos exatos
_, T_euler_1_pontos = metodo_euler_pontos(t_pontos, 1.0, T0, k, Ta)
_, T_euler_05_pontos = metodo_euler_pontos(t_pontos, 0.5, T0, k, Ta)
_, T_euler_01_pontos = metodo_euler_pontos(t_pontos, 0.1, T0, k, Ta)
_, T_rk4_2_pontos = runge_kutta_4_pontos(t_pontos, 2.0, T0, k, Ta)
_, T_rk4_1_pontos = runge_kutta_4_pontos(t_pontos, 1.0, T0, k, Ta)
_, T_rk4_05_pontos = runge_kutta_4_pontos(t_pontos, 0.5, T0, k, Ta)


# Verificação do tempo inicial
print("\nVerificação do tempo inicial t=0:")
print(f"Solução analítica: {solucao_analitica(0):.8f}")
print(f"Euler h=1.0: {T_euler_1_pontos[0]:.8f}")
print(f"RK4 h=2.0: {T_rk4_2_pontos[0]:.8f}")

# Imprime tabela de temperaturas para depuração
print("\nValores de Temperatura (para depuração):")
print(f"{'t (s)':<10} {'Analítica':<15} {'Euler h=1.0':<15} {'Euler h=0.5':<15} {'Euler h=0.1':<15} {'RK4 h=2.0':<15} {'RK4 h=1.0':<15} {'RK4 h=0.5':<15}")
for i, t in enumerate(t_pontos):
    T_ana = solucao_analitica(t)
    print(f"{t:<10} {T_ana:<15.8f} {T_euler_1_pontos[i]:<15.8f} {T_euler_05_pontos[i]:<15.8f} {T_euler_01_pontos[i]:<15.8f} "
          f"{T_rk4_2_pontos[i]:<15.8f} {T_rk4_1_pontos[i]:<15.8f} {T_rk4_05_pontos[i]:<15.8f}")

# Calcula e imprime tabela de desvios relativos
desvios = {
    't': t_pontos,
    'Euler h=1.0': [],
    'Euler h=0.5': [],
    'Euler h=0.1': [],
    'RK4 h=2.0': [],
    'RK4 h=1.0': [],
    'RK4 h=0.5': []
}

for i, t in enumerate(t_pontos):
    T_ana = solucao_analitica(t)
    desvios['Euler h=1.0'].append(desvio_relativo(T_euler_1_pontos[i], T_ana))
    desvios['Euler h=0.5'].append(desvio_relativo(T_euler_05_pontos[i], T_ana))
    desvios['Euler h=0.1'].append(desvio_relativo(T_euler_01_pontos[i], T_ana))
    desvios['RK4 h=2.0'].append(desvio_relativo(T_rk4_2_pontos[i], T_ana))
    desvios['RK4 h=1.0'].append(desvio_relativo(T_rk4_1_pontos[i], T_ana))
    desvios['RK4 h=0.5'].append(desvio_relativo(T_rk4_05_pontos[i], T_ana))

print("\nTabela de Desvio Relativo Percentual:")
print(f"{'t (s)':<10} {'Euler h=1.0':<15} {'Euler h=0.5':<15} {'Euler h=0.1':<15} {'RK4 h=2.0':<15} {'RK4 h=1.0':<15} {'RK4 h=0.5':<15}")
for i in range(len(t_pontos)):
    print(f"{t_pontos[i]:<10} {desvios['Euler h=1.0'][i]:<15.8f} {desvios['Euler h=0.5'][i]:<15.8f} {desvios['Euler h=0.1'][i]:<15.8f} "
          f"{desvios['RK4 h=2.0'][i]:<15.8f} {desvios['RK4 h=1.0'][i]:<15.8f} {desvios['RK4 h=0.5'][i]:<15.8f}")