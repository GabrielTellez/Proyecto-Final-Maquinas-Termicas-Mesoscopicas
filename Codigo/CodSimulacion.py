import numpy as np
import matplotlib.pyplot as plt

#Evolucion de la posicion a partir de la ec de Langevin
def ev_posicion(Delta_s, kappa, theta, x_old):
    Delta_W = np.random.normal(loc=0, scale=1)    
    x_new = x_old*(1 - kappa*Delta_s) + np.sqrt(2*theta*Delta_s)*Delta_W
    return x_new

#Proceso termodinamico
def proceso(l, Delta_s, kappa, theta, x0):
    x = np.zeros(l)
    x[0] = x0

    if(x0 == 1):
        for i in range(100):
            x[0] = ev_posicion(Delta_s, kappa[0], theta[0], x[0])

    for i in range(l - 1):
        x[i + 1] = ev_posicion(Delta_s, kappa[i + 1], theta[i + 1], x[i])
        
    return x


#Protocolo 1

#Proceso isotermico protocolo 1
def isoterma_protocolo_1(l, s_0, s_f, theta):
    s = np.linspace(s_0, s_f, l)
    
    kappa = (1 - 2*s/st)**2*(kappa_A - kappa_C) + kappa_C
    theta = np.ones(l)*theta
    y = theta/kappa
    
    return [s, kappa, theta, y]

#Proceso adiabatico protocolo 1
def adiabata_protocolo_1(l, s_0, s_f, theta_0, theta_f):
    s = np.linspace(s_0, s_f, l)
    
    kappa = (1 - 2*s/st)**2*(kappa_A - kappa_C) + kappa_C
    theta = (theta_f - theta_0)*(s - s_0)/(s_f - s_0) + theta_0
    y = theta/kappa
    
    return [s, kappa, theta, y]

#Ciclo protocolo 1
def protocolo_1(l, cant_iter, kappa_A, y_A, kappa_B, y_B, kappa_C, y_C, kappa_D, y_D, theta_AB, theta_CD):
    sA = 0
    sB = st/2*(1 - np.sqrt((kappa_B - kappa_C)/(kappa_A - kappa_C)))
    sC = st/2
    sD = st/2*(1 + np.sqrt((kappa_D - kappa_C)/(kappa_A - kappa_C)))

    Delta_s1 = (sB - sA)/l
    Delta_s2 = (sC - sB)/l
    Delta_s3 = (sD - sC)/l
    Delta_s4 = (st - sD)/l
    
    #Compresion isotermica
    var_1 = isoterma_protocolo_1(l, sA, sB, theta_AB)
    #Compresion adiabática
    var_2 = adiabata_protocolo_1(l, sB, sC, theta_AB, theta_CD)
    #Expansion isotermica
    var_3 = isoterma_protocolo_1(l, sC, sD, theta_CD)
    #Expansion adiabatica
    var_4 = adiabata_protocolo_1(l, sD, st, theta_CD, theta_AB)

    #[s, kappa, theta, y]
    variables = np.array([var_1, var_2, var_3, var_4])
    Deltas_s = [Delta_s1, Delta_s2, Delta_s3, Delta_s4]
    
    y_sto = []
    W_sto = []

    x0 = np.ones(cant_iter)*np.sqrt(variables[0][3, 0])

    for i in range(4):
        var = variables[i]
        y_proc = np.zeros(l)
        W_proc_stochastic = []

        for j in range(cant_iter):
            kappa_proc = var[1]
            theta_proc = var[2]

            x = proceso(l, Deltas_s[i], kappa_proc, theta_proc, x0[j])    

            U = (kappa_proc/2)*(x)**2
            W = 0
            x0[j] = x[-1]
                
            for k in range(1, l):
                W += (U[k] - U[k - 1])
                
            W_proc_stochastic.append(W)
            y_proc += (x)**2
            
        y_proc /= cant_iter

        y_sto.append(y_proc)
        W_sto.append(W_proc_stochastic)
    
    W_sto = np.array(W_sto)
    return variables, y_sto, W_sto


#Protocolo 2

#Proceso isotermico protocolo 2
def isoterma_protocolo_2(l, Delta_s, kappa_0, y_0, kappa_f, y_f, theta):
    s_f = l*Delta_s
    s = np.linspace(0, s_f, l)
    
    kappa = (kappa_0*kappa_f*s_f**2)/(np.sqrt(kappa_f)*s_f + (np.sqrt(kappa_0) - np.sqrt(kappa_f))*s)**2
    theta = np.ones(l)*theta
    y = theta/kappa
    
    return [s, kappa, theta, y]

#Proceso adiabatico protocolo 2
def adiabata_protocolo_2(l, Delta_s, kappa_0, y_0, theta_0, kappa_f, y_f, theta_f):
    s_f = l*Delta_s
    s = np.linspace(0, s_f, l)
    
    alpha = ((theta_f**2/kappa_f) - (theta_0**2/kappa_0))*(s/s_f) + (theta_0**2/kappa_0)
    theta = (theta_0*theta_f * s_f**2)/(np.sqrt(theta_f)*s_f + (np.sqrt(theta_0) - np.sqrt(theta_f))*s)**2
    kappa = theta**2/alpha
    y = theta/kappa
    
    return [s, kappa, theta, y]

#Ciclo protocolo 2
def protocolo_2(l, cant_iter, kappa_A, y_A, kappa_B, y_B, kappa_C, y_C, kappa_D, y_D, theta_AB, theta_CD):
    Delta_s1 = st/(4*l)
    Delta_s2 = st/(4*l)
    Delta_s3 = st/(4*l)
    Delta_s4 = st/(4*l)

    #Compresion isotermica
    var_1 = isoterma_protocolo_2(l, Delta_s1, kappa_A, y_A, kappa_B, y_B, theta_AB)
    #Compresion adiabática
    var_2 = adiabata_protocolo_2(l, Delta_s2, kappa_B, y_B, theta_AB, kappa_C, y_C, theta_CD)
    #Expansion isotermica
    var_3 = isoterma_protocolo_2(l, Delta_s3, kappa_C, y_C, kappa_D, y_D, theta_CD)
    #Expansion adiabatica
    var_4 = adiabata_protocolo_2(l, Delta_s4, kappa_D, y_D, theta_CD, kappa_A, y_A, theta_AB)

    #[s, kappa, theta, y]
    variables = np.array([var_1, var_2, var_3, var_4])
    Deltas_s = [Delta_s1, Delta_s2, Delta_s3, Delta_s4]
    y_sto = []
    W_sto = []

    x0 = np.ones(cant_iter)*np.sqrt(variables[0][3, 0])

    for i in range(4):
        var = variables[i]
        y_proc = np.zeros(l)
        W_proc_stochastic = []

        for j in range(cant_iter):
            kappa_proc = var[1]
            theta_proc = var[2]

            x = proceso(l, Deltas_s[i], kappa_proc, theta_proc, x0[j])  

            U = (kappa_proc/2)*(x)**2
            W = 0
            x0[j] = x[-1]
                
            for k in range(1, l):
                W += (U[k] - U[k - 1])
                
            W_proc_stochastic.append(W)
            y_proc += (x)**2
                
        y_proc /= cant_iter

        y_sto.append(y_proc)
        W_sto.append(W_proc_stochastic)
    
    return variables, y_sto, W_sto


#Protocolo 3

#Proceso isotermico protocolo 3
def isoterma_protocolo_3(l, Delta_s, kappa_0, y_0, kappa_f, y_f, theta_0):
    s_f = l*Delta_s
    s = np.linspace(0, s_f, l)
    
    theta = np.ones(l)*theta_0
    y = (np.sqrt(y_0) + (np.sqrt(y_f) - np.sqrt(y_0))*s/s_f)**2
    kappa = theta/y - (np.sqrt(y_f/y) - np.sqrt(y_0/y))/s_f

    kappa[0] = kappa_0
    kappa[-1] = kappa_f
    
    return [s, kappa, theta, y]

#Proceso adiabatico protocolo 3
def adiabata_protocolo_3(l, Delta_s, kappa_0, y_0, theta_0, kappa_f, y_f, theta_f):
    s_f = l*Delta_s
    s = np.linspace(0, s_f, l)
    
    y = y_0 + (y_f - y_0)*s/s_f
    theta = (y_0*theta_0 + (y_f*theta_f - y_0*theta_0)*s/s_f)/y
    kappa = theta/y - (y_f - y_0)/(2*y*s_f)

    kappa[0] = kappa_0
    kappa[-1] = kappa_f
    
    return [s, kappa, theta, y]

#Ciclo protocolo 3
def protocolo_3(l, cant_iter, kappa_A, y_A, kappa_B, y_B, kappa_C, y_C, kappa_D, y_D, theta_AB, theta_CD):
    sf_BC = (y_C - y_B)**2/(2*(y_C*theta_CD - y_B*theta_AB))
    sf_DA = (y_A - y_D)**2/(2*(y_A*theta_AB - y_D*theta_CD))

    sf_iso = (st - sf_BC - sf_DA)/2

    Delta_s1 = sf_iso/l
    Delta_s2 = sf_BC/l
    Delta_s3 = sf_iso/l
    Delta_s4 = sf_DA/l
    
    #Compresion isotermica
    var_1 = isoterma_protocolo_3(l, Delta_s1, kappa_A, y_A, kappa_B, y_B, theta_AB)
    #Compresion adiabática
    var_2 = adiabata_protocolo_3(l, Delta_s2, kappa_B, y_B, theta_AB, kappa_C, y_C, theta_CD)
    #Expansion isotermica
    var_3 = isoterma_protocolo_3(l, Delta_s3, kappa_C, y_C, kappa_D, y_D, theta_CD)
    #Expansion adiabatica
    var_4 = adiabata_protocolo_3(l, Delta_s4, kappa_D, y_D, theta_CD, kappa_A, y_A, theta_AB)

    variables = np.array([var_1, var_2, var_3, var_4])
    Deltas_s = [Delta_s1, Delta_s2, Delta_s3, Delta_s4]
    y_sto = []
    W_sto = []

    x0 = np.ones(cant_iter)*np.sqrt(variables[0][3, 0])

    for i in range(4):
        var = variables[i]
        y_proc = np.zeros(l)
        W_proc_stochastic = []

        for j in range(cant_iter):
            kappa_proc = var[1]
            theta_proc = var[2]

            x = proceso(l, Deltas_s[i], kappa_proc, theta_proc, x0[j])    

            U = (kappa_proc/2)*(x)**2
            W = 0
            x0[j] = x[-1]
                
            for k in range(1, l):
                W += (U[k] - U[k - 1])
                
            W_proc_stochastic.append(W)
            y_proc += (x)**2
                
        y_proc /= cant_iter

        y_sto.append(y_proc)
        W_sto.append(W_proc_stochastic)

    return variables, y_sto, W_sto


#Graficas
def graficas(variables, y_sto, W_sto, num):

    titulos = ["Compresión isotérmica", "Compresión adiabática", "Expansión isotérmica", "Expansión adiabática"]
    name_arch = []
    name_arch.append(["FigAB-1", "FigAB-2"])
    name_arch.append(["FigBC-1", "FigBC-2"])
    name_arch.append(["FigCD-1", "FigCD-2"])
    name_arch.append(["FigDA-1", "FigDA-2"])

    for i in range(4):
        kappa_proc = variables[i, 1]
        y_proc = variables[i, 3]
        y_sto_proc = y_sto[i]
        W_sto_proc = W_sto[i]

        plt.figure()
        plt.scatter(kappa_proc[0:-1:200], y_sto_proc[0:-1:200], c="red", label="Simulado")
        plt.plot(kappa_proc, y_proc, label="Esperado")
        plt.xlabel("$\kappa$")
        plt.ylabel("y")
        plt.legend()
        plt.title(titulos[i] + " P{}".format(num))
        plt.savefig(name_arch[i][0] + " P{}.png".format(num))
        plt.close()

        W_prom = np.mean(W_sto_proc)
        W_std = np.std(W_sto_proc)

        res_prom = " = {:.3f} $\pm$ {:.3f}".format(W_prom, W_std)
        
        plt.figure()
        plt.hist(W_sto_proc, bins=40, density=True)
        plt.axvline(W_prom, c="red", label="Promedio: $\mathcal{W}$" + res_prom)
        plt.title("Distribución trabajo estocástico Proceso {} P{}".format(i + 1, num))
        plt.xlabel("$\mathcal{W}_{\mathrm{sto}}$")
        plt.ylabel("Probabilidad")
        plt.legend()
        plt.savefig(name_arch[i][1] + " P{}.png".format(num))
        plt.close()

    #Ciclo

    plt.figure()

    for i in range(2):
        kappa_1 = variables[2*i, 1]
        kappa_2 = variables[2*i + 1, 1]
        y_proc_1 = variables[2*i, 3]
        y_proc_2 = variables[2*i + 1, 3]
        y_sto_proc_1 = y_sto[2*i]
        y_sto_proc_2 = y_sto[2*i + 1]

        plt.scatter(kappa_1[0:-1:200], y_sto_proc_1[0:-1:200], c="red", alpha=0.5)
        plt.scatter(kappa_2[0:-1:200], y_sto_proc_2[0:-1:200], c="green", alpha=0.5)

        plt.plot(kappa_1, y_proc_1, c="red")
        plt.plot(kappa_2, y_proc_2, c="green")

    plt.xlabel("$\kappa$")
    plt.ylabel("y")
    plt.title("Ciclo irreversible Protocolo {}".format(num))
    plt.savefig("FigCiclo-1P{}.png".format(num))
    plt.close()

    W_ciclo_sto = W_sto[0] + W_sto[1] + W_sto[2] + W_sto[3]
    W_ciclo_prom = np.mean(W_ciclo_sto)
    W_ciclo_std = np.std(W_ciclo_sto)

    res_ciclo_prom = " = {:.3f} $\pm$ {:.3f}".format(W_ciclo_prom, W_ciclo_std)

    plt.figure()
    plt.hist(W_ciclo_sto, bins=40, density=True)
    plt.axvline(W_ciclo_prom, c="red", label="Promedio: $\mathcal{W}$" + res_ciclo_prom)
    plt.title("Distribución trabajo estocástico del ciclo P{}".format(num))
    plt.xlabel("$\mathcal{W}_{\mathrm{sto}}$")
    plt.ylabel("Probabilidad")
    plt.legend()
    plt.savefig("FigCiclo-2P{}.png".format(num))
    plt.close()
    return None


#Constantes y parametros
chi = 0.6
nu = 0.6
c = 0.96
d = 1.03

l = 5000
cant_iter = 1000
st = 25

#Puntos A, B, C, D

kappa_A = 1
y_A = 1
kappa_B = chi
y_B = 1/chi
kappa_C = c*(nu**2)*chi
y_C = (c*nu*chi)**(-1)
kappa_D = d*nu**2
y_D = (d*nu)**(-1)

theta_AB = 1
theta_CD = nu

variables, y_sto, W_sto = protocolo_1(l, cant_iter, kappa_A, y_A, kappa_B, y_B, kappa_C, y_C, kappa_D, y_D, theta_AB, theta_CD)
graficas(variables, y_sto, W_sto, 1)

variables, y_sto, W_sto = protocolo_2(l, cant_iter, kappa_A, y_A, kappa_B, y_B, kappa_C, y_C, kappa_D, y_D, theta_AB, theta_CD)
graficas(variables, y_sto, W_sto, 2)

variables, y_sto, W_sto = protocolo_3(l, cant_iter, kappa_A, y_A, kappa_B, y_B, kappa_C, y_C, kappa_D, y_D, theta_AB, theta_CD)
graficas(variables, y_sto, W_sto, 3)
