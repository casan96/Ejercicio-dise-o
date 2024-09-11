import math

def Q(mp, cp, dt):
    return mp*cp*dt

def G(mp,area):
    return mp/area

def reynolds(g,d,miu):
    return d*g/miu

def coef_pelicula(reynolds, prandalt, k, diametro, miu, miu_w):
    return 0.36 * (reynolds**(0.55)) * (prandalt**(1/3)) * math.pow(miu/miu_w, 0.14) * k / diametro

def dmlt(dT1, dT2):
    return (dT2 - dT1) / (math.log(dT2/dT1))

def R(T1, T2, t1, t2):
    return (T1 - T2)/(t2 -t1)

def dmlt(dT1, dT2):
    return (dT1 - dT2) / (math.log(dT1/dT2))

def P(T1, t1, t2):
    return (t2 - t1)/(T1 - t1)

def Ft(R, P):
    return (math.sqrt(R**2 + 1)/(R-1))*math.log((1-P)/(1-R*P))/math.log(h_Ft(R,P))

def h_Ft(R, P): 
    return (2-P*(R + 1 - math.sqrt(R**2 + 1)))/(2-P*(R+1+math.sqrt(R**2 + 1)))

def area_shell(d, pt, c,b):
    return d*c*b/(pt)

def area_tubos(Nt,a_p,n):
    return Nt*a_p/(144*n)

def atp(Nt,d):
    return math.pi*(d**2)/4*Nt/2

def d_equivalente(pt, d0):
    return 4*((pt**2) - math.pi*(d0**2)/4)/(math.pi * d0)

def tw(Tc1,Tc2,Th1,Th2):
    return (1/2) * ((Tc1 + Tc2)/2 + (Th1 +Th2)/2)

def main():
    #Propiedades

    #etilenglicol
    # °C 
    T_avg_s = 72 
    # j/kg.°C 
    cp_s = 2606 
    # kg/m3 
    p_s = 1081.45
    # kg/m.s
    miu_s = 4.93e-3 
    # W/m.°C
    k_s = 0.26
    pr_s = 49.6

    #agua
    T_avg_t = 32.5 
    cp_t = 4176 
    p_t = 993.9 
    miu_t = 7.55e-4 
    k_t = 0.6245
    pr_t = 5

    #diametro tubo
    dt = 1 * (25.4e-3)

    #diametro inteno tubo
    d_i = 0.670 * (25.4e-3)

    #pitch 
    pt = 1.25 * (25.4e-3)

    #separacion de los tubos, pitch - diametro tubo
    cpr = pt - dt

    #diametro interno shell
    d_s = 17.25 * (25.4e-3)

    #separacion de baffle
    b = 0.2

    #area shell
    a_s = area_shell(d_s,pt,cpr, b)

    #flujo masico shell
    mp_s = 100000 / 3600

    #velocidad de masa shell
    G_s = G(mp_s, a_s)

    #diametro equivalente shell
    d_e = d_equivalente(pt, dt)

    #reynolds shell
    re_s = reynolds(G_s,d_e, miu_s)

    t_w = tw(25,40,97,47)
    print(f"tw = {t_w:.2f}")

    miu_w = 7.014e-3

    #coeficiente de pelicula h0
    h0 = coef_pelicula(re_s, pr_s, k_s, d_e, miu_s, miu_w)

    #numero de tubos
    Nt = 112

    #pasos por los tuvos
    n = 2

    #area tablas
    a_t_p = 0.3526 * (25.4e-3)**2

    #area total de los tubos
    a_t = atp(Nt,d_i)

    #flujo masico tubos
    mp_t = 200000 / 3600

    vm_t = mp_t/ (p_t * a_t)

    re_t = p_t * vm_t * d_i / miu_t

    f_t = (1.58 * math.log(re_t)-3.28)**(-2)

    nu_t = (f_t/2)*(re_t-1000)*pr_t/(1 + 12.7*((f_t/2)**(1/2))*((pr_t)**(2/3) - 1))

    hi = nu_t * k_t / d_i

    uf = 1/(dt/(d_i*hi) + dt*0.000176/d_i + dt*math.log(dt/d_i)/(2*60)+0.000176+1/h0)

    l = 4.5
    N_b = l/b-1
    f = math.exp(0.576 -0.19*math.log(re_s))
    dp_s = f *(G_s**2) * (N_b + 1) * d_s / (2 * p_s * d_e * ((miu_s / miu_w) ** 0.14))
    print(f"dp = {dp_s / 6894.75729:.4f}")

    Q_s = mp_s * cp_s * (40 - 25)
    dlt = dmlt(97-40,47-25)
    r = R(97,47,25,40)
    p = P(97,25,40)
    Fr = Ft(r,p)
    a_of = Q_s/(uf * Fr * dlt)
    print(f"Aof {a_of:.4f}")
    L_r = a_of/(math.pi * dt * Nt)
    print(f"l = {L_r:.6f}")

    dp_t = (4 * f_t * (l * n)/d_i + 4 * n)* p_t * (vm_t**2)/2
    print(f"dp t = {dp_t / 6894.75729:.4f}")



main()
# print(f"{d_equivalente(1.25, 1)/12}")
