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

def k_etil(T):
    return 0.088067 + 9.4712e-4*(T) - 1.3114e-6*(T**2)

def k_agua(T):
    return -0.432 + 0.0057255*(T) - 8.078e-6*(T**2) + 1.861e-9*(T**3)

def p_etil(T):
    return (1.315/(0.25125**(1+((1-T/720)**0.21868)))) * 62.068

def p_agua(T):
    tau = 1- T/647.096
    return (17.863  + 58.606*(tau**0.35) - 95.396*(tau**(2/3)) + 213.89*tau -141.26*(tau**(4/3))) * 18.015

def cp_etil(T):
    return (35540 + 436.78 * T - 0.18486 * (T**2)) / 62.068

def cp_agua(T):
    return (276370 - 2090.1 * T + 8.125 * (T**2) - 0.014116 * (T**3) + 9.3701e-6 * (T**4)) / 18.015

def miu_etil(T):
    return math.exp(-20.515 + 2468.5/T + 1.2435* math.log(T) + 2.4998e12*(T**(-5)))

def miu_agua(T):
    return math.exp(-52.843 + 3703.6/T + 5.866* math.log(T) - 5.879e-29*(T**10))

def main():
    #Propiedades

    #etilenglicol
    # °C 
    T_avg_s = 72 
    # j/kg.°C 
    cp_s = cp_etil(T_avg_s + 273.15) 
    # kg/m3 
    p_s = p_etil(T_avg_s + 273.15)
    # kg/m.s
    miu_s = miu_etil(T_avg_s + 273.15)
    # W/m.°C
    k_s = k_etil(T_avg_s + 273.15)
    print(f"k {k_s}")
    pr_s = cp_s * miu_s / k_s
    print(f"pr {pr_s}")

    #agua
    T_avg_t = 32.5 
    cp_t = cp_agua(T_avg_t + 273.15)
    p_t = p_agua(T_avg_t + 273.15)
    miu_t = miu_agua(T_avg_t + 273.15) 
    k_t = k_agua(T_avg_t + 273.15)
    pr_t = cp_t * miu_t / k_t

    #diametro tubo
    dt = 1 * (25.4e-3)
    print(f"do={dt:.5f}")

    #diametro inteno tubo BWG 14
    d_i = 0.834 * (25.4e-3)
    print(f"di={d_i:.5f}")

    #pitch 
    pt = 1.25 * (25.4e-3)
    print(f"pt={pt:.5f}")

    #separacion de los tubos, pitch - diametro tubo
    cpr = pt - dt
    print(f"C={cpr:.5f}")

    #diametro interno shell
    d_s = 19.25 * (25.4e-3)
    print(f"ds={d_s:.5f}")

    #separacion de baffle
    b = 0.18

    #area shell
    a_s = area_shell(d_s,pt,cpr, b)
    print(f"as={a_s:.5f}")

    #flujo masico shell
    mp_s = 100000 / 3600
    print(f"mps={mp_s:.5f}")

    #velocidad de masa shell
    G_s = G(mp_s, a_s)
    print(f"Gs={G_s:.5f}")

    #diametro equivalente shell
    d_e = d_equivalente(pt, dt)
    print(f"De={d_e:.5f}")

    #reynolds shell
    re_s = reynolds(G_s,d_e, miu_s)
    print(f"Res={re_s:.5f}")

    t_w = tw(25,40,97,47)
    print(f"tw = {t_w:.2f}")

    miu_w = miu_etil(t_w + 273.15)
    print(f"miuw={miu_w:.5f}")

    #coeficiente de pelicula h0
    h0 = coef_pelicula(re_s, pr_s, k_s, d_e, miu_s, miu_w)
    print(f"ho={h0:.5f}")

    #numero de tubos
    Nt = 132

    #pasos por los tuvos
    n = 2

    #area total de los tubos
    a_t = atp(Nt,d_i)
    print(f"at={a_t:.5f}")

    #flujo masico tubos
    mp_t = 200000 / 3600
    print(f"mpt={mp_t:.5f}")

    vm_t = mp_t/ (p_t * a_t)
    print(f"vmp={vm_t:.5f}")

    re_t = p_t * vm_t * d_i / miu_t
    print(f"Ret={re_t:.5f}")

    f_t = (1.58 * math.log(re_t)-3.28)**(-2)
    print(f"ft={f_t:.5f}")

    nu_t = (f_t/2)*(re_t-1000)*pr_t/(1 + 12.7*((f_t/2)**(1/2))*((pr_t)**(2/3) - 1))
    print(f"Nu={nu_t:.5f}")

    hi = nu_t * k_t / d_i
    print(f"hi={hi:.5f}")

    uf = 1/(dt/(d_i*hi) + dt*0.0003/d_i + dt*math.log(dt/d_i)/(2*45)+0.0002+1/h0)
    print(f"uc {uf}")

    l = 4
    N_b = l/b-1
    print(f"Nb={N_b:.5f}")
    f = math.exp(0.576 -0.19*math.log(re_s))
    print(f"f={N_b:.5f}")
    dp_s = f *(G_s**2) * (N_b + 1) * d_s / (2 * p_s * d_e * ((miu_s / miu_w) ** 0.14))
    print(f"dp = {dp_s / 6894.75729:.4f}")

    Q_s = mp_s * cp_s * (40 - 25)
    print(f"Qs={Q_s:.5f}")
    dlt = dmlt(97-40,47-25)
    print(f"dlmt={dlt:.5f}")
    r = R(97,47,25,40)
    print(f"R={r:.5f}")
    p = P(97,25,40)
    print(f"P={p:.5f}")
    Fr = Ft(r,p)
    print(f"Fr={Fr:.5f}")
    a_of = Q_s/(uf * Fr * dlt)
    print(f"Aof {a_of:.4f}")
    L_r = a_of/(math.pi * dt * Nt)
    print(f"lr = {L_r:.6f}")

    dp_t = (4 * f_t * (l * n)/d_i + 4 * n)* p_t * (vm_t**2)/2
    print(f"dp t = {dp_t / 6894.75729:.4f}")



main()
# print(f"{d_equivalente(1.25, 1)/12}")
