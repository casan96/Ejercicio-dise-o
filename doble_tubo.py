import math

def reynolds_f(mp,diametro,miu):
    return ( 4 * mp ) / ( math.pi * diametro * miu )

def dmlt_f(dT1, dT2):
    return (dT1 - dT2) / (math.log(dT1/dT2))

def Q_f(mp,cp,dt):
    return mp*cp*dt

def G_f(mp,area):
    return mp/area

def coef_pelicula(reynolds, prandalt, k, diametro):
    return 0.027 * math.pow(reynolds,0.8) * math.pow(prandalt,1/3) * k / diametro

def coef_global_limpio(coef_pelicula_i, coef_pelicula_o, area_i, area_o, diametro_ii, diametro_oi, l, k):
    return (1 / (coef_pelicula_i * area_i)) + (math.log(diametro_ii/diametro_oi)/(2*math.pi * l * k)) + (1 / (coef_pelicula_o * area_o))

def area_transferencia(diametro, l):
    return math.pi * diametro * l

def area_circun(diametro):
    return math.pow(diametro,2)*math.pi/4

def factor_friccion(reynolds):
    return 0.0035 + 0.264/math.pow(reynolds, 0.42)

def U_c(hi,ho):
    return hi * ho / (hi + ho)

diametro_interno_i = 3.068 * (25.4e-3)
diametro_interno_o = 3.5 * (25.4e-3)
diametro_anulo_i = 5.047 * (25.4e-3)

flujo_masico_i = 100000 / 3600
flujo_masico_a = 200000 / 3600

miu_i = 4.93e-3
miu_a = 7.55e-4

cp_etil = 2616
cp_agua = 4176

densidad_etil = 1081.45
densidad_agua = 993.5

k_i = 0.26
k_a = 0.6245

pr_i = 49.6
pr_a = 5

t_i1 = 97
t_i2 = 47
t_a1 = 25
t_a2 = 40

q_i = Q_f(flujo_masico_i,cp_etil, t_i1 - t_i2)
q_a = Q_f(flujo_masico_a, cp_agua, t_a2 - t_a1)
print(f"Q de etil.. = {q_i:.2f} \nQ de agua = {q_a:.2f}")

dmlt = dmlt_f(t_i1-t_a2, t_i2 - t_a1)
print(f"dmlt={dmlt:.2f} {diametro_interno_o} {diametro_anulo_i}")

reynolds_i = reynolds_f(flujo_masico_i, diametro_interno_o, miu_i)
reynolds_a = reynolds_f(flujo_masico_a, diametro_anulo_i - diametro_interno_o, miu_a)
print(f"Rei = {reynolds_i:.2f}\nReo = {reynolds_a:.2f}")

hi = coef_pelicula(reynolds_i, pr_i, k_i, diametro_interno_i)
ho = coef_pelicula(reynolds_a, pr_a, k_a, (math.pow(diametro_anulo_i,2) -math.pow(diametro_interno_o, 2)) / diametro_interno_o)
print(f"hi = {hi:.2f}\nho = {ho:.2f}")

uc = U_c(hi,ho)
print(f"Uc = {uc:.2f}")

rd = 0.0025
a = q_i/((uc / (1 + rd * uc)) * dmlt)
print(f"area diseño = {a:.2f}")
print(f"longitud de diseño {(a / 0.0929)/0.917*0.3048:.2f}")
#print(0.614 + (0.635 - 0.614)/(40-25)*(32.5-25))
