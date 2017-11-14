# Simulação computacional do pendulo duplo por Igor Fígueredo

# primeiramente é necessário importar as bibliiotecas com as quais se deseja trabalhar.

from numpy import sin, cos
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation

# vamos pensar sobre a física do bagulho, primeiramente resolvemos a equação do pendulo duplo utilizando formalismo lagrangeano. A referência utilizada para conferir os cálculos foi: http://nfist.pt/~pqueiroz/IC/pendulo/pendulo.pdf


g= 9.8
l1= 3.0
l2=1.0
m1=1.0
m2=1.0


def derivadas (estado,t):
	dydx=np.zeros_like(estado)
	dydx[0]=estado[1]

	diferenca = estado[0]-estado[2]

	denominador1 = (m1+m2)*l1 - m2*l1*cos(diferenca)*cos(diferenca)
	
	dydx[1] = ((-m2*l1*estado[1]*estado[1]*cos(diferenca)*sin(diferenca)+
		    m2*g*sin(estado[2])*cos(diferenca)-
		    m2*l2*estado[3]*estado[3]*sin(diferenca)-
		    (m1+m2)*g*sin(estado[0])))/denominador1

	dydx[2]=estado[3]

	denominador2 = (l2/l1)*denominador1
	dydx[3] = (m2*l2*estado[3]*estado[3]*sin(diferenca)*cos(diferenca)+
		   (m1+m2)*g*sin(estado[0])*cos(diferenca)+
		   (m1+m2)*l1*estado[1]*estado[1]*sin(diferenca)-
		   (m1+m2)*g*sin(estado[2]))/denominador2

	return dydx 

# contantes:.;

#criando uma lista de condições iniciais:

dt=0.05
t=np.arange(0.0,60,dt)

th1=120.0 #Angulo 1 inicial	
w1=0.0 # velocidade angular 1 inicial
th2= 180.0 # angulo 2 inicial
w2= 0.0 # velocidade angular 2 inicial

# estados iniciais:
estado = np.radians([th1,w1,th2,w2])

# Esta é a parte complicada, temos que resolver a edo numericamente, para isso o python dispoe de diversas bibliotecas, a que utilizaremos hoje é a scipy.

y = integrate.odeint(derivadas,estado,t)

# derivadas, é a edo que será resolvida, para isso é necessário crirar uma função que retorne o valor das derivadas em x e y.

# sabemos que para todo instante temos:

x1= l1*sin(y[:,0])
y1=-l1*cos(y[:,0])


x2= l2*sin(y[:,2]) + x1
y2 = -l2*cos(y[:,2]) + y1

fig = plt.figure()
ax = fig.add_subplot(111, autoscale_on=False, xlim=(-4, 4), ylim=(-4, 4))
ax.grid()

line, = ax.plot([], [], 'o-', lw=2)
time_template = 'time = %.1fs'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)


def init():
    line.set_data([], [])
    time_text.set_text('')
    return line, time_text


def animate(i):
    thisx = [0, x1[i], x2[i]]
    thisy = [0, y1[i], y2[i]]

    line.set_data(thisx, thisy)
    time_text.set_text(time_template % (i*dt))
    return line, time_text

ani = animation.FuncAnimation(fig, animate, np.arange(1, len(y)),
                              interval=25, blit=True, init_func=init)

#ani.save('Pendulo_duplo.mp4', fps=15)
plt.show()

