#--------------------------------------------------------------------
# PV PANEL MAXIMUM POWER POINT SOLVER
#--------------------------------------------------------------------
# PV Panel MPP Solver
# Author: Alvaro Augusto Volpato
# Description: this script intends to find a PV panel's MPP by solving its nonlinear equations.

# Header: importing libraries ---------------------------------------
import math
import numpy as np
from scipy.optimize import fsolve
linspace = np.linspace
logspace = np.logspace
import matplotlib.pyplot as pyplot
import matplotlib as mplot

mplot.rcParams['text.usetex'] = True
mplot.rcParams['text.latex.unicode'] = True


exp = np.exp
sqrt = np.sqrt

# Definig the parameters --------------------------------------------
IrrSrc = 8.487; 		# Irradiation curent at SRC
I0Src = 6.33*10**-9;		# Diode reverse current at SRC
Rs = 5.125*10**-3;		# Cell series resistance
Rp = 5.837;			# Cell parallel resistance
n = 1.149;			# Diode ideality factor
Ns = 50;			# Number of cells in series
Np = 1;				# Number of cells in parallel
GSrc = 1000;			# Irradiance at SRC
Eg0 = 1.17*1.6*10**-19;		# Silicon bandgap energy
alpha = 7.01*10**-4*1.6*10**-19;# Bandgap energy factor alpha
beta = 1108;			# Bandgap energy factor beta
k = 1.38*10**-23;		# Boltzmanns constant
ThetaSrc = 25 + 273.15;		# Temperature at SRC
mu = 0.05/100*IrrSrc;			# Short Circuit current temperature factor
q = 1.6*10**-19;		# Electron fundamental charge

Rs *= Ns/Np;			# Making the resistance transformation
Rp *= Ns/Np;			# Idem


# Defining the system nonlinear function ----------------------------

def FSim(V,I,G,Theta):		# Defining simulation nonlinear function
	return -I + Np*G/GSrc*(IrrSrc + mu*(Theta - ThetaSrc)) - Np*I0Src*(ThetaSrc/Theta)**(3/n)*exp(q/(n*k)*(Eg0 - alpha*Theta**2/(Theta + beta))*(1/ThetaSrc - 1/Theta))*( exp(q*(V + Rs*I)/(Ns*n*k*Theta)) - 1 ) - (V + Rs*I)/Rp;

def FMpp(z,G,Theta):		# Defining nonlinear function to find MPP
	V = z[0]
	I = z[1]

	F0 = -I + Np*G/GSrc*(IrrSrc + mu*(Theta - ThetaSrc)) - Np*I0Src*(ThetaSrc/Theta)**(3/n)*exp(q/(n*k)*(Eg0 - alpha*Theta**2/(Theta + beta))*(1/ThetaSrc - 1/Theta))*( exp(q*(V + Rs*I)/(Ns*n*k*Theta)) - 1 ) - (V + Rs*I)/Rp;
	F1 = I + (Rs*I - V)/Rp + Np*I0Src*(ThetaSrc/Theta)**(3/n)*exp(q/(n*k)*(Eg0 - alpha*Theta**2/(Theta + beta))*(1/ThetaSrc - 1/Theta))*exp(q*(V + Rs*I)/(Ns*n*k*Theta))*(q*(Rs*I - V)/(Ns*n*k*Theta));
	return [F0,F1]

# Function to set pretty minor and major ticks ----------------------

phi = (sqrt(5) + 1)/2

def prettyAxis(ax):
	ax.get_xaxis().set_minor_locator(mplot.ticker.AutoMinorLocator())
	ax.get_yaxis().set_minor_locator(mplot.ticker.AutoMinorLocator())
	# Seting axis lines x = 0 and y = 0 ----------------------------------
	ax.axhline(0,color = 'black', linewidth = 0.5)
	ax.axvline(0,color = 'black', linewidth = 0.5)
	# Setting major and minor grids --------------------------------------
	ax.grid(which = 'major', color = 'k', linestyle = '-', alpha = 0.5, linewidth = 0.5)
	ax.grid(which = 'minor', color = 'k', linestyle = ':', alpha = 0.5, linewidth = 0.5)
	# Setting axis and tick labels ---------------------------------------
	ax.tick_params(labelsize = 16)
	# Setting axis aspect ratio ------------------------------------------
	xleft, xright = ax.get_xlim();
	ybottom, ytop = ax.get_ylim();

	ax.set_aspect(abs(1/phi*(xright - xleft)/(ytop - ybottom)))
	return ax

# Initializing figures ----------------------------------------------
# I-V curve figure varying irradiance
fig1 = pyplot.figure();
ax1 = fig1.add_axes([0.1,0.1,0.8,0.8]);

# P-V curve figure varying irradiance
fig2 = pyplot.figure();
ax2 = fig2.add_axes([0.1,0.1,0.8,0.8]);

# I-V curve figure varying temperature
fig3 = pyplot.figure();
ax3 = fig3.add_axes([0.1,0.1,0.8,0.8]);

# P-V curve figure varying temperature
fig4 = pyplot.figure();
ax4 = fig4.add_axes([0.1,0.1,0.8,0.8]);

# I-V curve figure varying temperature and irradiance
fig5 = pyplot.figure();
ax5 = fig5.add_axes([0.1,0.1,0.8,0.8]);

# P-V curve figure varying temperature and irradiance
fig6 = pyplot.figure();
ax6 = fig6.add_axes([0.1,0.1,0.8,0.8]);

# Open-circuit voltage curve
fig7 = pyplot.figure();
ax7 = fig7.add_axes([0.1,0.1,0.8,0.8]);

# Efficiency curve
fig8 = pyplot.figure();
ax8 = fig8.add_axes([0.1,0.1,0.8,0.8]);

# -------------------------------------------------------------------
# MPP TRACKING VARYING IRRADIANCE
# -------------------------------------------------------------------
# Defining plot points and irradiance values ------------------------
size = 10**3;			# Number of plot points
g = linspace(10,1000,size,1);	# Array of irradiance values

# Preallocating solution vector
v = [0]*size;
i = [0]*size;
p = [0]*size;
guess = [25,4]		# MPP initial guess
for m in range(len(g)):
	
	# Defining funtion to be solved -----------------------------
	def f(z):
		G = g[m];
		Theta = ThetaSrc;
		return FMpp(z,G,Theta)

	sol = fsolve(f,guess)	# Solving
	v[m] = sol[0];
	i[m] = sol[1];
	p[m] = sol[0]*sol[1];
	guess = sol;	# Updating next current guess
		
ax1.plot(v,i, '--', linewidth = 1.5, color = 'black');
ax2.plot(v,p, '--', linewidth = 1.5, color = 'black');

# -------------------------------------------------------------------
# MPP TRACKING VARYING TEMPERATURE
# -------------------------------------------------------------------
# Defining plot points and temperature values -----------------------
size = 10**3 + 1;				# Number of plot points
t = linspace(-100 + 273.15,100 + 273.15,size,1);	# Array of irradiance values

# Preallocating solution vector
v = [0]*size;
i = [0]*size;
p = [0]*size;
guess = [25,4]		# MPP initial guess
for m in range(len(t)):

	def f(z): return FMpp(z,GSrc,t[m])	# Defining funtion to be solved

	sol = fsolve(f,guess)	# Solving
	v[m] = sol[0];
	i[m] = sol[1];
	p[m] = sol[0]*sol[1];
	guess = sol;	# Updating next current guess
		
ax3.plot(v,i, '--', linewidth = 1.5, color = 'black');
ax4.plot(v,p, '--', linewidth = 1.5, color = 'black');


# Defining plot points and irradiance values ------------------------
size = 100;			# Number of plot points
g = linspace(100,1000,10,1);	# Array of irradiance values

# -------------------------------------------------------------------
# IV-PV CURVE VARYING IRRADIATION
# -------------------------------------------------------------------
# Defining color map for plots --------------------------------------
cmap = mplot.cm.get_cmap('Spectral');
cindex = linspace(0,1,len(g),1);
colors = [0]*len(g);
for i in range(len(g)):
	colors[i] = cmap(cindex[i])

for m in range(len(g)):
	
	# Calculating open-circuit voltage for a given irradiance ---
	def f(V):
		I = 0;
		return FSim(V,I,g[m],ThetaSrc)

	Voc = fsolve(f,30);

	# Defining and preallocating vectors ------------------------
	v = linspace(0,Voc,size,1);	# Defining the voltage array for points
	i = [0]*size;			# Preallocating current solution vector
	p = [0]*size;			# Preallocating power solution vector
	guess = IrrSrc*g[m]/GSrc;	# Current initial guess

	# Solving equation through voltage list ---------------------
	for j in range(len(v)):
		def f(I): return FSim(v[j],I,g[m],ThetaSrc)	# Implicit function to solve

		i[j] = fsolve(f,guess)	# Solving
		p[j] = v[j]*i[j];	# Obtainig corresponding power
		guess = i[j];		# Updating next current guess
		
	ax1.plot(v,i, linewidth = 1.5, color = colors[m], label = format(g[m],'.0f'));
	ax2.plot(v,p, linewidth = 1.5, color = colors[m], label = format(g[m],'.0f'));

# -------------------------------------------------------------------
# IV-PV CURVE VARYING TEMPERATURE
# -------------------------------------------------------------------
# Defining plot points and temperature values -----------------------
size = 100;			# Number of plot points
t = linspace(-100 + 273.15,100 + 273.15,11,1);	# Array of irradiance values

# Defining color map for plots ---------------------------------------
cmap = mplot.cm.get_cmap('Spectral');
cindex = linspace(0,1,len(t),1);
colors = [0]*len(t);

for i in range(len(t)):
	colors[i] = cmap(cindex[i])

Voc = 30;		# Open-circuit voltage initial guess
for m in range(len(t)):	
	# Calculating open-circuit voltage for a given irradiance ----------
	def f(V): return FSim(V,0,GSrc,t[m])
	Voc = fsolve(f,Voc);

	# Defining and preallocating vectors -------------------------------	
	v = linspace(0,Voc,size,1);	# Defining the voltage array for points
	i = [0]*size;			# Preallocating current solution vector
	p = [0]*size;			# Preallocating power solution vector
	guess = IrrSrc;			# Current initial guess

	for j in range(len(v)): 
		def f(I): return FSim(v[j],I,GSrc,t[m])	# Implicit function to solve

		i[j] = fsolve(f,guess)	# Solving
		p[j] = v[j]*i[j];	# Obtainig corresponding power
		guess = i[j];		# Updating next current guess
		
	ax3.plot(v,i, linewidth = 1.5, color = colors[m], label = format(t[m] - 273.15,'.0f'));
	ax4.plot(v,p, linewidth = 1.5, color = colors[m], label = format(t[m] - 273.15,'.0f'));

# -------------------------------------------------------------------
# IV-PV CURVE VARYING IRRADIATION AND TEMPERATURE
# -------------------------------------------------------------------
t = linspace(-40 + 273.15,80 + 273.15,13,1);	# Array of temperature values
size = 10**3;
g = linspace(10,1000,size,1);			# Array of irradiance values

# Defining color map for plots --------------------------------------
cmap = mplot.cm.get_cmap('Spectral');
cindex = linspace(0,1,len(t),1);
colors = [0]*len(t);
for i in range(len(t)):
	colors[i] = cmap(cindex[i])

guess = [25,0];
for m in range(len(t)):

	# Defining and preallocating vectors -------------------------
	v = [0]*size;			# Defining the voltage array for points
	i = [0]*size;			# Preallocating current solution vector
	p = [0]*size;			# Preallocating power solution vectora
	
	for a in range(len(g)):
		def f(z): return FMpp(z,g[a],t[m]) # Defining function to be solved

		sol = fsolve(f,guess)	# Solving
		v[a] = sol[0];
		i[a] = sol[1];
		p[a] = sol[0]*sol[1];
		guess = sol;	# Updating next current guess

	ax5.plot(v,i, linewidth = 1.5, color = colors[m], label = format(t[m] - 273.15,'.0f'));
	ax6.plot(v,p, linewidth = 1.5, color = colors[m], label = format(t[m] - 273.15,'.0f'));

# --------------------------------------------------------------------
# Open-circuit voltage
# --------------------------------------------------------------------
# Defining plot points and temperature values ------------------------
size = 1000;			# Number of plot points
g = linspace(1,GSrc,size,1);			# Array of irradiance values
t = linspace(-40 + 273.15,80 + 273.15,13,1);	# Array of temperature values

# Defining color map for plots ---------------------------------------
cmap = mplot.cm.get_cmap('Spectral');
cindex = linspace(0,1,len(t),1);
colors = [0]*len(t);

for i in range(len(t)):
	colors[i] = cmap(cindex[i])

guess = 30;		# Open-circuit voltage initial guess

for m in range(len(t)):
	Voc = [0]*size;
	for a in range(len(g)):
		# Calculating open-circuit voltage for a given irradiance 
		def f(V): return FSim(V,0,g[a],t[m])
		Voc[a] = fsolve(f,guess);
		guess = Voc[a];
		
	ax7.plot(g,Voc, linewidth = 1.5, color = colors[m], label = format(t[m] - 273.15,'.0f'));

# -------------------------------------------------------------------
# EFFICIENCY
# -------------------------------------------------------------------
t = linspace(-40 + 273.15,80 + 273.15,13,1);	# Array of temperature values
size = 10**3;
g = linspace(1,1000,size,1);			# Array of irradiance values

# Defining color map for plots --------------------------------------
cmap = mplot.cm.get_cmap('Spectral');
cindex = linspace(0,1,len(t),1);
colors = [0]*len(t);
for i in range(len(t)):
	colors[i] = cmap(cindex[i])

guess = [25,0];
for m in range(len(t)):

	# Defining and preallocating vectors -------------------------
	v = [0]*size;			# Defining the voltage array for points
	i = [0]*size;			# Preallocating current solution vector
	p = [0]*size;			# Preallocating power solution vectora
	eta = [0]*size;
	for a in range(len(g)):
		def f(z): return FMpp(z,g[a],t[m]) # Defining function to be solved

		sol = fsolve(f,guess)	# Solving
		v[a] = sol[0];
		i[a] = sol[1];
		p[a] = sol[0]*sol[1];

		IL = g[a]/GSrc*( IrrSrc + mu*(t[m] - ThetaSrc))
		eta[a] = p[a]/(Np*IL*(v[a] + Rs*i[a]))
		guess = sol;	# Updating next current guess

	ax8.plot(g,eta, linewidth = 1.5, color = colors[m], label = format(t[m] - 273.15,'.0f'));
# --------------------------------------------------------------------
# PROCESSING FIGURES
# --------------------------------------------------------------------

pyplot.rc('text', usetex=True)
pyplot.rc('font', family='serif')

# Processing figure 1 ------------------------------------------------
ax1.set_xlabel('Voltage (V)', fontsize = 16)
ax1.set_ylabel('Current (A)', fontsize = 16)
ax1.set_title(r'Panel current versus voltage at $\theta = 25^{\circ}$C', fontsize = 16)
legend = ax1.legend(loc = 'upper right', shadow = True, title = 'Irradiation (W/m$^{-2}$)', fontsize = 12)
ax1 = prettyAxis(ax1)

# Processing figure 2 ------------------------------------------------
ax2.set_xlabel('Voltage (V)', fontsize = 16)
ax2.set_ylabel('Power (W)', fontsize = 16)
ax2.set_title(r'Panel power versus voltage at $\theta = 25^{\circ}$C', fontsize = 16)
legend2 = ax2.legend(loc = 'upper left', shadow = True, title = 'Irradiation (W/m$^{-2}$)', fontsize = 12)
ax2 = prettyAxis(ax2)

# Processing figure 3 ------------------------------------------------
ax3.set_xlabel('Voltage (V)', fontsize = 16)
ax3.set_ylabel('Current (A)', fontsize = 16)
ax3.set_title(r'Panel current versus voltage at $\phi = 1000$ W.m$^{-2}$', fontsize = 16)
legend3 = ax3.legend(loc = 'upper left', shadow = True, title = r'Temperature ($^{\circ}$C)', fontsize = 12)
ax3 = prettyAxis(ax3)

# Processing figure 4 ------------------------------------------------
ax4.set_xlabel('Voltage (V)', fontsize = 16)
ax4.set_ylabel('Power (W)', fontsize = 16)
ax4.set_title('Panel power versus voltage at $\phi = 1000$ W.m$^{-2}$', fontsize = 16)
legend4 = ax4.legend(loc = 'upper left', shadow = True, title = r'Temperature ($^{\circ}$C)', fontsize = 12)
ax4 = prettyAxis(ax4)

# Processing figure 5 ------------------------------------------------
ax5.set_xlabel('Voltage (V)', fontsize = 16)
ax5.set_ylabel('Current (A)', fontsize = 16)
ax5.set_title('Panel MPP curves', fontsize = 16)
legend5 = ax5.legend(loc = 'upper left', shadow = True, title = r'Temperature ($^{\circ}$C)', fontsize = 12)
ax5 = prettyAxis(ax5)

# Processing figure 6 ------------------------------------------------
ax6.set_xlabel('Voltage (V)', fontsize = 16)
ax6.set_ylabel('Power (W)', fontsize = 16)
ax6.set_title('Panel MPP curves', fontsize = 16)
legend6 = ax6.legend(loc = 'upper left', shadow = True, title = r'Temperature ($^{\circ}$C)', fontsize = 12)
ax6 = prettyAxis(ax6)

# Processing figure 7 ------------------------------------------------
ax7.set_xlabel(r'Irradiance $\phi$ (W.m$^{-2}$)', fontsize = 16)
ax7.set_ylabel(r'Open-circuit voltage $V_{OC}$ (V)', fontsize = 16)
ax7.set_title('Panel open-circuit voltage versus irradiance, parametrized by temperature', fontsize = 16)
legend7 = ax7.legend(loc = 'lower right', shadow = True, title = r'Temperature ($^{\circ}$C)', fontsize = 12)
ax7 = prettyAxis(ax7)

# Processing figure 8 ------------------------------------------------
ax8.set_xlabel(r'Irradiance $\phi$ (W.m$^{-2}$)', fontsize = 16)
ax8.set_ylabel(r'Panel efficiency $\eta$', fontsize = 16)
ax8.set_title('Panel efficiency versus irradiance, parametrized by temperature', fontsize = 16)
legend8 = ax8.legend(loc = 'lower right', shadow = True, title = r'Temperature ($^{\circ}$C)', fontsize = 12)
ax8 = prettyAxis(ax8)

# --------------------------------------------------------------------
# SHOW FIGURES 
# --------------------------------------------------------------------
pyplot.show();
