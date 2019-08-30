from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Rosario']

import matplotlib.pyplot as plt

import pylab

import numpy as np

# Pie chart
labels = ['Coal', 'Gas', 'Other thermal', 'Nuclear', 'Conventional Hydroelectric', 'Pumped storage', 'Wind', 'Solar', 'Biomass']
sizes = [3946, 188,89,213, 1144, 31, 241, 67, 65]
#colors
#colors = ['#ff9999','#66b3ff','#99ff99','#ffcc99']

#cm = pylab.get_cmap('Accent')
#colors = [cm(1.*i/11) for i in range(11)]
#print(colors)

fig1, ax1 = plt.subplots()
patches, texts, autotexts = ax1.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90, wedgeprops=dict(width=0.5, edgecolor='w'))
for text in texts:
    text.set_color('k')
for autotext in autotexts:
    autotext.set_color('k')


# Equal aspect ratio ensures that pie is drawn as a circle
ax1.axis('equal')  
plt.tight_layout()
plt.show()
