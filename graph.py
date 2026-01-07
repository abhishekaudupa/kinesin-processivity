import params
import numpy as np
import matplotlib.pyplot as plt

def binding_potential(chemical_state):
    limit = 20
    X = np.linspace(-limit, limit, 1000)
    ri = np.array([params.dx*i for i in range(-10*limit, 10*limit)])

    Vb = []
    match(chemical_state):
        case 'T':
            Vb0 = 16 * params.V0
        case 'D':
            Vb0 = 7 * params.V0
        case 'E':
            Vb0 = 10 * params.V0
    
    for x in X:
        Vb.append(Vb0 * np.sum(-np.exp(-(x-ri)**2)) - params.s*params.V0*x)

    return (X, np.array(Vb) * 1e6)

def binding_force(chemical_state):
    limit = 20
    X = np.linspace(-limit, limit, 1000)
    ri = np.array([params.dx*i for i in range(-10*limit, 10*limit)])

    Fb = []
    match(chemical_state):
        case 'T':
            Vb0 = 16 * params.V0
        case 'D':
            Vb0 = 7 * params.V0
        case 'E':
            Vb0 = 10 * params.V0
    
    for x in X:
        Fb.append(-2*Vb0 * np.sum(np.exp(-(x-ri)**2) * (x-ri)) + params.s * params.V0)

    return (X, np.array(Fb) * 1e6)

chemical_states = ['E', 'T', 'D']
description = ['Empty', 'ATP Bound', 'ADP Bound']
colors = ['blue', 'red', 'green']

fig, axes = plt.subplots(2, 1, figsize=(8, 4.5))
ax = axes.flatten()
for (chemical_state, desc, color) in zip(chemical_states, description, colors):
    X, Vb = binding_potential(chemical_state)
    X, Fb = binding_force(chemical_state)
    ax[0].plot(X, Vb, label = desc, color = color, lw=0.5)
    ax[0].set_xlabel('Microtubule axis, nm')
    ax[0].set_ylabel(r'Binding Potential, pN$\cdot$nm')
    ax[0].grid()
    ax[0].legend()
    
    ax[1].plot(X, Fb, label = desc, color = color, lw=0.5)
    ax[1].set_xlabel('Microtubule axis, nm')
    ax[1].set_ylabel(r'Binding Force, pN')
    ax[1].grid()
    ax[1].legend()

ax[0].set_xticks(range(-16, 17, 8))
ax[1].set_xticks(range(-16, 17, 8))
plt.tight_layout()
plt.savefig('potential.png', dpi=300, bbox_inches='tight')
