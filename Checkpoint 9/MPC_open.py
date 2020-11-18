# Nonlinear MPC with Casadi: Open-loop system
from casadi import *

T = 10.0;
N = 20;

# Declare model variables
x1 = MX.sym('x1')
x2 = MX.sym('x2')
x = vertcat(x1, x2)
x0 = MX.sym('x0', 2)
u = MX.sym('u')

# Model equations
xdot = vertcat(x2, (1 - x1**2)*x2 - x1 + u)

# Objective term
L = x1**2 + x2**2 + u**2

# Formulate discrete dynamics
dae = {'x': x, 'p': u, 'ode': xdot, 'quad': L}
opts = {'tf': T/N}
F = integrator('F', 'cvodes', dae, opts)

# Define the optimization problem
w, w0, lbw, ubw = [], [], [], []
g, lbg, ubg = [], [], []
J = 0

# Formulate the NLP
Xk = x0
for k in range(N):
    # Constraint on control
    Uk = MX.sym('U_' + str(k))
    w += [Uk]
    lbw += [-1]
    ubw += [+1]
    w0 += [0]
    # Move one sampling period
    Fk = F(x0 = Xk, p = Uk)
    Xk = Fk['xf']
    # Update the cost function
    J = J + Fk['qf']
    # Constraints on the state
    g += [Xk[1]]
    lbg += [-0.25]
    ubg += [inf]
    
prob = {'f': J, 'x': vertcat(*w), 'g': vertcat(*g), 'p': x0}

solver = nlpsol('solver', 'ipopt', prob, {'ipopt': {'max_iter': 200}})

x0 = [1, 0]
sol = solver(x0 = w0, lbx = lbw,
             ubx = ubw, lbg = lbg,
             ubg = ubg, p = x0)
w_opt = sol['x']
# Plot the solution
# Extract the optimal trajectory
u_opt = w_opt
x_opt = [x0];
for k in range(N):
    Fk = F(x0 = x_opt[-1], p = u_opt[k])
    x_opt += [Fk['xf'].full()]
x1_opt = [r[0] for r in x_opt]
x2_opt = [r[1] for r in x_opt]
# Plot the trajectories
tgrid = [T/N*k for k in range(N + 1)]
import matplotlib.pyplot as plt
plot = plt.figure()
plt.clf()
plt.plot(tgrid, x1_opt, '--')
plt.plot(tgrid, x2_opt, '-')
plt.grid(True)
plt.step(tgrid, vertcat(DM.nan(1), u_opt), '-.')
plt.xlabel('t')
plt.legend(['x1', 'x2', 'u'])
plt.title('Open-loop MPC results')
plt.show()    
# Save figure in .eps for best quality
plot.savefig('MPC_open.eps', format='eps') 