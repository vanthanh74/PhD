#!/usr/local/bin/python
import numpy as np

class ODE_solver:
    """
        Brusselator system:
            dx/dt = A + x^2 y - (B+1) x
            dy/dt = B x -x^2 y
        with Explicit Euler or explicit RK4
    """

    def __init__(self,
                init_state = [0., 1.],
                A = 1.,
                B = 3.,
                dt = 0.1,
                num_scheme="EE"):

        self.init_state = np.asarray(init_state, dtype=float)
        self.state = self.init_state.copy()
        self.A = float(A)
        self.B = float(B)
        self.dt = dt
        self.time = 0.
        self.num_scheme=num_scheme

    def f(self,s):
        x = s[0]
        y = s[1]
        return np.asarray([self.A+(x**2)*y-(self.B+1.)*x, self.B*x - (x**2)*y])


    def step(self):
        if (self.num_scheme == "EE"):
            # Update initial state
            self.init_state = self.state.copy()
            # Update time
            self.time += self.dt
            # Update state
            self.state += self.dt * self.f(self.init_state)
        else:
            if(self.num_scheme=="RK4"):
                # Update initial state
                self.init_state = self.state.copy()
                # Update time
                self.time += self.dt
                # Compute RK coefficients
                K1 = self.f(self.init_state)
                K2 = self.f(self.init_state + (self.dt / 2) * K1)
                K3 = self.f(self.init_state + (self.dt / 2) * K2)
                K4 = self.f(self.init_state + self.dt * K3)
                # Update state
                self.state += self.dt*K1/6+(self.dt/3)*(K2+K3)+self.dt*K4/6
            else:
                print("Error in N_Body.step: scheme not implemented")

    def propagate(self, state0, t0, nstep):
        self.init_state = np.asarray(state0, dtype=float)
        self.state = self.init_state.copy()
        self.time = t0
        for n in range(nstep):
            self.step()
        return self.state.copy()