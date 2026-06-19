
## Problem definition

Recall that the optimal control problem we are considering is:

```math
    \left\{ \begin{array}{ll}
    \displaystyle \min_{x,u} J(x,y,u) \coloneqq \frac{1}{2}\int_{0}^{T} \big(x^2(t) + y^2(t)\big) ~\mathrm dt \\[1em]
    \text{s.c.}~\dot x(t) = u(t), & t\in [t_0, t_f]~\mathrm{a.e.}, \\[0.5em]
    \displaystyle \phantom{\mathrm{s.c.}~}\dot y(t) = - u(t) - \frac{y(t)}{2}, & t\in [t_0, t_f]~\mathrm{a.e.}, \\[0.5em]
    \phantom{\mathrm{s.c.}~} -M \leq u(t) \leq M, & t\in [t_0, t_f], \\[0.5em]
    \phantom{\mathrm{s.c.}~} \big(x(0), y(0)\big) = (1, -0.5), \quad \big(x(T), y(T)\big) = (0.5, 0.5),
    \end{array} \right.
```
with ``M>0`` and ``T>0``. 

This optimal control problem is associated to the stationnary optimization problem 

```math
    \min{(x,y,u)} \big \{  x^2 + y^2 \mid (x,y,u) \in \mathbb R \times \mathbb R \times [-M, M], u = u - \frac{y}{2} = 0 \big\}.
```
The static solution is thus given by ``(x^\star, y^\star, u^\star) = (0, 0, 0)``. This solution corresponds to the static equilibrium ``(x,y,u)`` which minimizes the cost ``J(x,y,u)`` under the control constraint ``|u|\leq M``. Our goal is to numerically show that this problem has the `turnpike` property, i.e. if the final time ``T`` is large enough, the optimal trajectory has the following form : starting from the initial state ``(x(0), y(0)) = (1, -0.5)``, the trajectory try to reach as fast as possible the static equilibrium ``(x^\star, y^\star, u^\star) = (0, 0, 0)`` and to stay close to it until it has to leave it to reach the final state ``(x(T), y(T)) = (1/2, 1/2)``.


!!! info "Main goals"
    Solve the optimal control problem using the indirect multiple shooting method, and highlight the turnpike property by using a continuation with respect to the final time ``T``.

To do this, let us start by importing the needed package for this page, and defining the optimal control problem. 

```@example main 
# Import required packages for optimal control analysis
using Plots                                 # For plotting and visualization
using OptimalControl                        # Main package for optimal control problems
using NLPModelsIpopt                        # Interface for Ipopt nonlinear solver
using OrdinaryDiffEq                        # For solving differential equations
using MINPACK                               # For nonlinear equation solving (shooting method)
using DifferentiationInterface              # For automatic differentiation
using ForwardDiff                           # Forward mode automatic differentiation
using Ipopt, Optimization, OptimizationMOI  # Additional optimization tools

M = 1.                   
s0 = [0., 1, 0]            
xT = 0.5; yT = 0

F0(s) = [
    0.5*(s[2]^2 + s[3]^2);
    0;
    -s[3]/2
]

F1(s) = [
    0;
    1;
    -1
]

# Define the optimal control problem using OptimalControl.jl DSL
ocp(T) = @def begin
    t ∈ [0, T], time                                # Time
    s = (c, x, y) ∈ R³, state                       # State
    u ∈ R, control                                  # Control
    -M ≤ u(t) ≤ M                                   # Control's constraint
    s(0) == s0                                      # Initial state
    x(T) == xT                                      # Final state
    y(T) == yT                                      # Final state
    ṡ(t) == F0(s(t)) + u(t)*F1(s(t))                # Dynamics
    c(T) → min                                      # Cost
end

nothing; # hide
```

!!! info "Note"
    The optimal control problem parametrized with the time horizon `T` that can be adjusted in oreder to study the turnpike property.

## Indirect framework

Denoting ``s = (c, x, y)`` the augmented state, and ``p = (p_c, p_x, p_y)`` the costate, the pseudo-Hamiltonian is given by 

```math
H(s, p, u) = p_c \frac{1}{2}(x^2 + y^2) + p_x u + p_y (u - \frac{y}{2}).
```

Since the state dynamic is linear with respect to the control 

```math 
\dot s(t) = F_0(s(t)) + u(t) F_1(s(t)) \quad \text{for all } t \in [0, T]
```

where 

```math
F_0(s) = \begin{pmatrix} \frac{1}{2}(x^2 + y^2) \\ 0 \\ -\frac{y}{2} \end{pmatrix}, \quad F_1(s) = \begin{pmatrix} 0 \\ 1 \\ -1 \end{pmatrix},
```

the pseudo-Hamiltonian is linear with respect to the control and can be written as 

```math
H(s, p, u) = H_0(s, p) + u H_1(s, p) 
```

where 

```math
H_0(s, p) = p \cdot F_0(s) \quad \text{and} \quad H_1(s, p) = p \cdot F_1(s).
```

Thanks to the differential geometry tools provided by 'OptimalControl.jl' decribed in the [documentation](https://control-toolbox.org/OptimalControl.jl/stable/manual-differential-geometry.html), the two Hamitonians ``H_0`` and ``H_1`` can be directly defined as the lift of the vector fields `F0` and `F1`. 

```@example main
# Lift into (x,λ) space
H0 = Lift(F0)
H1 = Lift(F1)

# Pseudo-Hamiltonian
H(s,p,u) = H0(s,p) + u*H1(s,p)
```
