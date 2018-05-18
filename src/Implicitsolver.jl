using DifferentialEquations
using Plots
# pyplot()
function s(f,KK,MM,al,be,n,u0,t)
    C=al*MM+be*KK
    A=[C*Mi KK*Mi;eye(n) zeros(n,n)]
    b=[Mi*f.(t)...;zeros(n)...]
    g(u,p,t)=A*u+b
    prob=ODEProblem(g,u0,tspan)
    sol=solve(prob,ImplicitRKMil())
    h=plot(sol)
    display(h)
end
M=400; m=35
K=50000; k=200000
H=300
function f(t)
    [H*sin.(t);0]
end
MM=M*0.5*eye(2)+0.5*m*eye(2)
KK=[K+k -K;-K  K]
Mi=inv(MM)
# t=[0:.1:10.0;]
tspan=(0.0,2.0)
u0=zeros(4)
# f(t)=[H*sin(t);0]
al=0.5; be=0.02; n=2
s(f,KK,MM,al,be,n,u0,tspan)
