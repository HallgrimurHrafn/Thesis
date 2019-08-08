function [F,G]=StateMatrix(k_a,f_L,xi,permRatio,s_x,ReducedKlippel)
    ModelDetails;
    if ReducedKlippel
%         Reduced Klippel model
%         mechanical
        k_phi=@(x,i)2*u0*A*Fm^2/Deff^3*f_L(x,i);
        T=@(x,i) 2*u0*A*N*Fm/Deff^2 * f_L(x,i);
%         electrical
        Ledx=@(x,i) -2*u0*A*N^2/Deff*f_L(x,i)^2*xi(x,i)*2*s_x*i/Deff;
        Ledi=@(x,i) 2*u0*A*N^2/Deff*f_L(x,i)*(1-2*i*f_L(x,i)*xi(x,i));
%         T=@(x,i) 2*u0*A*N*Fm/Deff^2*f_L(x,i);
%         state matrices IGNORING MECHANICAL LOSS FACTOR
        F=@(x,i) [-R_e/Ledi(x,i)       0      -(T(x,i)+Ledx(x,i))/Ledi(x,i)
                       0               0                  1
                    T(x,i)/M    (k_phi(x,i)-k_a(x))/M     -0.146/M];
        G=@(x,i) [1/Ledi(x,i);0;0];
    else
%         Full klippel model

%         Mechanical System
%         mechanical transduction distortion
        T_med=@(x,i) 2*u0*A*N^2*i*x/Deff^3*f_L(x,i)^2;
%         mechanical transduction
        T_me=@(x,i) 2*u0*A*Fm*N/Deff^2*f_L(x,i)^2*(1+x^2/Deff^2+permRatio(x,i));
%         mechanical magnetic stiffness compensation
        k_phi=@(x,i)2*u0*A*Fm^2/Deff^3*f_L(x,i)^2*(1+permRatio(x,i));

        
%         Electrical System
%         Inductance derivative term with respect to current
        Ledi=@(x,i) 2*u0*A*N^2/Deff * f_L(x,i)*(1-2*i*f_L(x,i)*xi(x,i));
%         Inductance derivative term with respect to displacement
        Ledx=@(x,i) 4*u0*A*N^2/Deff^3 * f_L(x,i)^2*(x-D*s_x*xi(x,i));
%         Transduction derivative term with respect to current
        Tdi=@(x,i) -4*u0*A*Fm*N*x/Deff^2 * f_L(x,i)^2 *xi(x,i);
%         Transduction derivative term with respect to displacement
        Tdx=@(x,i) 2*u0*A*Fm*N/Deff^2 * f_L(x,i) *(1+f_L(x,i)*(2*x^2/Deff^2 ...
                            -2*x*s_x/Deff*xi(x,i)));
%         State Matrices
        F=@(x,i) [-R_e/(Ledi(x,i)+Tdi(x,i)),  0,    -(Tdx(x,i)+Ledx(x,i))/(Ledi(x,i)+Tdi(x,i))
                           0                  0                       1
                     (T_me(x,i)+T_med(x,i))/M,   (k_phi(x,i)-k_a(x))/M,          -0.146/M];        
        G=@(x,i) [1/(Ledi(x,i)+Tdi(x,i));0;0];   
    end
    
end