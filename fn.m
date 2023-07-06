classdef fn
    methods (Static)
        function output = cr3bp(t,state,u)
            rho1 = sqrt((state(1)+u)^2 + state(2)^2 + state(3)^2);
            rho2 = sqrt((state(1)-(1-u))^2 + state(2)^2 + state(3)^2);

            ax = -(1-u)*(state(1)+u)/rho1^3 - u*(state(1)-(1-u))/rho2^3 + state(1) + 2*state(5);
            ay = -(1-u)*state(2)/rho1^3 - u*state(2)/rho2^3 + state(2) - 2*state(4);
            az = -(1-u)*state(3)/rho1^3 - u*state(3)/rho2^3 ;

            output  = [state(4);state(5);state(6);ax;ay;az];
        end
        function C = jacobiconst(state,u)
            rho1 = sqrt((state(1)+u)^2 + state(2)^2 + state(3)^2);
            rho2 = sqrt((state(1)-(1-u))^2 + state(2)^2 + state(3)^2);
            v = [state(4),state(5),state(6)];
            C = (state(1)^2 + state(2)^2 + 2*(1-u)/rho1 + 2*u/rho2) - norm(v)^2;
        end
        function output = collinear(x)
            u = 0.0122;
            output = x -(((1-u)*(x+u)/(abs(x+u))^3)) - (u*(x-1+u)/(abs(x-1+u))^3);
        end
    



    end
end
