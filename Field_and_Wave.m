% Seyed Mohammad Kashani 9423422
% field and wave Dr. Azadi
% multilayer slub transmission problem 

mu_zero = 4*pi*10^(-7);
epsilon_zoro = 8.85*10^(-12);

% Enter TE OR TM
polarization = "TE";
wave_frequency = 2000000 ;

% Enter layers data
mu_r = mu_zero.*[1 1 1 1 1];
epsilon_r = epsilon_zoro.*[1 1 1 1 1];
layer_length = 100.*[ 2 3 2];
% Enter incident angle in radian
incident_angle_radian = .5;

if (length(mu_r) == length(epsilon_r)) && (length(layer_length) == length(mu_r)-2)
    for i = 1:length(mu_r)
       eta(i) = sqrt(mu_r(i) ./ epsilon_r(i));
    end
    
    for l = 1:length(mu_r)
        beta(l) = 2*pi*wave_frequency*sqrt(mu_r(l)*epsilon_r(l));
    end
    
    teta_i(1) = incident_angle_radian;
    for m = 1:length(mu_r)-1
        teta_t(m) = asin(beta(m)*sin(teta_i(m))/beta(m+1));
        teta_i(m+1) = teta_t(m);
    end
    
    for s = 2:length(mu_r)-1
%         beta_m(s-1) = beta(s)*sqrt(1-((mu(1)*epsilon(1))/(epsilon_zoro*mu_zero)) * (sin( incident_angle_radian))^2) ; 
          beta_m(s-1) = beta(s)*cos(teta_t(s));
    end
    if polarization == "TE"    
        for b = 1: length(mu_r)   
            z(b) = eta(b)/cos(teta_i(b));    
        end     
    elseif polarization == "TM"
        for p = 1: length(mu_r)   
            z(p) = eta(p).*cos(teta_i(p));    
        end
    else
        disp("choose TE or TM");
    end

    ABCD = [1 0;
            0 1];
    
    for o = 2: length(mu_r)-1   
        A_D = cos(beta_m(o-1)*layer_length(o-1));
        B = 1i*z(o)*sin(beta_m(o-1)*layer_length(o-1));
        C = (1i*sin(beta_m(o-1)*layer_length(o-1)))./ z(o);
        
        A (1,1) = A_D;
        A (1,2) = B;
        A (2,1) = C;
        A (2,2) = A_D;
        %A;
        ABCD = ABCD * A ;
    end
    
    ABCD
    R = (ABCD(1,1) + (ABCD(1,2)/z(length(mu_r))) - z(1)*(ABCD(2,1)+((ABCD(2,2))/(z(length(mu_r))))))/(ABCD(1,1) + (ABCD(1,2)/z(length(mu_r))) + z(1)*(ABCD(2,1)+((ABCD(2,2))/(z(length(mu_r))))))
    T = 2 /(ABCD(1,1) + (ABCD(1,2)/z(length(mu_r))) + z(1)*(ABCD(2,1)+((ABCD(2,2))/(z(length(mu_r))))))

    
else
    disp("Carefully enter the number of layers and the permeability and the permittivity");
end


