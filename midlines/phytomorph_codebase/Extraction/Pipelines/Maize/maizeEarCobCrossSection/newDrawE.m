function [X,Y,RAD,Xn,Yn,RADn] = newDrawE(major,minor,T,ROT,CEN,NS)

    stepSIZE = minor/NS;

    TH = linspace(-pi,pi,T);
    
    
    SZ = round([T NS]);
   
    [n1 n2] = ndgrid(linspace(-pi,pi,T),linspace(0,1,NS));
    
    
    RAD = [linspace(major-minor,major,NS)' linspace(1,minor,NS)'];
    RADn = [major*n2(1,:)'  minor*n2(1,:)'];
    
    Xn = major*n2.*cos(n1);
    Yn = minor*n2.*sin(n1);
    GR = [ROT*[Xn(:) Yn(:)]']';
    Xn = GR(:,1) + CEN(1);
    Yn = GR(:,2) + CEN(2);
    Xn = reshape(Xn,round(SZ));
    Yn = reshape(Yn,round(SZ));

  
    
    for e = 1:NS
       X(:,e) = major*cos(TH)';
       Y(:,e) = minor*sin(TH)';
       major = major - stepSIZE;
       minor = minor - stepSIZE;
    end
    X = flip(X,2);
    Y = flip(Y,2);
    GR = [ROT*[X(:) Y(:)]']';
    X = GR(:,1) + CEN(1);
    Y = GR(:,2) + CEN(2);
    X = reshape(X,round(SZ));
    Y = reshape(Y,round(SZ));
    
    
    
end