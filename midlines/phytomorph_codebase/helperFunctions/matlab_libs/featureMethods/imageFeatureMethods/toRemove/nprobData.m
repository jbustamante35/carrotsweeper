function [p] = nprobData(data,model)
    dm = size(data,2);    
    p = probData(data,model);
    mx = probData(model(1:dm),model);
    p = p/mx;
end