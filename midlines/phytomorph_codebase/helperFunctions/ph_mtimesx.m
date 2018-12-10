function [] = ph_mtimesx ()
    try
        C = mtimesx(M,E);
        C = M*E;
    end
    catch
    end
end