function [] = findLabel(I,tv)
    H = rgb2hsv_fast(I,'','H');
    S = rgb2hsv_fast(I,'','S');
    
    M = (H < tv(1) | H > tv(2)) & S > tv(3);
    
end