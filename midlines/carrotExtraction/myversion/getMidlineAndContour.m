function [mline, cntr] = getMidlineAndContour(msk, vis)

try
    out = isolate_carrot_Roots(msk, vis, [], []);
%     out = isolateRoots(msk, vis, [], []);
    mline = out(1).midlines.data';
    
    rm = mline(:,1) < 300;
    mline(rm,:) = [];
    mline(:,1) = mline(:,1) - 300;
    
    cntr = out(1).contours.data';
    rm = cntr(:,1) < 300;
    cntr(rm,:) = [];
    cntr(:,1) = cntr(:,1) - 300;
    
catch e
    fprintf(2, 'Error extracting Midline and Contour\n%s\n', e.getReport);
    mline = [];
    cntr = [];
end

end
