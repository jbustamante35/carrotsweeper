function [G] = issueGradesOverData(result,gidxI,icE,icU,iPDF_u,iPDF_c,gidxE,ecE,ecU,ePDF_u,ePDF_c,rEI,rIE,idU,idE,myiPDF,myePDF)
    for e = 1:numel(result)
        fprintf(['***************************************************************\n']);
        fprintf(['* Start Grade \n'])
        fprintf(['***************************************************************\n']);
        G{e} = issueGrades(result{e}(:,1:end-4),gidxI,icE,icU,iPDF_u,iPDF_c,gidxE,ecE,ecU,ePDF_u,ePDF_c,rEI,rIE,idU,idE,myiPDF,myePDF);
        fprintf(['***************************************************************\n']);
        fprintf(['* End Grade \n'])
        fprintf(['***************************************************************\n']);
    end
end