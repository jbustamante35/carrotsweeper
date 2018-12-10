function [table] = addPointToTable(table,type,index,pointName,point)
    switch type
        case 'point'
            table{index,[pointName '_1']} = point(1);
            table{index,[pointName '_2']} = point(2);
        case 'boundingBox'
            for e = 1:4
                table{index,[pointName '_1']} = point(1);
                table{index,[pointName '_2']} = point(2);
            end
    end
    
end