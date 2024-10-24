p0byrng = p0byrngT2ss(:, 11:17) ;

SM = zeros(size(p0byrng,1),1);
SMvar = zeros(size(p0byrng,1),1);
p0 = 3.75:0.01:3.90;
p0 = p0';
%now, we want to find averages for each row
for i = 1:size(p0,1)
    row = p0byrng(i,:);
    rowdata = row(find(row));
    if (length(rowdata) > 4)
        SM(i) = mean(rowdata);
        SMvar(i) = std(rowdata) / sqrt(length(rowdata)) ;
    else
        SM(i) = NaN;
        SMvar(i) = NaN;
    end
end