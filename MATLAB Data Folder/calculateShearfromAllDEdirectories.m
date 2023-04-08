
dirs = ls('*dE');
for i=1:size(dirs, 1)
    dir_label = dirs(i, :);
    namesplit = split(dir_label, "_");
    if size(namesplit, 1) < 3
        prefix = namesplit{1};
    else
        prefix = strcat(namesplit{1},'_', namesplit{2});
    end

    if prefix(2:2) == 'B'
        if prefix(1:1) == 'B'
            Lx = 500;
            Ly = 500;
        elseif prefix(1:1) == 'L'
            Lx = 250;
            Ly = 250;
        end
    else
        Lx = 353.553391;
        Ly = 353.553391; 
    end

    if (strcmp(prefix,'L2Qu') || strcmp(prefix,'M2Qu'))
        p0num = 31;
    else
        p0num = 16;
        rngnum = 20;
        ShearModulusCalculator(prefix, Lx, Ly, p0num, rngnum);
    end
    close(gcf);
end