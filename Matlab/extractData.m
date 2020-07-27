function A = extractData(folderName, sizeA)

lst = dir(folderName);
names = {lst.name};

A = zeros(sizeA,156,61);
i = 1;
for name = names
    if contains(name, 'subject')
        fileID = fopen(string(folderName) + '/' + string(name));
        tmp = fread(fileID,[156 sizeA],'double');
        A(:,:,i) = tmp';
        i = i+1;       
    end

end

end