function A = extractData(folderName, sizeA)

lst = dir(folderName);
names = {lst.name};

A = zeros(sizeA,121,80);
i = 1;
for name = names
    if contains(name, 'subject')
        fileID = fopen(string(folderName) + '/' + string(name));
        tmp = fread(fileID,[121 sizeA],'double');
        A(:,:,i) = tmp';
        i = i+1;       
    end

end

end