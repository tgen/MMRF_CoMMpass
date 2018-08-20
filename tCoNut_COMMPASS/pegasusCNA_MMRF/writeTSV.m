function writeTSV(fileName,cnaMat,col3Header)
%writeTSV(fileName,cnaMat)
%
% Write TSV files for Copy Number Analysis (CNA)
%
% Jessica Aldrich
% Tgen 
% January 9, 2014

fid = fopen(fileName,'w+');

fprintf(fid,'%s\t%s\t%s\n','Chr','Position',col3Header);
for i=1:size(cnaMat,1)
    fprintf(fid,'%d\t%d\t%f\n',cnaMat(i,:));
end

fclose(fid);