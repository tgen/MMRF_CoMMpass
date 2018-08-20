function x = readTSV(tsvFile)


fid = fopen(tsvFile);
x = textscan(fid,'%f\t%f\t%f','headerlines',1);
fclose(fid);

x= cell2mat(x);