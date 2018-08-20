function baf = readBAF(bafTXT)

fid = fopen(bafTXT);
baf = textscan(fid,'%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f','headerlines',1);
fclose(fid);
baf =cell2mat(baf);