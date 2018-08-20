function plotLinearAbsBAF(Y, CID, chromIDs,sample,colorSpec)
%plotLinear(X,Y, CID, chromIDs, sample)
% Plot the entile genome
n = numel(chromIDs);
chr_endIdx = zeros(1, n);
chr_data_len = zeros(1,n);
for i = 1:n
    tmp = CID == chromIDs(i);
    chr_endIdx(i) = find(tmp, 1, 'last');
    chr_data_len(i) = length(find(tmp));
end

x_lims = [0 chr_endIdx(n)];
y_lims = [0, 0.5];
%y_lims = [min(min(Y), 0), ceil(max(max(Y), 1))];

% Draw a vertical bar at the end of a chromosome to indicate the border
x_vbar = repmat(chr_endIdx, 2, 1);
y_vbar = repmat(y_lims', 1, n);
offset = double(y_lims(1))-0.05;
% Label the autosome with their chromosome numbers
x_label = chr_endIdx - ceil(chr_data_len/2);
% % y_label = repmat(offset, 1, length(x_label));
y_label = zeros(1, length(x_label))+ offset;
chr_labels = cellstr(num2str(chromIDs'));
chr_labels = strrep(chr_labels, '23', 'X');
%chr_labels = strrep(chr_labels, '24', 'Y');
% A gray zero line
%x_zero = [0; length(Y)];
%y_zero = [0;0.5];

absY = abs(0.5-Y);

%== Plot
hold on
h_ratio = plot(absY, '.','color',colorSpec);

line(x_vbar, y_vbar, 'color', [0.8 0.8 0.8]);
%line(x_zero, y_zero, 'Color', [1 0.3 0.3], 'Linewidth', 1.5);
text(x_label, y_label, chr_labels,...
    'Fontsize', 10, 'HorizontalAlignment', 'Center');
h_axis = get(h_ratio, 'parent');
set(h_axis, 'xtick', [], 'ygrid', 'on', 'box', 'on',...
            'xlim', x_lims, 'ylim', y_lims)

title(sample, 'Interpreter', 'none','fontsize',14)
xlabel({'', 'Chromosome'})
ylabel('Absolute BAF','fontsize',16)
hold off
end % end of function