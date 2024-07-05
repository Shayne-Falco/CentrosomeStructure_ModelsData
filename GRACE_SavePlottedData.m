clear all
fig = openfig('n1n2.fig');
fig.Visible = true;
a = get(gca,'Children');
xdata = get(a, 'XData');
ydata = get(a, 'YData');


data = [xdata' ydata']

save('my_data.out', 'data', '-ascii');