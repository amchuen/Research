function int = simpInteg(xdata, ydata)


dx = xdata(2) - xdata(1);
int = ydata(1) + ydata(end) + dx/3 .*(4.*sum(ydata(2:2:end)) + 2.*sum(ydata(3:2:end-1)));

end