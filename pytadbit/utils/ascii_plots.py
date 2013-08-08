"""
06 Aug 2013


"""


from math import log, exp


def arange(beg, end, step):
    return [beg + i * step for i in xrange(int(abs(beg-end)/step+.5))]


def ascii_plot (ydata, xdata=None, logscale=False, pch='o', title='plot',
                xlabel='X', ylabel='Y', width=72, height=50):
    """
    Curve (ASCII format).

    :param ydata: list of values to be plotted
    :param None xdata: x coordinate corresponding to ydata. If None will range
       between 1 and the length of ydata.
    :param False logscale: display data with logarithmic Y axis
    :param 'o' pch: string for points (whatever + = - * etc...)
    :param 'plot' title: string for title of the plot
    :param 'X' xlabel: label for the X axis
    :param 'Y' ylabel: label for the Y axis
    :param 100 width: width in term of characters
    :param 100 height: height in term of characters

    :returns: string corresponding to plot

    **Example:**
    
    print ascii_plot([0,5,9,18,7], width=60, height=10)

 plot
 ----

Y
18.0000+
       |                                             o               
       |                                                             
       |                                                             
       |                                                             
       |                                                             
9.0000 +                              o                              
       |                                                            o
       |               o                                             
       |                                                             
       |                                                             
0.0000 +o                                                            
     0 +---------+---------+---------+---------+---------+---------+
       1.000000  1.666667  2.333333  3.000000  3.666667  4.333333  5.000000  

                                    X                              

    """
    if not xdata:
        xdata = range(1, len(ydata)+1)
    yydata = []
    logf = log if logscale else lambda x: x
    expf = exp if logscale else lambda x: x
    for i in ydata:
        try:
            yydata.append(logf(i))
        except ValueError:
            yydata.append(float('-inf'))
    ydiff = float(abs(float(min(yydata)) - max(yydata))/(height * 2))
    y_arange = [(i - ydiff, i + ydiff) for i in
                sorted(arange(min(yydata), max(yydata) + ydiff, ydiff * 2), reverse=True)]
    xdiff = float(abs(float(min(xdata)) - max(xdata)))/(width * 2)
    x_arange = [(i-xdiff, i+xdiff) for i in
                sorted(arange(float(min(xdata)), max(xdata) + xdiff, xdiff * 2))]
    graph = ' ' + str(title)
    graph += '\n'
    graph += ' ' + '-' * len(title)
    graph += '\n\n'
    graph += ylabel
    graph += '\n'
    val = 6 - max([len('{0:.0f}'.format(y)) for _, y in y_arange])
    form = '{' + ':<7.{}f'.format(val) + '}+'
    graph += form.format (expf(max(yydata)))
    for yval, (y1, y2) in enumerate(y_arange):
        graph+='\n'
        if not (yval)%5 and yval != 0:
            graph += form.format (expf((y1+y2)/2))
        else:
            graph += ' ' * 7 + '|'
        pos = 0
        for x1, x2 in x_arange:
            for i in xrange(pos, len(yydata)):
                if (y1 < yydata[i] <= y2 and
                    x1 < xdata[i]  <= x2):
                    graph += pch
                    pos += 1
                    break
            else:
                graph += ' '
    graph += '\n'
    if logscale:
        graph += ' 1/inf ' + ''.join(
            ['+' if not x%10 else '-' for x in xrange(width+1)]) + '\n'
    else:
        graph += '     0 ' + ''.join(
            ['+' if not x%10 else '-' for x in xrange(width+1)]) + '\n'
        
    val = 7 - max([len('{0:.0f}'.format(y)) for _, y in x_arange])
    form = '{' + ':<7.{}f'.format(val) + '}  '
    graph += ' '*7 + ''.join(
        [form.format(float(sum(x_arange[x])/2)) for x in xrange(0,width,10)]
    ) + ('' if width % 10 else form.format(float(sum(x_arange[-1])/2)))+ '\n\n'
    graph += ' ' * 7 + '{0:^{1}}'.format(xlabel, width)
    return graph

