#!/usr/bin/python

from optparse import OptionParser, OptionGroup
from numpy import recarray
import matplotlib.pylab as plt
import numpy
import string

# Catalogue of planets to observe
def ptable(recarray, columns=None):
    from prettytable import PrettyTable
    x = PrettyTable(recarray.dtype.names)

    if columns is None: recarray.sort(order=recarray.dtype.names)
    else: recarray.sort(order=columns)
    for row in recarray:
        x.add_row(row)
    print x

# Read cvs file exported from Excel
def readStatFile(filename, separator=','):
    import re, string

    data = None
    try:
        fin = open(filename, 'r')
        data = fin.readlines()
        fin.close()
    except: raise RuntimeError('Cannot read file %s' % filename)
    headers = data[0].strip().split(separator)
    units = data[1].strip().split(separator)

    desc = ''
    formatlist = []
    for idx, val in enumerate(data[2].strip().split(separator)):
        if len(val) < 1: val = 0
        try:
            float(val)
            desc = ','.join((desc, 'float'))
            formatlist.append(float)
        except:
            desc = ','.join((desc, 'object'))
            formatlist.append(str)
    desc = numpy.dtype(desc[1:])
    desc.names = headers

    info = numpy.array([])
    for line in data[2:]:
        infovec = line.strip().split(separator)
        infovec = [(element if element else str(0)) for element in infovec]
        if len(info) < 1:
            info = infovec
        else:
            info = numpy.vstack([info, infovec])

    [nrow, ncol] = info.shape
    information = recarray((nrow,),dtype=desc)
    for idx in range(ncol):
        information[desc.names[idx]]=numpy.array(info[:,idx],dtype=formatlist[idx])

    return [units, information]

if __name__ == '__main__':

    parser = OptionParser(usage='python %prog [options] -f <SATS DATA file>.csv', version="%prog 1.0")
    parser.add_option('-f', '--file',
                      dest='infile',
                      type=str,
                      default=None,
                      help='XML file containing data points')

    group = OptionGroup(parser, 'Data Grouping')
    group.add_option('-s', '--sort',
                      dest='sortorder',
                      type=str,
                      default=None,
                      help='Index list for sorting indicators, e.g. "0,1,2"')
    group.add_option("--subregions",
                      dest='subregions',
                      action="store_true",
                      default=False,
                      help="List of subregions in data")
    group.add_option("--regions",
                      dest='regions',
                      action="store_true",
                      default=False,
                      help="List of regions in data")
    group.add_option('--region',
                      dest='region',
                      type=str,
                      default=None,
                      help='Name of region to display data')
    group.add_option('--subregion',
                      dest='subregion',
                      type=str,
                      default=None,
                      help='Name of subregion to display data')
    group.add_option('--colidx',
                      dest='colidx',
                      type=int,
                      default=None,
                      help='Index of variable column selected')
    parser.add_option_group(group)

    group = OptionGroup(parser, 'Column Correlation')
    group.add_option('--ycol',
                      dest='ycol',
                      type=int,
                      default=None,
                      help='Column to plot on the Y-axis')
    group.add_option('--xcol',
                      dest='xcol',
                      type=int,
                      default=None,
                      help='Column to plot on the X-axis')
    parser.add_option_group(group)

    group = OptionGroup(parser, 'Output/Display Options')
    group.add_option("--indicators",
                      dest='indicators',
                      action="store_true",
                      default=False,
                      help="Output indicators to screen")
    group.add_option("--columns",
                      dest='columns',
                      action="store_true",
                      default=False,
                      help="Output headers of variable columns to screen")
    group.add_option("--table",
                      dest='showdata',
                      action="store_true",
                      default=False,
                      help="Display data in file to screen")
    parser.add_option_group(group)

    parser.add_option('-o', "--output",
                      dest='savegraph',
                      action="store_true",
                      default=False,
                      help="Save generated graph to PNG format")
    parser.add_option('-v', "--verbose",
                      dest='verbose',
                      action="store_true",
                      default=False,
                      help="Display graphs")
    (opts, args) = parser.parse_args()

    if opts.infile is None:
        print "Data file must be provided"
        raise SystemExit(parser.print_usage())

    if (opts.xcol is not None and opts.ycol is None) or (opts.xcol is None and opts.ycol is not None):
        raise RuntimeError('Both --xcol and --ycol parameters must be specified')

    [units, data] = readStatFile(opts.infile, separator=';')

    # Sort in order to see grouping in scatter plot
    sortorder = ['Region', 'Sub Region', 'Country']
    if opts.sortorder is not None:
        sortidx = numpy.asarray(opts.sortorder.split(','), dtype=numpy.int).astype(int)
        if numpy.any(numpy.array([string.lower(unit) for unit in numpy.array(units)[sortidx]]) != 'indicator'):
            raise RuntimeError('Non-indicator index selected')
        sortorder = numpy.array(data.dtype.names)[sortidx].tolist()
    data.sort(order=sortorder)

    # Extract selection of data
    if opts.region is not None:
        regionidx = numpy.nonzero(data['Region']==opts.region)[0]
        data = data[regionidx]
    if opts.subregion is not None:
        subregionidx = numpy.nonzero(data['Sub Region']==opts.subregion)[0]
        data = data[subregionidx]
    if opts.colidx is not None:
        varname = data.dtype.names[opts.colidx]
        if string.lower(units[opts.colidx]) == 'indicator':
            raise RuntimeError('"%s" is indicator column, not variable column'%varname)
        print('Showing only variable "%s"'%varname)
        countries = data['Country']
        fig = plt.figure(figsize=[19,11], facecolor='white')
        ax = fig.add_subplot(111)
        ax.set_title(varname)
        ax.set_ylabel(units[opts.colidx])
        ax.tick_params(axis='x', which='major', labelsize=6)
        plt.xticks(range(len(data[varname])), countries, rotation='vertical')
        plt.plot(range(len(data[varname])), data[varname], '.')
        plt.ylim(numpy.min(data[varname])-1,numpy.max(data[varname])+1)
        if opts.verbose: plt.show()
        if opts.savegraph: plt.savefig('Column_%s.png'%'_'.join(varname.split()),dpi=300)
        import sys
        sys.exit(0)


    # Correlator variables
    if opts.xcol is not None and opts.ycol is not None:
        xvar = data.dtype.names[opts.xcol]
        xdata = data[xvar]
        yvar = data.dtype.names[opts.ycol]
        ydata = data[yvar]
        fig = plt.figure(figsize=[19,11], facecolor='white')
        ax = fig.add_subplot(111)
        ax.set_ylabel('%s [%s]' % (yvar, units[opts.ycol]))
        ax.set_xlabel('%s [%s]' % (xvar, units[opts.xcol]))
        plt.plot(xdata, ydata, '.')
        plt.plot(numpy.unique(xdata[ydata!=0]),
                numpy.poly1d(numpy.polyfit(xdata[ydata!=0], ydata[ydata!=0], 1))(numpy.unique(xdata[ydata!=0])), 'r--')
        if opts.verbose: plt.show()
        if opts.savegraph: plt.savefig('%s_vs_%s.png'%('_'.join(xvar.split()), '_'.join(yvar.split())),dpi=300)
        import sys
        sys.exit(0)

    # Display all data from file neatly
    if opts.showdata:
        ptable(data,sortorder)
        import sys
        sys.exit(0)

    # Display indicators and variable
    if opts.indicators:
        import string
        print 'indicators:'
        for idx, header in enumerate(data.dtype.names):
            if string.lower(units[idx]) == 'indicator':
                print '\t %d - %s' % (idx, header)
    if opts.columns:
        print 'columns:'
        for idx, header in enumerate(data.dtype.names):
            if string.lower(units[idx]) != 'indicator':
                print '\t %d - %s' % (idx, header)
    if opts.indicators or opts.columns:
        import sys
        sys.exit(0)

    if opts.regions:
        regions = map(str.strip, data['Region'])
        for region in numpy.unique(regions):
            print region
        import sys
        sys.exit(0)
    if opts.subregions:
        subregions = map(str.strip, data['Sub Region'])
        for subregion in numpy.unique(subregions):
            print subregion
        import sys
        sys.exit(0)


    # Colors per subregion
    subregion = numpy.unique(data['Sub Region'])
    from matplotlib import colors
    clrs=colors.cnames.keys()[::-1]
    clrs=[color for color in clrs if 'dark' in color]
    # Marker numbers per region
    regions = numpy.unique(data['Region'])
    import matplotlib.markers
    mkrs=matplotlib.markers.MarkerStyle.markers.keys()[-len(regions):]
    numbers = numpy.arange(len(regions))+1

    varidx = numpy.nonzero(numpy.array([string.lower(unit) for unit in units])!='indicator')[0]
    variables = numpy.array(data.dtype.names)[varidx]
    countries = data['Country']
    for colidx,name in enumerate(variables):
        print 'Variable %s' % name
        fig = plt.figure(figsize=[19,11], facecolor='white')
        ax = fig.add_subplot(111)
        ax.set_title(name)
        ax.set_ylabel(units[varidx[0]+colidx])
        ax.tick_params(axis='x', which='major', labelsize=6)
        plt.xticks(range(len(data[name])), countries, rotation='vertical')
        plt.hold(True)
        for rowidx, val in enumerate(data[name]):
            clr = clrs[numpy.nonzero(subregion == data['Sub Region'][rowidx])[0][0]%len(clrs)]
            mkr = mkrs[numpy.nonzero(regions == data['Region'][rowidx])[0][0]]
            if opts.region is not None or opts.subregion is not None:
                plt.plot(rowidx, val, marker=mkr, color=clr, mew=1, ms=7)
            else:
                plt.plot(rowidx, val, marker=mkr, color=clr, alpha=0.5, markersize=0.1)
                marker = numbers[numpy.nonzero(regions == data['Region'][rowidx])[0][0]]
                plt.text(rowidx, val, marker, {'fontsize':10}, color=clr)
        plt.hold(False)
        if opts.savegraph:
            if opts.region is not None:
                plt.savefig('Region_%s_Variable_%s.png'%('_'.join(opts.region.split()), '_'.join(name.split())), dpi=300)
            elif opts.subregion is not None:
                plt.savefig('SubRegion_%s_Variable_%s.png'%('_'.join(opts.subregion.split()), '_'.join(name.split())), dpi=300)
            else: plt.savefig('Variable_%s.png'%('_'.join(name.split())), dpi=300)

    if opts.verbose: plt.show()


# -fin-
