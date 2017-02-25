#!/usr/bin/python

from optparse import OptionParser, OptionGroup
from numpy import recarray
import matplotlib.pylab as plt
import numpy
import string
import os

other_regions={
               'African Space Nations':['Algeria', 'Angola', 'Egypt', 'Kenya', 'Morocco', 'Nigeria', 'South Africa', 'Tunisia'],
               'SKA Africa':['Botswana', 'Ghana', 'Kenya', 'Madagascar', 'Mauritius', 'Mozambique', 'Namibia', 'South Africa', 'Zambia'],
               'SKA Organisation':['Australia', 'Canada', 'China', 'India', 'Italy', 'Netherlands', 'New Zealand', 'South Africa', 'Sweden', 'United Kingdom'],
              }

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
        if not any(infovec[infovec=='']):
            continue # empty line
        infovec = [(element.strip() if element else str(0)) for element in infovec]
        if len(info) < 1:
            info = infovec
        else:
            info = numpy.vstack([info, infovec])

    [nrow, ncol] = info.shape
    information = recarray((nrow,),dtype=desc)
    for idx in range(ncol):
        information[desc.names[idx]]=numpy.array(info[:,idx],dtype=formatlist[idx])

    return [units, information]

# Extract relevant region
def get_region(data, region):
    regionidx = numpy.nonzero(data['Region']==region)[0]
    if len(regionidx)<1:
        regionidx = []
        if other_regions.has_key(region):
            for cntry in other_regions[region]:
                if len(numpy.nonzero([cntry in s for s in data['Country']])[0]) < 1:
                    print('Country %s not in data set, or check spelling' % cntry)
                    continue
                regionidx.append(numpy.nonzero([cntry in s for s in data['Country']])[0][0])
        else:
            raise RuntimeError('Region %s not classified' % region)
    regionidx = numpy.array(regionidx)
    return data[regionidx]

# Plot functions
def plot_colidx(units, data, colidx, label=None):
    varname = data.dtype.names[opts.colidx]
    if string.lower(units[opts.colidx]) == 'indicator':
        raise RuntimeError('"%s" is indicator column, not variable column'%varname)
    print('Showing only variable "%s"'%varname)
    countries = data['Country']
    ax.set_title(varname)
    ax.set_ylabel(units[opts.colidx])
    ax.tick_params(axis='x', which='major', labelsize=6)
    plt.xticks(range(len(data[varname])), countries, rotation='vertical')
    if label is None:
        plt.plot(range(len(data[varname])), data[varname], '.')
    if label is not None:
        plt.plot(range(len(data[varname])), data[varname], '.', label=label)
    plt.ylim(numpy.min(data[varname])-1,numpy.max(data[varname])+1)
    return varname
def plot_xcorr(units, data, xcol, ycol, label=None):
    xvar = data.dtype.names[xcol]
    xdata = data[xvar]
    yvar = data.dtype.names[ycol]
    ydata = data[yvar]
    ax.set_ylabel('%s [%s]' % (yvar, units[ycol]))
    ax.set_xlabel('%s [%s]' % (xvar, units[xcol]))
    if label is None:
        plt.plot(xdata, ydata, '.')
        plt.plot(numpy.unique(xdata[ydata!=0]),
                numpy.poly1d(numpy.polyfit(xdata[ydata!=0], ydata[ydata!=0], 1))(numpy.unique(xdata[ydata!=0])), 'r--')
    if label is not None:
        plt.plot(xdata, ydata, '.', label=label)
    return [xvar, yvar]
def plot_regions(units, data, name, alpha=1., label=None, all_flag=False):
    # Colors per subregion
    subregions = numpy.unique(data['Sub Region'])
    from matplotlib import colors
    clrs=colors.cnames.keys()[::-1]
    clrs=[color for color in clrs if 'dark' in color][:len(subregions)]
    # Marker numbers per region
    regions = numpy.unique(data['Region'])
    import matplotlib.markers
    mkrs=matplotlib.markers.MarkerStyle.markers.keys()[-len(regions):]
    numbers = numpy.arange(len(regions))+1

    for rowidx, val in enumerate(data[name]):
        if label is not None and rowidx < 1:
            if not all_flag:
                plt.plot(rowidx, val, marker=mkrs[rowidx%len(mkrs)], color=clrs[rowidx%len(clrs)], mew=1, ms=7, alpha=alpha, label=label)
            else:
                plt.plot(rowidx, val, marker=mkrs[rowidx%len(mkrs)], color=clrs[rowidx%len(clrs)], markersize=0.1, alpha=alpha, label=label)
                marker = numbers[numpy.nonzero(regions == data['Region'][rowidx])[0][0]]
                plt.text(rowidx, val, marker, {'fontsize':10}, color=clrs[rowidx%len(clrs)])
        else:
            if not all_flag:
                plt.plot(rowidx, val, marker=mkrs[rowidx%len(mkrs)], color=clrs[rowidx%len(clrs)], mew=1, ms=7, alpha=alpha)
            else:
                plt.plot(rowidx, val, marker=mkrs[rowidx%len(mkrs)], color=clrs[rowidx%len(clrs)], markersize=0.1, alpha=alpha)
                marker = numbers[numpy.nonzero(regions == data['Region'][rowidx])[0][0]]
                plt.text(rowidx, val, marker, {'fontsize':10}, color=clrs[rowidx%len(clrs)])


if __name__ == '__main__':

    parser = OptionParser(usage='python %prog [options] <SATS DATA file>.csv', version="%prog 1.0")
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

    if len(args) < 1:
        print "At least one data file must be provided"
        raise SystemExit(parser.print_usage())

    if (opts.xcol is not None and opts.ycol is None) or (opts.xcol is None and opts.ycol is not None):
        raise RuntimeError('Both --xcol and --ycol parameters must be specified')

    if len(args) > 1:
        data_dict = {}
        for infile in args:
            data_date = os.path.splitext(os.path.basename(infile))[0].split('_')[-1]
            [units, data] = readStatFile(infile, separator=';')
            data_dict[data_date]={'units' : units, 'data' : data}
    else:
        [units, data] = readStatFile(args[0], separator=';')

    # Sort in order to see grouping in scatter plot
    sortorder = ['Region', 'Sub Region', 'Country']
    if opts.sortorder is not None:
        if len(args) > 1:
            units = data_dict[data_dict.keys()[0]]['units']
            data = data_dict[data_dict.keys()[0]]['data']
        sortidx = numpy.asarray(opts.sortorder.split(','), dtype=numpy.int).astype(int)
        if numpy.any(numpy.array([string.lower(unit) for unit in numpy.array(units)[sortidx]]) != 'indicator'):
            raise RuntimeError('Non-indicator index selected')
        sortorder = numpy.array(data.dtype.names)[sortidx].tolist()
    data.sort(order=sortorder)

    # Extract selection of data
    if opts.region is not None:
        if len(args) > 1:
            units = data_dict[data_dict.keys()[0]]['units']
            data = data_dict[data_dict.keys()[0]]['data']
            for data_date in data_dict.keys():
                data_dict[data_date]['data'] = get_region(data_dict[data_date]['data'], opts.region)
        else:
            data = get_region(data, opts.region)

    if opts.subregion is not None:
        if len(args) > 1:
            units = data_dict[data_dict.keys()[0]]['units']
            data = data_dict[data_dict.keys()[0]]['data']
        subregionidx = numpy.nonzero(data['Sub Region']==opts.subregion)[0]
        if len(args) > 1:
            for data_date in data_dict.keys():
                data_dict[data_date]['data'] = data_dict[data_date]['data'][subregionidx]
        else:
            data = data[subregionidx]

    if opts.colidx is not None:
        fig = plt.figure(figsize=[19,11], facecolor='white')
        ax = fig.add_subplot(111)
        if len(args) > 1:
            for data_date in data_dict.keys():
                units = data_dict[data_date]['units']
                data = data_dict[data_date]['data']
                varname = plot_colidx(units, data, opts.colidx, label=data_date)
            ax.legend(loc=0)
        else:
            varname = plot_colidx(units, data, opts.colidx)
        if opts.verbose: plt.show()
        if opts.savegraph: plt.savefig('Column_%s.png'%'_'.join(varname.split()),dpi=300)
        import sys
        sys.exit(0)


    # Correlator variables
    if opts.xcol is not None and opts.ycol is not None:
        fig = plt.figure(figsize=[19,11], facecolor='white')
        ax = fig.add_subplot(111)
        if len(args) > 1:
            for data_date in data_dict.keys():
                units = data_dict[data_date]['units']
                data = data_dict[data_date]['data']
                [xvar, yvar] = plot_xcorr(units, data, opts.xcol, opts.ycol, label=data_date)
            ax.legend(loc=0)
        else:
            [xvar, yvar] = plot_xcorr(units, data, opts.xcol, opts.ycol)
        if opts.verbose: plt.show()
        if opts.savegraph: plt.savefig('%s_vs_%s.png'%('_'.join(xvar.split()), '_'.join(yvar.split())),dpi=300)
        import sys
        sys.exit(0)

    # Display all data from file neatly
    if opts.showdata:
        if len(args) > 1:
            raise RuntimeError('Display table functionality only available for single input file')
        ptable(data,sortorder)
        import sys
        sys.exit(0)

    # Display indicators and variable
    if opts.indicators or opts.columns:
        if len(args) > 1:
            raise RuntimeError('Variable display functionality only available for single input file')
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

    # Display indicators and variable
    if opts.regions or opts.subregions:
        if len(args) > 1:
            units = data_dict[data_dict.keys()[0]]['units']
            data = data_dict[data_dict.keys()[0]]['data']
    if opts.regions:
        regions = map(str.strip, data['Region'])
        for region in numpy.unique(regions):
            print region
    if opts.subregions:
        subregions = map(str.strip, data['Sub Region'])
        for subregion in numpy.unique(subregions):
            print subregion
    if opts.regions or opts.subregions:
        import sys
        sys.exit(0)


    if len(args) > 1:
        units = data_dict[data_dict.keys()[0]]['units']
        data = data_dict[data_dict.keys()[0]]['data']
    varidx = numpy.nonzero(numpy.array([string.lower(unit) for unit in units])!='indicator')[0]
    variables = numpy.array(data.dtype.names)[varidx]
    countries = data['Country']

    all_flag=True
    if opts.region is not None or opts.subregion is not None:
        all_flag=False
    for colidx,name in enumerate(variables):
        print 'Variable %s' % name
        try: data[name]
        except: continue # skip over heading you cannot deal with
        fig = plt.figure(figsize=[19,11], facecolor='white')
        ax = fig.add_subplot(111)
        ax.set_title(name)
        ax.set_ylabel(units[varidx[0]+colidx])
        ax.tick_params(axis='x', which='major', labelsize=6)
        plt.xticks(range(len(data[name])), countries, rotation='vertical')
        if len(args) > 1:
            plt.hold(True)
            for date_idx, data_date in enumerate(numpy.sort(data_dict.keys())):
                units = data_dict[data_date]['units']
                data = data_dict[data_date]['data']
                plot_regions(units, data, name,
                        alpha=float(date_idx+1)/float(len(data_dict.keys())),
                        label=data_date, all_flag=all_flag)
            plt.hold(False)
            ax.legend(loc=0)
        else:
            plot_regions(units, data, name, all_flag=all_flag)
        if opts.savegraph:
            if opts.region is not None:
                plt.savefig('Region_%s_Variable_%s.png'%('_'.join(opts.region.split()), '_'.join(name.split())), dpi=300)
            elif opts.subregion is not None:
                plt.savefig('SubRegion_%s_Variable_%s.png'%('_'.join(opts.subregion.split()), '_'.join(name.split())), dpi=300)
            else: plt.savefig('Variable_%s.png'%('_'.join(name.split())), dpi=300)

    if opts.verbose: plt.show()


# -fin-
