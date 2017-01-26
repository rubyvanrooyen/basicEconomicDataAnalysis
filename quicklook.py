#!/usr/bin/python

from optparse import OptionParser
from numpy import recarray
import numpy
import matplotlib.pylab as plt

# Catalogue of planets to observe
def ptable(recarray, columns=None):
    from prettytable import PrettyTable
    x = PrettyTable(recarray.dtype.names)

    if columns is None: recarray.sort(order=recarray.dtype.names)
    else: recarray.sort(order=columns)
    for row in recarray:
        x.add_row(row)
    print x

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
    parser.add_option("--countries",
                      dest='countries',
                      action="store_true",
                      default=False,
                      help="Graphs showing all countries")
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
# python quicklook.py -f ESF_DATA_ORIG.csv --countries -o -v

    if opts.infile is None:
        print "Data file must be provided"
        raise SystemExit(parser.print_usage())

    [units, data] = readStatFile(opts.infile, separator=';')
    # ptable(data,['Sub Region', 'Region', 'Country'])

    # Colors per subregion
    subregion = numpy.unique(data['Sub Region'])
    from matplotlib import colors
    clrs=colors.cnames.keys()#[-len(subregion):]
    clrs=[color for color in clrs if 'dark' in color]

    # Marker numbers per region
    regions = numpy.unique(data['Region'])
    import matplotlib.markers
    mkrs=matplotlib.markers.MarkerStyle.markers.keys()[-len(regions):]
    numbers = numpy.arange(len(regions))+1

    # Sort in order to see grouping in scatter plot
    data.sort(order=['Region', 'Sub Region', 'Country'])

    # Display graphs of variables for all countries
    if opts.countries:
        countries = data['Country']
        for colidx,name in enumerate(data.dtype.names[3:]):
            fig = plt.figure(figsize=[19,11], facecolor='white')
            ax = fig.add_subplot(111)
            ax.set_title(name)
            ax.set_ylabel(units[3+colidx])
            ax.tick_params(axis='x', which='major', labelsize=6)
            plt.xticks(range(len(data[name])), countries, rotation='vertical')
            plt.hold(True)
            for rowidx, val in enumerate(data[name]):
                clr = clrs[numpy.nonzero(subregion == data['Sub Region'][rowidx])[0][0]%len(clrs)]
                mkr = mkrs[numpy.nonzero(regions == data['Region'][rowidx])[0][0]]
                plt.plot(rowidx, val, marker=mkr, color=clr, alpha=0.3,
                        markersize=0.1)
                marker = numbers[numpy.nonzero(regions == data['Region'][rowidx])[0][0]]
                plt.text(rowidx, val, marker, {'fontsize':10}, color=clr)
            plt.hold(False)
            if opts.savegraph: plt.savefig('_'.join(name.split())+'.png')

    # data.sort(order=['Sub Region', 'Region', 'Country'])
    # # data.sort(order=['Region', 'Sub Region', 'Country'])
    # countries = data['Country']
    # for colidx,name in enumerate(data.dtype.names[3:]):
    #     for rowidx, val in enumerate(data[name]):
    #         clr = clrs[numpy.nonzero(subregion == data['Sub Region'][rowidx])[0][0]%len(clrs)]
    #         mkridx = numpy.nonzero(regions == data['Region'][rowidx])[0][0]
    #         mkr = mkrs[mkridx]
    #         fig = plt.figure(colidx*(mkridx+1), figsize=[19,11], facecolor='white')
    #         ax = fig.add_subplot(111)
    #         ax.set_title(name)
    #         ax.set_ylabel(units[3+colidx])
    #         ax.tick_params(axis='x', which='major', labelsize=6)
    #         plt.xticks(range(len(data[name])), countries, rotation='vertical')
    #         plt.hold(True)
    #         plt.plot(rowidx, val, marker=mkr, color=clr, alpha=0.3,
    #                 markersize=0.1)
    #         marker = numbers[numpy.nonzero(regions == data['Region'][rowidx])[0][0]]
    #         plt.text(rowidx, val, marker, {'fontsize':10}, color=clr)
    #         plt.hold(False)


    if opts.verbose: plt.show()


# -fin-
