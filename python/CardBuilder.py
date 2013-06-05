"""

A mini-language for defining data cards.

Authors: Evan K. Friis, Roger Wolf, Felix Frensch

"""

import collections
import fnmatch
import logging
import os

from HiggsAnalysis.CombinedLimit.DatacardParser import Datacard
from HiggsAnalysis.HiggsToTauTau.TFileCache import get_tfile

log = logging.getLogger(__name__)


def format_signal_mass(mass):
    """ Helper function to format a floating-point signal mass in HCG style

    >>> format_signal_mass(110.0)
    '110'
    >>> format_signal_mass(125.0)
    '125'
    >>> format_signal_mass(125.1)
    '125.1'
    """
    if abs(mass - int(mass)) < 1E-6:
        return str(int(mass))
    else:
        return "%0.1f" % mass


class Process(object):
    """ Represents a specific process (Zjets) in a category

    One can inplace-multiply a nuisance tuple to register
    that nuisance affects this process.

    >>> proc = Process('a_process', False)
    >>> proc.is_signal
    False
    >>> proc.name
    'a_process'
    >>> my_nuisance = Nuisance('a_syst', 'lnN')
    >>> proc.apply_systematic(my_nuisance(1.05))
    >>> proc.systematics[my_nuisance]
    1.05
    """
    def __init__(self, name, is_signal):
        self.name = name
        self.is_signal = is_signal
        self.systematics = {}

    def apply_systematic(self, other):
        syst_name, syst_value = other.systematic, other.value
        if syst_name in self.systematics:
            raise KeyError("Systematic %s already entered into process %s"
                           % (syst_name, self.name))
        self.systematics[syst_name] = syst_value

    def __eq__(self, other):
        """Two processes are only equal if they point to the same object """
        return self is other


class ProcessGroup(object):
    """ Represents a collection of processes

    Multiplex nuisance application across them.  The input arguments are
    flattened.

    >>> proc1 = Process('proc1', False)
    >>> proc2 = Process('proc2', False)
    >>> proc3 = Process('proc3', False)
    >>> procs = ProcessGroup([proc1, proc2], proc3)

    Apply a nuisance uncertainty.

    >>> my_syst = Nuisance('a_syst', 'lnN')
    >>> procs.apply_systematic(my_syst(1.05))

    The nuisance is now applied to both processes.

    >>> proc1.systematics[my_syst]
    1.05
    >>> proc1.systematics[my_syst]
    1.05
    """
    def __init__(self, *processes):

        # http://stackoverflow.com/questions/2158395/\
        # flatten-an-irregular-list-of-lists-in-python
        def flatten(l):
            for el in l:
                if isinstance(el, collections.Iterable) and not isinstance(
                        el, basestring):
                    for sub in flatten(el):
                        yield sub
                else:
                    yield el

        self.processes = list(flatten(processes))

    def apply_systematic(self, syst):
        for proc in self.processes:
            proc.apply_systematic(syst)

    def __iter__(self):
        for proc in self.processes:
            yield proc


_DEFAULT_SHAPE_MAP = {
    '*': ['$CATEGORY/$PROCESS', '$CATEGORY/$PROCESS_$SYSTEMATIC'],
    'qqH': ['$CATEGORY/$PROCESS$MASS', '$CATEGORY/$PROCESS_$SYSTEMATIC'],
    'ggH': ['$CATEGORY/$PROCESS$MASS', '$CATEGORY/$PROCESS_$SYSTEMATIC'],
    'VH': ['$CATEGORY/$PROCESS$MASS', '$CATEGORY/$PROCESS_$SYSTEMATIC'],
    'WH': ['$CATEGORY/$PROCESS$MASS', '$CATEGORY/$PROCESS_$SYSTEMATIC'],
    'WH_htt': ['$CATEGORY/$PROCESS$MASS', '$CATEGORY/$PROCESS_$SYSTEMATIC'],
    'WH_hww': ['$CATEGORY/$PROCESS$MASS', '$CATEGORY/$PROCESS_$SYSTEMATIC'],
    'ZH': ['$CATEGORY/$PROCESS$MASS', '$CATEGORY/$PROCESS_$SYSTEMATIC'],
    'ZH_htt': ['$CATEGORY/$PROCESS$MASS', '$CATEGORY/$PROCESS_$SYSTEMATIC'],
    'ZH_hww': ['$CATEGORY/$PROCESS$MASS', '$CATEGORY/$PROCESS_$SYSTEMATIC'],
}


class Category(object):
    """ Represents a category (bin) in a data card

    One can add signal and background processes like so:

    """

    def __init__(self, name, shape_map=None, data_name='data_obs'):
        self.name = name
        self.signals = []
        self.backgrounds = []
        self.data = data_name
        self.shape_map = shape_map if shape_map else _DEFAULT_SHAPE_MAP

    def __iter__(self):
        """ So the interface is the same as CategoryGroup

        Just yields itself.
        """
        yield self

    def add_signal(self, name):
        """ Add a signal process to this category

        Returns the :class:`Process` created.

        """
        self.signals.append(Process(name, True))
        return self.signals[-1]

    def add_bkg(self, name):
        """ Add a background process to this category

        Returns the :class:`Process` created.

        """
        self.backgrounds.append(Process(name, False))
        return self.backgrounds[-1]

    def __mod__(self, processgroup):
        """ Return the processes in a group which belong to this category """
        if not isinstance(processgroup, ProcessGroup):
            raise TypeError("__contains__ only supports ProcessGroups")
        output = []
        for process in processgroup:
            if process in self.signals or process in self.backgrounds:
                output.append(process)
        return ProcessGroup(*output)

    def build_card_dict(self, mass, shape_file):
        """ Build a dictionary suitable of output using DataCardWriter.

        This whole mess it to take the configuration of the samples + systs,
        figure out the rates for histogram by looking at the shape file, then
        put everything in the same format that is produced by the HCG
        DatacardParser.py tool.

        :param: mass - the mass to substitute into the histogram names
        :param: shape_file - a :class:`ROOT.TFile` containing the histograms
            for this category. They will be used to get the rates for each
            process.
        """

        def get_yield(proc_name):
            # Find the most specific matching shape matcher
            best_key = None
            best_mapping = None
            for key, mapping in self.shape_map.iteritems():
                if fnmatch.fnmatch(proc_name, key):
                    if not best_key or '*' not in key:
                        best_key = key
                        best_mapping = mapping
            path = best_mapping[0].replace('$CATEGORY', self.name).replace(
                '$PROCESS', proc_name).replace('$MASS', str(mass))
            the_histo = get_tfile(shape_file).Get(path)
            if not the_histo:
                raise IOError("Could not find %s for cat: %s proc: %s" %
                              (path, self.name, proc_name))
            return the_histo.Integral()

        # Build the output datacard in the correct form.  Yes, some of these
        # things are redundant, this is a wacky thing. :/
        output = Datacard()
        output.bins = [self.name]
        for signal in self.signals:
            output.isSignal[signal.name] = True
        for bkg in self.backgrounds:
            output.isSignal[bkg.name] = False
        output.signals = [x.name for x in self.signals]
        output.keyline = [(self.name, x.name, True) for x in self.signals] + \
            [(self.name, x.name, False) for x in self.backgrounds]
        output.exp = {
            self.name:  dict((x.name, get_yield(x.name)) for x in
                             self.signals + self.backgrounds)
        }
        output.hasShapes = True
        output.obs = {self.name: get_yield(self.data)}
        output.processes = output.exp[self.name].keys()
        output.shapeMap = {'*': {}}
        for key, mapping in self.shape_map.iteritems():
            output.shapeMap['*'][key] = [
                os.path.basename(shape_file), mapping[0], mapping[1]]
        # Now build systematics - first find the total union of all
        all_systematics = set()
        for proc in self.signals + self.backgrounds:
            all_systematics.update(proc.systematics.keys())
        for systematic in all_systematics:
            # Figure out the effect for each systematic
            process_map = dict(
                (x.name, x.systematics.get(systematic, 0.0))
                for x in self.signals + self.backgrounds)
            output.systs.append((
                systematic.name,
                False,  # no idea what this does
                systematic.pdf,
                [systematic.sideband] if systematic.sideband else [],
                {self.name: process_map}
            ))
        return output


class CategoryGroup(object):
    """ Represents a group of categories.

    The methods are applied to all members.
    """

    def __init__(self, *members):
        self.members = members

    def __iter__(self):
        for member in self.members:
            yield member

    def add_signal(self, name):
        """ Add a signal process to all of the members of this group.

        Returns a :class:`ProcessGroup` of the new processes.

        :param: name - name of the process (i.e. TH1F)
        """
        return [cat.add_signal(name) for cat in self.members]

    def add_bkg(self, name):
        """ Add a background process to all of the members of this group.

        Returns a :class:`ProcessGroup` of the new processes.

        :param: name - name of the process (i.e. TH1F)
        """
        return ProcessGroup(*[cat.add_bkg(name) for cat in self.members])


class Nuisance(object):
    """ Represents a systematic uncertainty

    Keeps track of the name (key), the type, and in case of the gmN, the
    sideband count.

    """
    def __init__(self, name, pdf='lnN', sideband=None):
        self.name = name
        self.pdf = pdf
        self.sideband = sideband
        if pdf == 'gmN' and sideband is None:
            raise ValueError("You must specify a sideband count for gmN systs")

    def __repr__(self):
        if self.sideband is None:
            return "Nuisance('%s', '%s')" % (self.name, self.pdf)
        else:
            return "Nuisance('%s', '%s', %s)" % (
                self.name, self.pdf, repr(self.sideband))

    def __hash__(self):
        # So we can use this as a key in a dictionary
        return hash(repr(self))

    def __call__(self, value=1.):
        """ Scale this systematic by the nuisance uncertainty

        Returns a (:class:`Systematic`, name) which links this nuisance to a
        value.

        >>> procs = Process('fakes', False)
        >>> fake_syst = Nuisance('fake_rate_error', 'lnN')
        >>> fake_syst(1.05)
        Systematic(Nuisance('fake_rate_error', 'lnN'), 1.05)

        This is just sugar so you can do things like:

        >>> fake_syst(1.05) >> procs
        >>> procs.systematics[fake_syst]
        1.05

        """
        return Systematic(self, value)


class Systematic(object):
    """ A systematic is a nuisance with value applied """
    def __init__(self, systematic, value):
        self.systematic = systematic
        self.value = value

    def __repr__(self):
        return 'Systematic(%s, %0.2f)' % (repr(self.systematic), self.value)

    def __rshift__(self, process):
        process.apply_systematic(self)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
