"""

Cache open TFiles to prevent redundant opening and closing operations.

"""

import ROOT

_OPEN_TFILES = {}


def get_tfile(path):
    return _OPEN_TFILES.setdefault(path, ROOT.TFile.Open(path, 'READ'))
