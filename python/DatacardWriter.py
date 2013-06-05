"""

The analog of HiggsAnalysis/CombinedLimit/python/DatacardParser.py

Takes in a set of dictionaries in DatacardParser style, and writes an
output .txt data card.

Based on the logic in combineCards.py

"""

import sys

from HiggsAnalysis.CombinedLimit.DatacardParser import isVetoed


def write_datacard(datacards, outfd=sys.stdout, labels=None, column_width=5,
                   channelVetos=tuple(), shape=True, dirname=''):
    cmax = column_width
    obsline = []
    obskeyline = []
    keyline = []
    expline = []
    systlines = {}
    signals = []
    backgrounds = []
    shapeLines = []
    paramSysts = {}
    flatParamNuisances = {}

    for ich, DC in enumerate(datacards):
        label = "ch%d" % (ich + 1)
        # check if labels were specified by hand
        if labels is not None:
            label = labels[ich]

        singlebin = len(DC.bins) == 1
        if label == ".":
            label = DC.bins[0] if singlebin else ""
        elif not singlebin:
            label += "_"
        for b in DC.bins:
            bout = label if singlebin else label + b
            if isVetoed(bout, channelVetos):
                continue
            obskeyline.append(bout)
            # so that we get only DC.processes contributing to this bin
            for (p, e) in DC.exp[b].items():
                if not DC.isSignal[p]:
                    continue
                #print "in DC.exp.items:b,p", b,p
                expline.append("%.4f" % e)
                keyline.append((bout, p, DC.isSignal[p]))
                # so that we get only DC.processes contributing to this bin
            for (p, e) in DC.exp[b].items():
                if DC.isSignal[p]:
                    continue
                #print "in DC.exp.items:b,p", b,p
                expline.append("%.4f" % e)
                keyline.append((bout, p, DC.isSignal[p]))

        # systematics
        for (lsyst, nofloat, pdf, pdfargs, errline) in DC.systs:
            systeffect = {}
            if pdf == "param":
                if lsyst in paramSysts:
                    if paramSysts[lsyst] != pdfargs:
                        raise RuntimeError(
                            "Parameter uncertainty %s mismatch between cards."
                            % lsyst)
                else:
                    paramSysts[lsyst] = pdfargs
                continue
            for b in DC.bins:
                bout = label if singlebin else label + b
                if isVetoed(bout, channelVetos):
                    continue
                if not bout in systeffect:
                    systeffect[bout] = {}
                # so that we get only DC.processes contributing to this bin
                for p in DC.exp[b].keys():
                    r = str(errline[b][p])
                    if type(errline[b][p]) == list:
                        r = "%.3f/%.3f" % (errline[b][p][0], errline[b][p][1])
                    elif type in ("lnN", 'gmM'):
                        r = "%.3f" % errline[b][p]
                    if errline[b][p] == 0:
                        r = "-"
                    # get max col length, as it's more tricky
                    # do do it later with a map
                    if len(r) > cmax:
                        cmax = len(r)
                    systeffect[bout][p] = r
            if lsyst in systlines:
                (otherpdf, otherargs, othereffect, othernofloat) = \
                    systlines[lsyst]
                if otherpdf != pdf:
                    if (pdf == "lnN" and otherpdf.startswith("shape")):
                        if systlines[lsyst][0][-1] != '?':
                            systlines[lsyst][0] += '?'
                        for b, v in systeffect.items():
                            othereffect[b] = v
                    elif (pdf.startswith("shape") and otherpdf == "lnN"):
                        if pdf[-1] != '?':
                            pdf += '?'
                        systlines[lsyst][0] = pdf
                        for b, v in systeffect.items():
                            othereffect[b] = v
                    elif (pdf == otherpdf + "?") or (pdf + "?" == otherpdf):
                        systlines[lsyst][0] = pdf.replace("?", "") + "?"
                        for b, v in systeffect.items():
                            othereffect[b] = v
                    else:
                        raise RuntimeError(
                            "Card #%i defines systematic %s "
                            "as using pdf %s, while a previous file "
                            "defines it as using %s"
                            % (ich, lsyst, pdf, otherpdf))
                else:
                    if pdf == "gmN" and int(pdfargs[0]) != int(otherargs[0]):
                        raise RuntimeError(
                            "Card %i defines systematic %s as using gamma "
                            "with %s events in sideband, while a previous "
                            "file has %s"
                            % (ich, lsyst, pdfargs[0], otherargs[0]))
                    for b, v in systeffect.items():
                        othereffect[b] = v
            else:
                pdfargs = [str(x) for x in pdfargs]
                systlines[lsyst] = [pdf, pdfargs, systeffect, nofloat]
        # flat params
        for K in DC.flatParamNuisances.iterkeys():
            flatParamNuisances[K] = True
        # put shapes, if available
        if len(DC.shapeMap):
            for b in DC.bins:
                bout = label if singlebin else label + b
                if isVetoed(bout, channelVetos):
                    continue
                p2sMap = DC.shapeMap[b] if b in DC.shapeMap else {}
                p2sMapD = DC.shapeMap['*'] if '*' in DC.shapeMap else {}
                for p, x in p2sMap.items():
                    xrep = [xi.replace("$CHANNEL", b) for xi in x]
                    if xrep[0] != 'FAKE' and dirname != '':
                        xrep[0] = dirname + "/" + xrep[0]
                    shapeLines.append((p, bout, xrep))
                for p, x in p2sMapD.items():
                    if p in p2sMap:
                        continue
                    xrep = [xi.replace("$CHANNEL", b) for xi in x]
                    if xrep[0] != 'FAKE' and dirname != '':
                        xrep[0] = dirname + "/" + xrep[0]
                    shapeLines.append((p, bout, xrep))
        elif shape:
            for b in DC.bins:
                bout = label if singlebin else label + b
                shapeLines.append(('*', bout, ['FAKE']))
        # combine observations, but remove line if any of the datacards
        # doesn't have it
        if len(DC.obs) == 0:
            obsline = None
        elif obsline is not None:
            for b in DC.bins:
                bout = label if singlebin else label + b
                if isVetoed(bout, channelVetos):
                    continue
                obsline += [str(DC.obs[b])]

    bins = []
    for (b, p, s) in keyline:
        if b not in bins:
            bins.append(b)
        if s:
            if p not in signals:
                signals.append(p)
        else:
            if p not in backgrounds:
                backgrounds.append(p)

    def write_ln(*xs):
        """Write a line to the output stream"""
        outfd.write(' '.join(xs))
        outfd.write('\n')

    #write_ln("Combination of", "  ".join(args))
    write_ln("imax %d number of bins" % len(bins))
    write_ln("jmax %d number of processes minus 1"
             % (len(signals) + len(backgrounds) - 1))
    write_ln("kmax %d number of nuisance parameters"
             % (len(systlines) + len(paramSysts)))
    write_ln("-" * 130)

    if shapeLines:
        chmax = max([max(len(p), len(c)) for p, c, x in shapeLines])
        cfmt = "%-" + str(chmax) + "s "
        shapeLines.sort(
            lambda x, y: cmp(x[0], y[0]) if x[1] == y[1] else cmp(x[1], y[1]))
        for (process, channel, stuff) in shapeLines:
            write_ln("shapes", cfmt % process,
                     cfmt % channel, ' '.join(stuff))
        write_ln("-" * 130)

    if obsline:
        cmax = max(
            [cmax] + [len(l) for l in obskeyline] + [len(x) for x in obsline])
        cfmt = "%-" + str(cmax) + "s"
        write_ln("bin         ", "  ".join([cfmt % x for x in obskeyline]))
        write_ln("observation ", "  ".join([cfmt % x for x in obsline]))

    write_ln("-" * 130)

    pidline = []
    signals = []
    backgrounds = []
    tmpsignals = []
    for (b, p, s) in keyline:
        if s:
            if p not in tmpsignals:
                tmpsignals.append(p)
    for (b, p, s) in keyline:
        if s:
            if p not in signals:
                signals.append(p)
            pidline.append(signals.index(p) - len(tmpsignals) + 1)
        else:
            if p not in backgrounds:
                backgrounds.append(p)
            pidline.append(1 + backgrounds.index(p))
    cmax = max([cmax] +
               [max(len(p), len(b)) for p, b, s in keyline] +
               [len(e) for e in expline])
    hmax = max([10] +
               [len("%-12s[nofloat]  %s %s" % (l, p, a)) for l, (p, a, e, nf)
                in systlines.items()])
    cfmt = "%-" + str(cmax) + "s"
    hfmt = "%-" + str(hmax) + "s  "
    write_ln(hfmt % "bin",     "  ".join([cfmt % p for p, b, s in keyline]))
    write_ln(hfmt % "process", "  ".join([cfmt % b for p, b, s in keyline]))
    write_ln(hfmt % "process", "  ".join([cfmt % x for x in pidline]))
    write_ln(hfmt % "rate",    "  ".join([cfmt % x for x in expline]))

    write_ln("-" * 130)

    sysnamesSorted = systlines.keys()
    sysnamesSorted.sort()
    for name in sysnamesSorted:
        (pdf, pdfargs, effect, nofloat) = systlines[name]
        if nofloat:
            name += "[nofloat]"
        systline = []
        for b, p, s in keyline:
            try:
                systline.append(effect[b][p])
            except KeyError:
                systline.append("-")
        write_ln(hfmt % ("%-21s   %s  %s" % (name, pdf, " ".join(pdfargs))),
                 "  ".join([cfmt % x for x in systline]))
    for (pname, pargs) in paramSysts.items():
        write_ln("%-12s  param  %s" % (pname, " ".join(pargs)))

    for pname in flatParamNuisances.iterkeys():
        write_ln("%-12s  flatParam" % pname)
