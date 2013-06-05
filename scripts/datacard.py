#!/usr/bin/env python

"""

Create a datacard from a shape file and a card specification file(s).

"""

from RecoLuminosity.LumiDB import argparse
import imp
import logging
import os
import shutil

from HiggsAnalysis.HiggsToTauTau.CardBuilder import Category, CategoryGroup,\
    format_signal_mass
from HiggsAnalysis.HiggsToTauTau.DatacardWriter import write_datacard

log = logging.getLogger(__name__)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('output_dir', metavar='output/dir/path',
                        help="Directory where output cards will go")

    parser.add_argument('shape_file', metavar='input_shapes.root',
                        help="ROOT file containing shapes.")

    parser.add_argument('cards_specs', nargs='+', metavar='cat.py',
                        help="Card setup .pys for each channel")

    parser.add_argument('--masses', type=float, nargs='+',
                        default=[float(x) for x in range(110, 165, 5)],
                        help="Signal masses.  Default: %(default)s")

    parser.add_argument('--label', help="Append label to output")

    parser.add_argument('--verbose', help="Increase logging output")

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING)

    if not os.path.exists(args.shape_file):
        raise SystemExit("Couldn't open shape file %s" % args.shape_file)

    if not os.path.isdir(args.output_dir):
        raise SystemExit("Output path %s is not a directory!"
                         % args.output_dir)

    logging.info("Copying shape file %s => %s")
    shutil.copy(args.shape_file, args.output_dir)
    shape_file_path = os.path.join(args.output_dir,
                                   os.path.basename(args.shape_file))

    for card in args.cards_specs:
        modname = os.path.basename(card).replace('.py', '')
        log.info("Loading card specification from %s", card)
        module = imp.load_source(modname, card)

        if not hasattr(module, 'NAME'):
            raise SystemExit("Card @ %s doesn't have a NAME definition" % card)

        if not hasattr(module, 'CARD'):
            raise SystemExit("Card @ %s doesn't have a CARD definition" % card)

        card_name = module.NAME

        card = module.CARD

        if not isinstance(card, (Category, CategoryGroup)):
            raise SystemExit("The CARD variable must point to a Category or"
                             " CategoryGroup instance")

        log.info("Creating data card: %s", card_name)

        for mass in args.masses:
            card_dicts = []
            for category in card:
                log.info("Building category: %s", category.name)
                card_dicts.append(category.build_card_dict(
                    format_signal_mass(mass), shape_file_path))

            card_path = os.path.join(
                args.output_dir,
                ''.join([card_name, '-', format_signal_mass(mass), '.txt']))

            with open(card_path, 'w') as cardfd:
                write_datacard(card_dicts, cardfd)
