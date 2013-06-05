

# Define categories
eet = Category('eet')
emt = Category('emt')
mmt = Category('mmt')

# Group them for convenience
llt = CategoryGroup(eet, emt, mmt)

# Register processes
wh = llt.add_signal('WH{MASS}')
wh_hww = llt.add_signal('WH_hww{MASS}')

fakes = llt.add_bkg('fakes')
zz = llt.add_bkg('zz')
wz = llt.add_bkg('wz')

llt.add_data('data_obs')

charge_fakes = eet.add_bkg('charge_fakes')

mc_procs = ProcessGroup(wh, wh_hww, wz, zz)

# Define systematics
lumi = Systematic('lumi_8TeV', 'lnN')
mc_procs *= lumi(1.04)


