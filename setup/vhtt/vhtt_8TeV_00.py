from HiggsAnalysis.HiggsToTauTau.CardBuilder import Category, CategoryGroup, \
    Nuisance, ProcessGroup

# WH -> LLT channels

NAME = 'vhtt_0_8TeV'

# Define categories
emt = Category('emt')
mmt = Category('mmt')

# Group them for convenience
llt = CategoryGroup(emt, mmt)

CARD = llt

# Register processes
wh = llt.add_signal('WH')
wh_hww = llt.add_signal('WH_hww')

fakes = llt.add_bkg('fakes')
zz = llt.add_bkg('zz')
wz = llt.add_bkg('wz')

# Group different processes together
simulated_samples = ProcessGroup(wh, wh_hww, wz, zz)
signals = ProcessGroup(wh, wh_hww)

# Define systematics
lumi = Nuisance('lumi_8TeV')
eff_e = Nuisance('CMS_eff_e')
eff_m = Nuisance('CMS_eff_m')
eff_t = Nuisance('CMS_eff_t')

scale_t = Nuisance('CMS_scale_t')

fake_b = Nuisance('CMS_fake_b_8TeV')
emt_fakes = Nuisance('CMS_emt_fakes_8TeV')
mmt_fakes = Nuisance('CMS_mmt_fakes_8TeV')

pdf_qqbar = Nuisance('pdf_qqbar')
scale_diboson = Nuisance('QCDscale_VV')
scale_vh = Nuisance('QCDscale_VV')

# apply systematics with different values
lumi(1.04) >> simulated_samples

# This means apply systematic to only simulated_samples which live in emt
eff_e(1.04) >> emt % simulated_samples
eff_m(1.02) >> mmt % simulated_samples
eff_t(1.08) >> simulated_samples

scale_t(1.03) >> simulated_samples

pdf_qqbar(1.033) >> zz
scale_diboson(1.023) >> zz

pdf_qqbar(1.04) >> wz
scale_diboson(1.04) >> wz

pdf_qqbar(1.034) >> signals
scale_vh(1.004) >> signals

emt_fakes(1.15) >> emt % fakes
mmt_fakes(1.38) >> mmt % fakes
