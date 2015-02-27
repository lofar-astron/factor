LoSoTo.Steps = [smooth, norm]
LoSoTo.Solset = [{{ solset }}]
LoSoTo.Soltab = [{{ solset }}/{{ soltab_amp }}, {{ solset }}/{{ soltab_phase }}]

LoSoTo.Steps.smooth.Operation = SMOOTH
LoSoTo.Steps.smooth.FWHM = [{{ smoothing_window }}]
LoSoTo.Steps.smooth.Axes = [time]
LoSoTo.Steps.smooth.Mode = runningmedian

LoSoTo.Steps.norm.Operation = NORM
LoSoTo.Steps.norm.Soltab = [{{ solset }}/{{ soltab_amp }}]
LoSoTo.Steps.norm.NormVal = 1.0
LoSoTo.Steps.norm.NormAxis = time
