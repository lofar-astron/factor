Strategy.InputColumn = {{ incol }}
Strategy.ChunkSize   = 200
Strategy.UseSolver   = F
Strategy.Steps       = [correct]

Step.correct.Model.Sources                 = []
Step.correct.Model.CommonScalarPhase.Enable= F
Step.correct.Model.Cache.Enable            = T
Step.correct.Model.DirectionalGain.Enable  = F
Step.correct.Model.Gain.Enable             = T
Step.correct.Model.Phasors.Enable          = F
Step.correct.Operation                     = CORRECT
Step.correct.Output.Column                 = {{ outcol }}
Step.correct.Model.Beam.Enable             = F
Step.correct.Output.WriteCovariance        = F
