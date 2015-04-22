Strategy.InputColumn = {{ incol }}
Strategy.ChunkSize   = 200
Strategy.UseSolver   = F
Strategy.Steps       = [subtract, correct]

Step.subtract.Model.Sources                   = [{{ sources }}]
Step.subtract.Model.Cache.Enable              = T
Step.subtract.Model.Phasors.Enable            = F
Step.subtract.Model.DirectionalGain.Enable    = F
Step.subtract.Model.Gain.Enable               = T
Step.subtract.Model.Rotation.Enable           = F
Step.subtract.Model.CommonScalarPhase.Enable  = T
Step.subtract.Model.TEC.Enable                = T
Step.subtract.Model.CommonRotation.Enable     = F
Step.subtract.Operation                       = SUBTRACT
Step.subtract.Model.Beam.Enable               = F
Step.subtract.Output.WriteCovariance          = F
Step.subtract.Output.Column                   = {{ outcol }}

Step.correct.Model.Sources                 = []
Step.correct.Model.DirectionalGain.Enable  = F
Step.correct.Model.Gain.Enable             = T
Step.correct.Model.Phasors.Enable          = F
Step.correct.Operation                     = CORRECT
Step.correct.Output.Column                 = {{ outcol2 }}
Step.correct.Model.Beam.Enable             = F
Step.correct.Output.WriteCovariance        = F
Step.correct.Model.Beam.UseChannelFreq     = F
