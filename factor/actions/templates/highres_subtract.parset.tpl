Strategy.InputColumn = {{ incol }}
Strategy.ChunkSize   = 200
Strategy.UseSolver   = F
Strategy.Steps       = [subtract, correct1]

Step.subtract.Model.Sources                   = [{{ sources }}]
Step.subtract.Model.Cache.Enable              = T
Step.subtract.Model.Phasors.Enable            = F
Step.subtract.Model.DirectionalGain.Enable    = F
Step.subtract.Model.Gain.Enable               = T
Step.subtract.Model.Rotation.Enable           = F
Step.subtract.Model.CommonScalarPhase.Enable  = F
Step.subtract.Model.CommonRotation.Enable     = F
Step.subtract.Operation                       = SUBTRACT
Step.subtract.Model.Beam.Enable               = F
Step.subtract.Output.Column                   = {{ outcol1 }}

Step.correct1.Model.Sources                 = []
Step.correct1.Model.DirectionalGain.Enable  = F
Step.correct1.Model.Gain.Enable             = T
Step.correct1.Model.Phasors.Enable          = F
Step.correct1.Operation                     = CORRECT
Step.correct1.Output.Column                 = {{ outcol2 }}
Step.correct1.Model.Beam.Enable             = F
Step.correct1.Output.WriteCovariance        = F
Step.correct1.Model.Beam.UseChannelFreq     = F
