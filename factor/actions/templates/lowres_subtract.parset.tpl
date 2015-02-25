Strategy.InputColumn = {{ incol }}
Strategy.ChunkSize   = 200
Strategy.UseSolver   = F
Strategy.Steps       = [subtract]

Step.subtract.Model.Sources                   = []
Step.subtract.Model.Cache.Enable              = T
Step.subtract.Model.Phasors.Enable            = F
Step.subtract.Model.DirectionalGain.Enable    = F
Step.subtract.Model.Gain.Enable               = T
Step.subtract.Model.Rotation.Enable           = F
Step.subtract.Model.CommonScalarPhase.Enable  = F
Step.subtract.Model.CommonRotation.Enable     = F
Step.subtract.Operation                       = SUBTRACT
Step.subtract.Model.Beam.Enable               = F
Step.subtract.Output.WriteCovariance          = F
Step.subtract.Output.Column                   = {{ outcol }}
