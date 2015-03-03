Strategy.InputColumn = {{ incol }}
Strategy.ChunkSize   = 200
Strategy.UseSolver   = F
Strategy.Steps       = [add]

Step.add.Model.Sources                   = [{{ sources }}]
Step.add.Model.Cache.Enable              = T
Step.add.Model.Phasors.Enable            = F
Step.add.Model.DirectionalGain.Enable    = F
Step.add.Model.Gain.Enable               = T
Step.add.Model.Rotation.Enable           = F
Step.add.Model.CommonScalarPhase.Enable  = T
Step.add.Model.CommonRotation.Enable     = F
Step.add.Model.TEC.Enable                = T
Step.add.Operation                       = ADD
Step.add.Model.Beam.Enable               = F
Step.add.Output.WriteCovariance          = F
Step.add.Output.Column                   = {{ outcol }}
