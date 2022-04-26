function bn(agent::AbstractAgent, model::ABM)
 return false
end

function optphens(site::Tuple{Int, Int}, model::ABM)
  return [1.5]
end

function optphens2(site::Tuple{Int,Int}, model::ABM)
  opt1 = [0.8214794136627462, 0.49335401288317304, 0.3435556249656141, 0.25050033180075226]
  opt2 = opt1 .+ 0.1
  ss = EvoDynamics.coord2vertex(site, model)
  if model.step[1] > 6
    return [opt1[ss]]
  else
    return [opt2[ss]]
  end
end

function optphens3(site::Tuple{Int,Int}, model::ABM)
  phenotypic_values = [[0.7614758101208934 2.3911353663016586; 2.2091361414343313 1.7858661540288683; 0.7920974352892358 0.7630236717263839; 2.587205421882114 2.311211631439866],
    [0.6437641305315445 2.097629262577973; 0.9954545836033715 1.1314248848669466; 2.469792530385348 1.490299526620522; 1.6158867433451882 0.2566056477862022]]
  agent_site = (site[2] - 1) * size(model.space)[2] + (site[1])
  if model.step[1] > 7
    return phenotypic_values[1][agent_site, :]
  else
    return phenotypic_values[2][agent_site, :]
  end
end

env_resources(time::Int) = [200]
env_resources2(time::Int) = [200, 220, 183, 190]