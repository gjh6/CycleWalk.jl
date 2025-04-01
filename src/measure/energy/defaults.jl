function update_energy_data!(
	eData::ED,
	partition::LinkCutPartition,
	update::Update{T}
) where {ED <: AbstractEnergyData, T<:Int}
	# @show "in default energy data", typeof(eData)
	# by default, do nothing with the energy data unless the update is specified 
	return
end