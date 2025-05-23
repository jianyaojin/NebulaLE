#ifndef __TRAJECTORIES_CPU_DRIVER_H_
#define __TRAJECTORIES_CPU_DRIVER_H_

#include "cpu_driver.h"

namespace nbl { namespace drivers {

/**
 * \brief CPU driver for tracking energy deposits during scattering.
 * \see cpu_driver
 */

template<
	typename scatter_list_t,
	typename intersect_t,
	typename geometry_manager_t>
class trajectories_cpu_driver
	: public cpu_driver<scatter_list_t, intersect_t, geometry_manager_t>
{
public:
	using base_t = cpu_driver<scatter_list_t, intersect_t, geometry_manager_t>;

	using typename base_t::particle_index_t;
	using typename base_t::material_t;
	using typename base_t::material_manager_t;
	using typename base_t::seed_t;

	/// Constructor
	trajectories_cpu_driver(
		intersect_t intersect,
		material_manager_t const & materials,
		geometry_manager_t const & geometry,
		real min_energy, real max_energy,
		seed_t seed = util::random_generator<false>::default_seed)
	: base_t(intersect, materials, geometry, min_energy, max_energy, seed)
	{}

	/**
	 * \brief Perform a single iteration of the simulation for all particles.
	 *
	 * \param function Callback function to be called if a particle loses energy
	 *                 in a scattering event. Should have signature
	 *                 `void(particle const & before, particle const & after, uint32_t tag)`.
	 */
	template<typename deposit_function>
	void do_iteration(deposit_function function)
	{
		const auto particle_count = this->_particles.get_total_count();
		for (particle_index_t particle_idx = 0; particle_idx < particle_count; ++particle_idx)
		{
			if (!this->_particles.active(particle_idx))
				continue;

			this->init(particle_idx);

			// We want to track only scattering events, not interface events
			bool track = this->_particles.next_scatter(particle_idx);
			const particle before = this->_particles[particle_idx];

			this->intersect(particle_idx);
			this->scatter(particle_idx);

			const particle after = this->_particles[particle_idx];
			if (track && before.kin_energy != after.kin_energy)
				function(before, after,
					this->_particles.get_primary_tag(particle_idx));
		}
	}

	/// Keep simulating until there are no particles left.
	template<typename deposit_function>
	void simulate_to_end(deposit_function&& func)
	{
		for (particle_index_t particle_idx = 0;
			particle_idx < this->_particles.get_total_count(); ++particle_idx)
		{
			while (this->_particles.active(particle_idx))
			{
				this->init(particle_idx);

				// Variables that check for the status of the particle
				bool track_scatter = this->_particles.next_scatter(particle_idx);
				bool track_terminated = this->_particles.terminated(particle_idx);

				const particle before = this->_particles[particle_idx];
				const int parent_edge = static_cast<int>(this->_particles.get_edge_tag(particle_idx));

				// Get the number of particles in the simulation to check later if a secondary was generated.
				const size_t num_before = this->_particles.get_total_count(); 

				this->intersect(particle_idx);
				this->scatter(particle_idx);

				const particle after = this->_particles[particle_idx];
				const int child_edge_one = static_cast<int>(this->_particles.get_edge_tag(particle_idx));
				const size_t num_after = this->_particles.get_total_count(); 

				// If no secondary generated, child_edge_two is equal to child_edge_one. 
				// If secondary generated, child_edge_two will have the value of that edge.
				const int child_edge_two =
    				(num_before == num_after)
      				? child_edge_one
      				: static_cast<int>(this->_particles.get_edge_tag(num_after - 1));

				if (track_scatter)
					func(before, after,
						this->_particles.get_primary_tag(particle_idx),
						parent_edge, child_edge_one, child_edge_two 
					);
				else if (track_terminated)
					func(before, after,
						this->_particles.get_primary_tag(particle_idx),
						parent_edge, -999, -999
					);
			}
		}
	}
};

}} // namespace nbl::drivers

#endif // __TRAJECTORIES_CPU_DRIVER_H_
