#ifndef __KIEFT_QUASI_ELASTIC_H_
#define __KIEFT_QUASI_ELASTIC_H_

namespace nbl { namespace scatter {

/**
 * \brief Phonon scattering, Kieft & Bosch's implementation.
 *
 * In the original version of Nebula this is used for elastic scattering at low energies (< 100 eV) but here we separate it and use it as its own
 * scattering channel.
 *
 *   - See doi:10.1088/0022-3727/41/21/215310 (Kieft paper)
 *   - See doi:10.4233/uuid:f214f594-a21f-4318-9f29-9776d60ab06c (Verduin thesis)
 *
 * \tparam gpu_flag             Is the code to be run on a GPU?
 */
template<bool gpu_flag>
class schreiber_phonon
{
public:
	/**
	 * \brief Indicate that this class never generates secondary electrons.
	 */
	constexpr static bool may_create_se = false;

	/**
	 * \brief Print diagnostic info
	 */
	static void print_info(std::ostream& stream)
	{
		stream << std::boolalpha <<
			" * Kieft & Bosch quasi_elastic model\n";
	}

	/**
	 * \brief Get maximal energy, in eV
	 */
	PHYSICS real get_max_energy() const
	{
		return expr(_log_imfp_table.get_scalemax());
	}

	/**
	 * \brief Sample a random free path length
	 */
	template<typename particle_t>
	inline PHYSICS real sample_path(particle_t const & this_particle, util::random_generator<gpu_flag> & rng) const
	{
		// Get inverse mean free path for this kinetic energy
		const real imfp = expr(_log_imfp_table.get(logr(this_particle.kin_energy)));

		// Draw a distance
		return rng.exponential(1 / imfp);
	}

	/**
	 * \brief Perform a scattering event
	 */
	template<typename particle_manager>
	inline PHYSICS void execute(
		particle_manager& particle_mgr,
		typename particle_manager::particle_index_t particle_idx,
		util::random_generator<gpu_flag>& rng) const
	{
		// Retrieve current particle from global memory
		auto this_particle = particle_mgr[particle_idx];

		real cos_theta, sin_theta;
		{// draw a random elastic scatter angle by interpolating tables
			const float x = logr(particle_mgr[particle_idx].kin_energy);
			const float y = rng.unit();
			cos_theta = clampr(_icdf_table.get(x, y), -1, 1);
			sin_theta = sqrtr(1 - cos_theta*cos_theta);
		}

		// normalize current direction
		normalise(this_particle.dir);

		// find a random normal vector to the current direction of flight and normalize
		vec3 normal_dir = normalised(make_normal_vec(this_particle.dir, rng.phi()));

		// determine the scattered direction
		this_particle.dir = this_particle.dir * cos_theta + normal_dir * sin_theta;

		// Lose energy via phonon loss
		this_particle.kin_energy -= minr(50e-3f, _phonon_loss);

		// Store the scattered particle in memory
		particle_mgr[particle_idx] = this_particle;

		// Update the new particle's edge tag and increment the edgetag count in cascades
		particle_mgr.update_edge_tag(particle_idx);

		// Update the counter for scattering events
		particle_mgr.update_scatter_counter(particle_idx);
	}

	/**
	 * \brief Create, given a material file
	 */
	static CPU schreiber_phonon create(hdf5_file const & mat)
	{
		if (!mat.exists("/kieft/quasi_elastic"))
			throw std::runtime_error("Kieft quasi elastic model not found in "
				"material file " + mat.get_filename());

		schreiber_phonon qel;

		{
			qel._log_imfp_table = mat.fill_table1D<real>("/kieft/quasi_elastic/imfp");
			auto K_range = mat.get_log_dimscale("/kieft/quasi_elastic/imfp", 0, qel._log_imfp_table.width());
			qel._log_imfp_table.set_scale(
				(real)std::log(K_range.front()/units::eV), (real)std::log(K_range.back()/units::eV));
			qel._log_imfp_table.mem_scope([&](real* imfp_vector)
			{
				const real unit = real(mat.get_unit("/kieft/quasi_elastic/imfp") * units::nm);
				for (size_t x = 0; x < K_range.size(); ++x)
				{
					imfp_vector[x] = std::log(imfp_vector[x] * unit);
				}
			});
		}

		{
			qel._icdf_table = mat.fill_table2D<real>("/kieft/quasi_elastic/costheta_icdf");
			auto K_range = mat.get_log_dimscale("/kieft/quasi_elastic/costheta_icdf", 0, qel._icdf_table.width());
			auto P_range = mat.get_lin_dimscale("/kieft/quasi_elastic/costheta_icdf", 1, qel._icdf_table.height());
			qel._icdf_table.set_scale(
				(real)std::log(K_range.front()/units::eV), (real)std::log(K_range.back()/units::eV),
				(real)P_range.front(), (real)P_range.back());
			qel._icdf_table.mem_scope([&](real** icdf_vector)
			{
				for (size_t x = 0; x < K_range.size(); ++x)
				for (size_t y = 0; y < P_range.size(); ++y)
				{
					icdf_vector[x][y] = std::max(-1._r, std::min(1._r, icdf_vector[x][y]));
				}
			});
		}

		qel._phonon_loss = static_cast<real>(mat.get_property_quantity("phonon_loss") / units::eV);

		return qel;
	}

	/**
	 * \brief Clone from another instance.
	 */
	template<bool source_gpu_flag>
	static CPU schreiber_phonon create(schreiber_phonon<source_gpu_flag> const & source)
	{
		schreiber_phonon target;

		target._phonon_loss = source._phonon_loss;

		target._log_imfp_table = util::table_1D<real, gpu_flag>::create(source._log_imfp_table);
		target._icdf_table = util::table_2D<real, gpu_flag>::create(source._icdf_table);

		return target;
	}

	/**
	 * \brief Deallocate data held by an instance of this class.
	 */
	static CPU void destroy(schreiber_phonon & qel)
	{
		util::table_1D<real, gpu_flag>::destroy(qel._log_imfp_table);
		util::table_2D<real, gpu_flag>::destroy(qel._icdf_table);
	}

private:
	/**
	 * \brief Table storing the inverse mean free path as function of energy.
	 *
	 * Actually, stores `log(inverse mean free path / nm^-1)` as function of
	 * `log(kinetic energy / eV)`.
	 */
	util::table_1D<real, gpu_flag> _log_imfp_table;

	/**
	 * \brief Table storing the probability distribution for the scattering angle.
	 *
	 * "ICDF" is short for Inverse Cumulative Distribution Function, which is
	 * what this table stores.
	 *
	 * Specifically, stores `cos(theta)` as function of
	 *   - x axis: `log(kinetic energy / eV)`
	 *   - y axis: cumulative probability (between 0 and 1)
	 */
	util::table_2D<real, gpu_flag> _icdf_table;

	real _phonon_loss;  ///< Amount of energy lost in a phonon event (eV)

	template<bool>
	friend class schreiber_phonon;
};

}} // namespace nbl::scatter

#endif // __KIEFT_QUASI_ELASTIC_H_
