/*
 * DatabaseConnection.h
 *
 *  Created on: 6 Feb 2012
 *      Author: stefans
 */

#ifndef DATABASECONNECTION_H_
#define DATABASECONNECTION_H_

#include <string>
#include <vector>
#include <sqlite3.h>
#include "../../Particle.h"
#include "../../bounds/BoundaryConditions.h"
#include "../../Configuration.h"

namespace treecode {

namespace output {

template <class Vec>
class DatabaseConnection {
public:

	DatabaseConnection(std::string dbname) :
		particles_init_stmt_(NULL),
		particles_pos_stmt_(NULL),
		particles_vel_stmt_(NULL),
		energies_stmt_(NULL){
		sqlite3_open(dbname.c_str(), &db_);
	}

	~DatabaseConnection() {
		sqlite3_finalize(particles_init_stmt_);
		sqlite3_finalize(particles_pos_stmt_);
		sqlite3_finalize(particles_vel_stmt_);
		sqlite3_finalize(energies_stmt_);
		sqlite3_close(db_);
	}

	void write_init_particles(const std::vector<Particle<Vec>*>& parts){
		if(particles_init_stmt_ == NULL)
			sqlite3_prepare_v2(db_, "INSERT INTO particles (id, charge, mass) "
						"VALUES (:id, :charge, :mass)", -1, &particles_init_stmt_, NULL);

		sqlite3_exec(db_, "BEGIN", NULL, NULL, NULL);
		for(Particle<Vec> *p : parts){
			sqlite3_bind_int(particles_init_stmt_, 1, p->getId());
			sqlite3_bind_int(particles_init_stmt_, 2, p->getCharge());
			sqlite3_bind_double(particles_init_stmt_, 3, p->getMass());
			sqlite3_step(particles_init_stmt_);
			sqlite3_reset(particles_init_stmt_);
		}
		sqlite3_exec(db_, "COMMIT", NULL, NULL, NULL);
	}
	void write_timestep_particles(unsigned int timestep, const std::vector<Particle<Vec>*>& parts){
		if(particles_pos_stmt_ == NULL){
			sqlite3_prepare_v2(db_, "INSERT INTO positions (particle_id, timestep, x, y, z) "
				"VALUES (:id, :ts, :x, :y, :z)", -1, &particles_pos_stmt_, NULL);
		}if(particles_vel_stmt_ == NULL){
			sqlite3_prepare_v2(db_, "INSERT INTO velocities (particle_id, timestep, x, y, z) "
				"VALUES(:id, :ts, :x, :y, :z)", -1, &particles_vel_stmt_, NULL);
		}

		sqlite3_exec(db_, "BEGIN", NULL, NULL, NULL);
		for(Particle<Vec> *p : parts){
			sqlite3_bind_int(particles_pos_stmt_, 1, p->getId());
			sqlite3_bind_int(particles_vel_stmt_, 1, p->getId());
			sqlite3_bind_int(particles_pos_stmt_, 2, timestep);
			sqlite3_bind_int(particles_vel_stmt_, 2, timestep);

			for(int i=0;i<p->getPosition().rows();i++){
				double r = p->getPosition()[i];
				double v = p->getVelocity()[i];
				sqlite3_bind_double(particles_pos_stmt_, i+3, r);
				sqlite3_bind_double(particles_vel_stmt_, i+3, v);
			}
			sqlite3_step(particles_pos_stmt_);
			sqlite3_step(particles_vel_stmt_);
			sqlite3_reset(particles_pos_stmt_);
			sqlite3_reset(particles_vel_stmt_);
		}
		sqlite3_exec(db_, "COMMIT", NULL, NULL, NULL);
	}
	void write_sim_params(unsigned int num_parts, const Configuration<Vec>& config, const BoundaryConditions<Vec>& bc){
		sqlite3_stmt* stmt;
		sqlite3_prepare(db_, "INSERT INTO params ("
				"num_particles, length, dimensions, plasma_param, timestep, max_time, force_softening) VALUES "
				"(:num, :L, :D, :gamma, :dt, :max_t, :fs)", -1, &stmt, NULL);
		sqlite3_bind_int(		stmt, 1, num_parts);
		sqlite3_bind_double(	stmt, 2, bc.getSize());
		sqlite3_bind_int(		stmt, 3, config.getNumDimensions());
		sqlite3_bind_double(	stmt, 4, config.getDensity());
		sqlite3_bind_double(	stmt, 5, config.getTimestep());
		sqlite3_bind_double(	stmt, 6, config.getMaxTime());
		sqlite3_bind_double(	stmt, 7, config.getForceSoftening());
		sqlite3_bind_double(	stmt, 8, config.getTheta());
		sqlite3_step(stmt);
		sqlite3_finalize(stmt);
	}

	void writeEnergies(unsigned int timestep, double ke, double pe){
		if(energies_stmt_ == NULL){
			sqlite3_prepare_v2(db_, "INSERT INTO energies (timestep, kinetic, potential)"
					"VALUES (:ts, :ke, :pe);", -1, &energies_stmt_, NULL);
		}

		sqlite3_bind_int(energies_stmt_, 1, timestep);
		sqlite3_bind_double(energies_stmt_, 2, ke);
		sqlite3_bind_double(energies_stmt_, 3, pe);
		sqlite3_step(energies_stmt_);
		sqlite3_reset(energies_stmt_);
	}

	void init_database(unsigned int dims){
		const char *sql =
				"CREATE TABLE params ("
				"num_particles INTEGER, "
				"length REAL, "
				"dimensions INTEGER, "
				"plasma_param REAL, "
				"timestep REAL, "
				"max_time REAL, "
				"force_softening REAL,"
				"theta REAL);"

				"CREATE TABLE particles ("
				"id INTEGER PRIMARY KEY, "
				"charge INTEGER, "
				"mass REAL);"

				"CREATE TABLE positions ("
				"particle_id INTEGER, "
				"timestep INTEGER, "
				"x REAL, y REAL, z REAL);"

				"CREATE TABLE velocities ("
				"particle_id INTEGER,"
				"timestep INTEGER, "
				"x REAL, y REAL, z REAL);"

				"CREATE TABLE energies ("
				"timestep INTEGER, "
				"kinetic REAL, "
				"potential REAL);";
		sqlite3_exec(db_, sql, NULL, NULL, NULL);
	}
	void clear_database(){
		const char *sql =
				"DELETE FROM params;"
				"DELETE FROM particles;"
				"DELETE FROM positions;"
				"DELETE FROM velocities;"
				"DELETE FROM energies";
		sqlite3_exec(db_, sql, NULL, NULL, NULL);
	}

private:
	sqlite3* db_;

	sqlite3_stmt *particles_init_stmt_, *particles_pos_stmt_;
	sqlite3_stmt *particles_vel_stmt_, *energies_stmt_;
};

} /* namespace output */
} /* namespace treecode */
#endif /* DATABASECONNECTION_H_ */
