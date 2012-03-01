#!/usr/bin/env ruby

require 'pp'
require 'sqlite3'
require 'gnuplot'
require 'optparse'

class Options
	attr_accessor :terminal, :output, :file, :plot_type, :bin_width, :timestep

	@@plot_types = [:energy, :speed, :spherical_density]

	def self.parse
		return_opts = Options.new

		#Parse optional arguments
		opts = OptionParser.new do |opts|
			opts.banner = "Usage: #{$0} database [OPTIONS]\n\t\"database\" must be an sqlite database, with the correct schema."

			opts.on("-t", "--terminal TERMINAL", "Gnuplot terminal to use."){|t|
				return_opts.terminal = t
			}

			opts.on("-o", "--output OUTPUT", "Gnuplot output file. Useful for saving plots."){|o|
				return_opts.output = o
			}

			opts.on("-b", "--bin-width BINWIDTH", "Width of histogram bins, if plotting a histogram.") {|b|
				return_opts.bin_width = b.to_f
			}

			opts.on("-s", "--timestep TIMESTEP", "Plot a snapshot of the system at this timestep.") {|s|
				return_opts.timestep = s.to_f
			}

			opts.on("-p", "--plot-type PLOT_TYPE", 
					"Plot type:",
					"Must be one of \"#{@@plot_types.join(", ")}\"",
				   "Defaults to \"energy\""){|p|

				return_opts.plot_type = p.to_sym
			}
	
			opts.on_tail("-h", "--help", "This help text."){|h|
				puts opts
				exit
			}
		end
		opts.parse!(ARGV)

		#Check a bin width is specified if we are plotting an hgram
		if return_opts.plot_type == :speed && return_opts.bin_width == nil
			STDERR.puts opts
			STDERR.puts "Note that you must specify a bin width when plotting a histogram!"
			exit 1
		end	

		#Check that we have a timestep specified if taking a snapshot
		if return_opts.plot_type == :speed && return_opts.timestep == nil
			STDERR.puts opts
			STDERR.puts "You must specify a time at which to snapshot the system!"
			exit 1
		end
			

		#Now, if we have more than one arg, failboat.
		#The remaining arg should be the database file
		if ARGV.length != 1
			STDERR.puts opts
			exit 1
		end
		return_opts.file = ARGV[0]
		return return_opts
	end
end


def get_energies(db)
	db.results_as_hash = true
	dt = db.execute("SELECT timestep FROM params")[0][0]
	rows = db.execute("SELECT timestep*#{dt} as timestep, kinetic, potential FROM energies")
	return rows.map {|r| {:timestep => r['timestep'], :kinetic => r['kinetic'], :potential => r['potential']} }
end

def get_coord_data(timestep, db, particle_spec, table, limit)
	db.results_as_hash = true
	if particle_spec == nil
		rows = db.execute("SELECT particle_id, x, y, z FROM #{table} WHERE timestep = :ts", :ts => timestep)
	else
		query = "SELECT particle_id, x, y, z FROM #{table} " + 
					"JOIN particles ON particles.id = #{table}.particle_id " +
					"WHERE timestep = :ts AND charge = :q AND mass BETWEEN :m1 AND :m2"

		query << " LIMIT #{limit}" if limit

		rows = db.execute(query,
						:ts => timestep, :q => particle_spec[:charge], 
						:m1 => particle_spec[:mass] - 1E-10,
						:m2 => particle_spec[:mass] + 1E-10)
	end
	return rows.map{|r| {:id => r['particle_id'], table.to_sym => Vector[r['x'], r['y'], r['z']]} }
end

def get_velocities(timestep, db, particle_spec = nil, limit = nil)
	get_coord_data(timestep, db, particle_spec, 'velocities', limit)
end

def get_positions(timestep, db, particle_spec = nil, limit = nil)
	get_coord_data(timestep, db, particle_spec, 'positions', limit)
end

def plot_animation(n, db, gp)
end

def plot_animation(dt, opts, db, gp)
	Gnuplot::SPlot.new(gp) do |plot|
		plot.terminal "pdf"
		plot.output "test.pdf"
		plot.xrange "[0:1]"
		plot.yrange "[0:1]"
		plot.zrange "[0:1]"
		
		timestep = 0
		pos = get_positions(timestep, db)
		#while timestep < 100
			timestep += dt
			positions = []
			pos.each{ |r| positions << r[:positions].to_a}
			splot = positions
			ds = Gnuplot::DataSet.new([1,1,1]){}
			plot.data << ds
		#end
		
	end
end




def plot_spherical_density(bin_width, data, opts, gp)
	bins = {}
	data.each do |row|
		bin = (row[:positions].r / bin_width).floor
		bins.include?(bin) ? bins[bin] += 1 : bins[bin] = 1
	end	
	bins = bins.sort
	radii, densities = bins.map{|r,d| 
		vol = 4.0/3.0 * Math::PI * bin_width**3 * ((r+1)**3 - r**3)
		[r * bin_width + bin_width/2.0, d / vol]
	}.transpose

	Gnuplot::Plot.new(gp) do |plot|
		plot.title "Density vs Radius"
		plot.xlabel "Radius"
		plot.ylabel "Density"
		plot.yrange "[0:*]"
		plot.xrange "[0:*]"
		plot.boxwidth bin_width


		density = Gnuplot::DataSet.new([radii,densities]) do |ds|
			ds.title = "Density"
			ds.with = "boxes"
		end
		plot.data << density
	end
end
def plot_energies(data, opts, gp)
	Gnuplot::Plot.new(gp) do |plot|
		plot.title "Energy vs Time"
		plot.xlabel "Time"
		plot.ylabel "Energy"

		plot.terminal opts.terminal if opts.terminal
		plot.output opts.output if opts.output

		timesteps = data.collect{|rec| rec[:timestep]}
		kinetic   = data.collect{|rec| rec[:kinetic]}
		potential = data.collect{|rec| rec[:potential]}
		total     = data.collect{|rec| rec[:kinetic] + rec[:potential]}

		kinetic = Gnuplot::DataSet.new([timesteps, kinetic]) {|ds|
			ds.with = "lines"
			ds.title = "Kinetic"
		}
		potential = Gnuplot::DataSet.new([timesteps, potential]){|ds|
			ds.with = "lines"
			ds.title = "Potential"
		}
		#The pile of crap in the arguments here is just a coefficient-wise sum
		total = Gnuplot::DataSet.new([timesteps, total]){|ds|
			ds.with = "lines"
			ds.title = "Total"
		}
		plot.data << kinetic << potential << total
	end
end

def plot_speed_dist(bin_width, data, opts, gp)
	bins = {}
	mean = 0
	data.each do |row|
		mean += row[:velocities].r / data.length
		bin = (row[:velocities].r / bin_width).floor
		bins.include?(bin) ? bins[bin] += 1 : bins[bin] = 1
	end	
	bins = bins.sort

	Gnuplot::Plot.new(gp) do |plot|
		plot.title "Speed Distribution"
		plot.xlabel "Speed"
		plot.ylabel "Frequency"
		plot.yrange "[0:*]"
		plot.xrange "[0:*]"
		plot.boxwidth bin_width

		speeds,freqs = bins.to_a.transpose
		#Remap speeds to actual speeds, correct for offcentredness
		speeds.map! {|s| s*bin_width + bin_width / 2} 

		histogram = Gnuplot::DataSet.new([speeds, freqs]){|ds|
			ds.with = "boxes"
			ds.title = "Frequency"
		}

		a = Math::sqrt(Math::PI / 8) * mean
		fit = Gnuplot::DataSet.new(
			"#{data.length * bin_width}*sqrt(2.0 / pi) * x**2 * exp(-x**2 / (2*#{a}**2)) / #{a}**3"){|ds|
			ds.title = "Matched Maxwell Dist."
		}

		plot.data << histogram << fit
	end
end

def plot_collision_angle_rms(db, sample_size, step_size, particle_spec)
	initial_vel_data = get_velocities(0, db, particle_spec, sample_size)

	times = []
	rms_angle = []

	time = step_size
	while true
		rms = 0
		initial_vel_data.each do |r|
			row = db.execute("SELECT x,y,z FROM velocities WHERE particle_id = :id AND timestep = :ts",
							 :id => r[:id], :ts => time)[0]
			curr_vel = Vector[row[0], row[1], row[2]]
			angle = Math::acos((curr_vel/curr_vel.r).inner_product(r[:velocities]/r[:velocities].r))
			#puts angle * 180.0 / Math::PI / 2
			rms += angle*angle / sample_size
		end
		rms = Math::sqrt(rms)
		puts "#{time}\t#{rms * 180.0 / (Math::PI)}"
		times << time
		rms_angle << rms
		time += step_size
	end
end
	
#$VERBOSE = true
Gnuplot::open do |gp|
	opts = Options::parse
	db = SQLite3::Database.new(opts.file)

	case opts.plot_type
	when :positions
		coords = get_coord_data(opts.timestep, db, nil, 'positions', nil)
		coords.each{|rec| puts rec[:positions].to_a.join(",")}	
	when :e_positions
		coords = get_coord_data(opts.timestep, db, {:charge => -1, :mass => 1}, 'positions', nil)
		coords.each{|rec| puts rec[:positions].to_a.join(",")}	
	when :animate
		plot_animation(10, opts, db, gp)
	when :speed
		velocities = get_velocities(opts.timestep, db, {:charge => -1, :mass => 1})
		plot_speed_dist(opts.bin_width, velocities, opts, gp)
	when :spherical_density
		positions = get_positions(opts.timestep, db, {:charge => -1, :mass => 1})
		plot_spherical_density(opts.bin_width, positions, opts, gp)
	when :collision_angle
		plot_collision_angle_rms(db, 50, 10, {:charge => -1, :mass => 1})
	else
		data = get_energies(db)
		plot_energies(data, opts, gp)
	end
end
