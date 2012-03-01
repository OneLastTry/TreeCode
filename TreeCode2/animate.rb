#!/usr/bin/env ruby
require 'erb'
require 'tempfile'

#Argument 0 should be the name of the data file
#Argument 1 should be the initial gnuplot command (eg, plot, splot)
#Argument 2 should be the repeated part of the command (eg,  u 1:2 w l) in ruby templating format
#Argument 3 should be the number of repetition (eg, for 25 particles, there would be 25 sets of coords)
#Argument 4 should be the step size (ie, every <step size> lines)

gnuplot_file = File.open(ARGV[0], 'r')
gnuplot_cmd  = ARGV[1]
template     = ARGV[2]
num_records  = ARGV[3].to_i
step_size    = ARGV[4].to_i
x_min 		 = ARGV[5]
x_max        = ARGV[6]
y_min        = ARGV[7]
y_max        = ARGV[8]
z_min        = ARGV[9]
z_max        = ARGV[10]

file = File.absolute_path('animate.tmp')
cmd_str = "#{gnuplot_cmd} "
templ = ERB.new(template)

num_records.times do |i| 
	b = binding
	cmd_str << templ.result(b)
	cmd_str << ', ' unless i == num_records - 1
end
puts step_size

index = 0
file_index = 0
gnuplot_file.each_line do |line|
	next unless (index+=1) % step_size == 0
	tmp_file = File.new('animate.tmp', 'w')
	tmp_file.puts line
	tmp_file.close
	

	IO.popen("/usr/bin/gnuplot", "w") do |gp|
		gp.puts("set xrange[#{x_min}:#{x_max}]")
		gp.puts("set yrange[#{y_min}:#{y_max}]")
		gp.puts("set zrange[#{z_min}:#{z_max}]")
		gp.puts("set key off")
		gp.puts("set terminal png size 1024,1024")
		gp.puts("set output \"%06d.png\"" % (file_index+=1))
		gp.puts(cmd_str)
		gp.puts("set output")
	end

	print "\r#{file_index}"
end	
puts


