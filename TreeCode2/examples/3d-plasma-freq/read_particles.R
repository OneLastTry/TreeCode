#Takes the result of read.table
read_particles <- function(dims, pos=NULL, vel=NULL, skip_last=F, timerange=NULL){
	#Must supply /some/ data...
	if(is.null(pos) && is.null(vel))
		stop("Must supply at least one of pos or vel.")

	#Calculate the numberof timesteps and number of particles
	if(!is.null(pos)){
		max_time <- length(pos[,1])
		num_particles <- length(pos[1,])
	}
	if(!is.null(vel)){
		max_time <- length(pos[,1])
		num_particles <- length(pos[1,])
	}
	#This can be useful if we have trailing whitespace.
	if(skip_last)
		num_particles <- num_particles - 1
	#Default time range to EVERYTHING!
	if(is.null(timerange))
		timerange <- 1:max_time


	timesteps <- list()
	for(t in timerange){
		particles <- list()
		for(i in (1:(floor(num_particles/dims)))){
			#Read in dims entries, add to list, add to particles, add to timesteps
			if(!is.null(pos)){
				particle.pos <- pos[t,seq(i*dims,length=dims)]
				names(particle.pos) <- NULL
			}
			else
				particle.pos <- NULL	
			if(!is.null(vel)){
				particle.vel <- velocities[t,seq(i*dims,length=dims)]
				names(particle.vel) <- NULL
			}
			else
				particle.vel <- NULL

			particles[[i]] <- list(pos=particle.pos, vel=particle.vel)
		}
		timesteps[[t]] <- particles
	}
	return(timesteps)
}

mag <- function(vec){
	return(sqrt( (vec %*% vec)[1,1] ))
}
