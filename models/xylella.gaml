/**
* Name: Xylella MAS
* Author: Elisa Fasanelli, Nicolò Gozzi, Sarah Perrone
* Description: Multi-Agent Simulation of xyella bacterium diffusion in Apulia
* Tags: epidemiology, xylella
*/

model xylella


global {
	
	
	/* 
	 * 
	 * Variables
	 * 
	 */
	 
	//path to shapefile
	file shape_file_olive <- file("../includes/shp/shp_ulivi/ulivi_le_ta_br_95_centroidi.shp");
	file shape_file_lecce <- file("../includes/shp/shp_puglia/lecce/lecce_poli.shp");
	file shape_file_taranto <- file("../includes/shp/shp_puglia/taranto/taranto_poli.shp");
	file shape_file_brindisi <- file("../includes/shp/shp_puglia/brindisi/brindisi_poli.shp");
	file shape_file_boundaries <- file("../includes/shp/shp_puglia/puglia/puglia_meridionale_boundaries.shp");
	
	//bounds within xylella swarms move and interact
	geometry bounds;
	geometry shape <- envelope(shape_file_boundaries);
	
	//if true log files are generated
	bool do_log <- true;	
	
	//max speed of swarms (i.e. the maximum distance covered daily)
	float max_speed <- 100.0;
	
	//maximum radius of swarms
	float max_radius <- 400.0;
	
	//total number of xylella swarms
	int nb_xylella <- 10000;
	
	//probability of swarms delocalization
	float relocation_probability <- 0.99;
	
	//amplitude of swarm movement
	float movement_amplitude <- 22.5;
	
	//lenght of the simulation in cycles (1cycle = 1day)
	int length_simulation <- 1826;
	
	//year when start cutting procedure 
	int start_cutting_year <- 2013;
	
	//coordinates of initial infected zone
	float x_start_inf_area <- 112945.0;
	float y_start_inf_area <- 90560.0;
	float radius_start_inf_area <- 15000.0;
	
	//path to log files
	string id <- string(machine_time);
	string path_to_log_xylella <- "../logs/log_xylella" + id + ".csv";
	string path_to_log_olive <- "../logs/log_olive" + id + ".csv"; 
	string path_to_log_olive_total <- "../logs/log_olive_total" + id + ".csv";
	string path_to_log_infection_network <- "../logs//log_infection" + id + ".csv";
	
	//starting date, season and current date
	date start <- date([2010,1,1]);
	int season <- 1;
	date today <- start;
	
	//number of olive agents checked during a single checking period
	int olive_checked <- 0;
	
	//number of olive agents cut during a single checking period
	int olive_cutted <- 0;
	
	//percentage of infected (xylella and olive) in different provinces
	float perc_infected_lecce_x <- 0.0;
	float perc_infected_taranto_x <- 0.0;
	float perc_infected_brindisi_x <- 0.0; 
	float perc_infected_lecce_o <- 0.0;
	float perc_infected_taranto_o <- 0.0;
	float perc_infected_brindisi_o <- 0.0; 
	float mean_radius_lecce <- 0.0;
	float mean_radius_brindisi <- 0.0;
	float mean_radius_taranto <- 0.0;

	//infective power of xylella bacterium
	float beta <- 0.0;
	
	//maps of wind directions
	map<int,float> headings <- [0::11.25,
								1::33.75,
								2::56.25,
								3::78.75,
								4::101.25,
								5::123.75,
								6::146.25,
								7::168.75,
								8::191.25,
								9::213.75,
								10::236.25,
								11::258.75,
								12::281.25,
								13::303.75,
								14::326.25,
								15::348.75
								];
								
	//maps of swarms density in different months			
	map<int,float> density_map <- [1::4/9,
							   	   2::5/9,
							   	   3::6/9,
							       4::7/9,
							       5::8/9,
							       6::1.0,
							       7::8/9,
							       8::7/9,
							       9::6/9,
						           10::5/9,
						           11::4/9,
						           12::3/9
							      ];
	
	//network of infection contacts
	graph infection_graph <- graph<string, string>([]);
	
	//degree map
	map<string,float> degrees_map;
	
	//if true is performed a cut proportional to degree
	float cut_hubs <- true;
	
	//define the cutoff percentile
	float cutoff_perc <- 0.9;
	
	//current cut_off 
	float cut_off <- 0.0;
	
	
	/*
	 * 
	 * Initialization
	 * 
	 */
	 
	init {

		//initialize beta
		beta <- get_beta(today.month);
		
		//initialize log files 
		if do_log {
			
			//xylella
			save ("cycle;nb;mean_perc_sane;mean_perc_infected;mean_perc_infected_le;mean_perc_infected_ta;mean_perc_infected_br;tot_area")
		   		to: path_to_log_xylella type: "text";
		   		
		   		
		   	//olive
			save ("cycle;mean_perc_sane;mean_perc_infected;mean_perc_infected_le;mean_perc_infected_ta;mean_perc_infected_br;mean_radius_le;mean_radius_ta;mean_radius_br;mean_radius;total_area")
		   		to: path_to_log_olive type: "text";
		   		
		    //olive (complete)
		    save ("cycle;name;perc_sane;perc_infected;radius")
		   		to: path_to_log_olive_total type: "text";
		   		
		   	//infections 
		   	save("cycle;from;to;perc_new_infected")
		   		to: path_to_log_infection_network type: "text";
		}
		
		//create province agents from shp files
		//lecce
		create province from: shape_file_lecce {
			file csv <- csv_file('../includes/frequenza/lecce_dict.csv');
			wind_direction_frequency <- matrix(csv);
			name <- "lecce";
		}
		
		//taranto
		create province from: shape_file_taranto {
			file csv <- csv_file('../includes/frequenza/taranto_dict.csv');
			wind_direction_frequency <- matrix(csv);
			name <- "taranto";
		}
		
		//brindisi
		create province from: shape_file_brindisi {
			file csv <- csv_file('../includes/frequenza/brindisi_dict.csv');
			wind_direction_frequency <- matrix(csv);
			name <- "brindisi";
		}
		
		//create olive agents from shp file
		create olive from: shape_file_olive 
			with: [area::float(read('SHAPE_AREA')), 
				   name::string(read('orig_ogc_f')), 
				   color::#green, 
				   perc_infected::0.0, 
				   perc_sane::1.0,
				   prob_checking:: rnd(1.0)
			];
			
		//initialize location province of each olive agent
		ask olive {
			ask province {
				if inter(self, myself) != nil {
					myself.in_province <- self;
				}
			}	
		}
		
		//initialize bounds for swarms
		bounds <- union(province);
		
		//create xylella agents
		create xylella number: nb_xylella {
			
			//spatial attributes
			speed <- rnd(max_speed);
			radius <- rnd(1.0) * max_radius;
			area <- radius * radius * #pi;
			density <- density_map[today.month];
			
			//graphics
			color <- #blue;
			location <- any_location_in(one_of(province));
			geom <- circle(radius);
			
			//infection information
			perc_sane <- 1.0;
			perc_infected <- 0.0;
			
			//initialize current province of swarm
			ask province {
				if inter(self, myself) != nil {
					myself.current_province <- self;
				}
			}
			
			//initialize initial infected area
			geometry gallipoli_area <- circle(radius_start_inf_area,{x_start_inf_area,y_start_inf_area});
			
			//infect part of agents in initial infected zone
			//xyella
			ask xylella {
				geometry inters <- inter(self.geom, gallipoli_area);
				if inters != nil {
					if flip(0.5) {
						self.perc_infected <- rnd(1.0);
						self.perc_sane <- 1 - self.perc_infected;
						self.color <- #red;
					}
				}
			}
			
			//olive
			ask olive {
				geometry inters <- inter(self.geom, gallipoli_area);
				if inters != nil {
					if flip(0.5){
						self.perc_infected <- rnd(1.0);
						self.perc_sane <- 1 - self.perc_infected;
						self.color <- #red;
					}
				}
			}
		}	
	}
	
	//function that evaluates mean percentage of infected in a province
	action get_infected(string pr, string species_){
		
		float total <- 0.0; 
		int n <- 0; 
		
		if species_ = 'xylella' {
			ask xylella {
				if self.current_province.name = pr {
					total <- total + self.perc_infected; 
					n <- n + 1; 
				}
			}
		}
		
		else if species_ = 'olive' {
			ask olive {
				if self.in_province.name = pr {
					total <- total + self.perc_infected; 
					n <- n + 1; 
				}
			}
		}
		
		if n != 0 {
			return total/n; 
		}
		
		else {
			return 0;
		}
	}
	
	//action that returns mean radius of groves in different provinces 
	action get_mean_radius(string pr) {
		float total <- 0.0;
		int n <- 0;
		ask olive {
			if self.in_province.name = pr {
				total <- total + self.radius;
				n <- n + 1;
			} 
		}
		
		if n=0 {
			return 0;
		}
		else {
			return total/n;
		}
	}
	
	
	//function that evaluates beta given the month
	action get_beta(int month) {
		float max_beta <- 0.4;
		int mean_month <- 6; 
		float stdev <- 1.0;
		return max_beta*exp((-1)*(month-mean_month)^2/(2*stdev^2));
	}
	
	//reflex that updates time information
	reflex update_time {
		
		//current date
		today <- today + 1#day;
		
		//if checking period is finished set checked to false for olive and olive_cutted and olive_checked to zero
		if (today.month = 5 or today.month = 9 or today.month = 1) and today.day = 1{
			ask olive {
				self.checked <- false;
			}		
					
			olive_cutted <- 0; 
			olive_checked <- 0;
		}
		
		//update season (0: autumn, 1:winter, 2:spring, 3:summer)
		if today.month in [9,10,11] {
			season <- 0;
		}
		
		else if today.month in [12,1,2] {
			season <- 1;
		}
		
		else if today.month in [3,4,5] {
			season <- 2;
		}
		
		else if today.month in [6,7,8] {
			season <- 3;
		}
		
		//update degree map
		loop v over: infection_graph.vertices {
			if not('xylella' in v) {
				degrees_map[v] <- degree_of(infection_graph, v);
			}
		}
		
		//update cut-off degree if necessary 
		if cut_hubs = true {
			list ordered <- degrees_map sort_by each;
			int cutoff_index <- int(cutoff_perc*length(ordered))-1;
			if cutoff_index != -1 {
				cut_off <- ordered[cutoff_index];
			}
		}
		
		//update beta
		beta <- get_beta(today.month);
		
		//update mean percentae of infected in different provinces
		perc_infected_lecce_x <- get_infected('lecce', 'xylella');
		perc_infected_taranto_x <- get_infected('taranto', 'xylella');
		perc_infected_brindisi_x <- get_infected('brindisi', 'xylella'); 
		perc_infected_lecce_o <- get_infected('lecce', 'olive');
		perc_infected_taranto_o <- get_infected('taranto', 'olive');
		perc_infected_brindisi_o <- get_infected('brindisi', 'olive');
		mean_radius_lecce <- get_mean_radius("lecce");
		mean_radius_brindisi <- get_mean_radius("brindisi");
		mean_radius_taranto <- get_mean_radius("taranto");
	}
	
	
	//reflex that write in log files 
	reflex log when: do_log {

		//log xylella
		save(""+cycle+";"+length(xylella)+";"+mean(xylella collect each.perc_sane)+";"+mean(xylella collect each.perc_infected)+";"
						 +perc_infected_lecce_x+";"+perc_infected_taranto_x+";"+perc_infected_brindisi_x+";"
						 +sum(xylella collect each.geom.area)
		)
			to: path_to_log_xylella type: "text" rewrite: false;
			
		//log olive
		save(""+cycle+";"+mean(olive collect each.perc_sane)+";"+mean(olive collect each.perc_infected)+";"
			             +perc_infected_lecce_o+";"+perc_infected_taranto_o+";"+perc_infected_brindisi_o+";"
			             +mean_radius_lecce+";"+mean_radius_taranto+";"+mean_radius_brindisi+";"
			             +mean(olive collect each.radius)+";"+sum(olive collect each.area)
		)
			to: path_to_log_olive type: "text" rewrite: false;
		
		//olive total: every 100 cycles log the state of each olive
		ask olive {
			if (cycle/100 = float(div(cycle,100))) {
				save(""+cycle+";"+name+";"+perc_sane+";"+perc_infected+";"+radius)
					to: path_to_log_olive_total type: "text" rewrite: false;
			}
		}
	}
	
	
	//reflex that stops the simulation 
	reflex stop_simulation when: cycle = length_simulation {
		do pause ;
	} 
}


/*
 * 
 * Agents
 * 
 */
 
//province species
species province {
	
	rgb color <- #gray;                     //color
	string name;						    //name
	matrix<float> wind_direction_frequency; //matrix of wind direction frequencies
	
	//aspect of agents
	aspect base {
		draw shape color: color;
	}	
}


//xylella species
species xylella skills: [ moving ]{
	
	rgb color;  		                 //color
	float radius;	                	 //swarm radius
	float perc_infected;                 //percentage of infected
	float perc_sane;	                 //percentage of sane
	float density; 		                 //swarm density
	float area <- radius * radius * #pi; //swarm area
	geometry geom;		                 //swarm geometry
	province current_province;           //current province of swarm
	
	//aspect of agents
	aspect base {
		draw circle(radius) color: color;
	}

	
	//reflex that defines swamrs movement
	reflex move {
		
		speed <- rnd(max_speed); //velocità movimento (=lunghezza percorsa)
		
		//with 1% swarms relocate inside their current province
		if flip(relocation_probability) {
			
			//random number needed for wind direction sampling
			float rnd_sample <- rnd(1.0); 
				
			//sample a wind direction (taking into account current province of the swarm)	
			loop i from:0 to:current_province.wind_direction_frequency.columns-1 {
				if current_province.wind_direction_frequency[i,season] > rnd_sample{
					heading <- headings[i-1] - 90.0; //since 0° of GAMA is east 
					break;
				}
				else if i = current_province.wind_direction_frequency.columns-1 {
					heading <- headings[i] - 90.0;   //since 0° of GAMA is east 
					break;
				}
			}
			
			//move in sampled direction with a 22.5° amplitude cone
			do wander bounds: bounds amplitude: movement_amplitude;
			
			//reinitialize swarm geometry (center changed)
			geom <- circle(radius);
			
		}
		else {
			//relocate inside the current province
			location <- any_location_in(one_of(province));
			geom <- circle(radius);
		}
		
		//reinitialize current province of swarm
		ask province {
			if inter(self, myself) != nil {
				myself.current_province <- self;
			}
		}
	}
	
	//reflex that defines infection process (xylella --> olive)
	reflex infects when: perc_infected != 0.0 {
		
		//ask overlapping olive
		ask olive at_distance(radius) {
			
			//find intersection between swarm and olive
			geometry inters <- inter(myself.geom, self.geom);
			
			if inters != nil and myself.geom != nil {
				
				//evaluate probability of infection
				float perc <- beta*myself.density*myself.perc_infected*(inters.area / myself.geom.area);
				
					if flip(perc){
						
						//infect a new percentage of olive
						float new_infected <- self.perc_sane*(inters.area / self.geom.area);
						self.perc_infected <- self.perc_infected + new_infected;
						self.perc_sane <- 1 - self.perc_infected;
						self.color <- #red;
						
						//add a link to the infection network
						infection_graph <- infection_graph add_edge(self.name :: myself.name);
						
						//log
						if do_log {
							save(""+cycle+";"+myself.name+";"+self.name+";"+new_infected)
								to: path_to_log_infection_network type: "text" rewrite: false;
					}		
				}	
			}
		}
	}
	
	//reflex that updates density
	reflex update_density {
		
		if density != density_map[today.month] {
			
			//if month changed ...
			if density_map[today.month] > density {
				
				//density is getting higher since new offsprings are more than dead
				float alpha <- density_map[today.month]/density;
				perc_sane <- (perc_sane+alpha-1)/alpha;
				perc_infected <- 1 - perc_sane;
			}
			
			density <- density_map[today.month];
		}
	}
}


//olive species
species olive {
	
	float area; 						//grove area
	rgb color; 							//color
	float radius <- sqrt(area / #pi); 	//grove radius
	geometry geom <- circle(radius); 	//grove geometry
	float perc_infected; 				//percentage of infected
	float perc_sane; 					//percentage of sane
	bool checked <- false; 				//if true it has already been checked
	float prob_checking; 				//probability of being checked
	province in_province; 				//provice where the grove is located
		
	//agents aspect
	aspect base {
		draw circle(radius) color: color;
	}
	
	//reflex that defines infection process (olive --> xylella)
	reflex infects when: perc_infected != 0.0 {
		
		//ask overlapping xylella
		ask xylella at_distance(radius) {
			
			//find intersection between olive and swarmr
			geometry inters <- inter(myself.geom, self.geom);
			if inters != nil and myself.geom != nil {
				
				//evaluate probability of infetion
				float perc <- beta*myself.perc_infected*(inters.area / myself.geom.area);
				
				if flip(perc){
					
					//infect a new percentage of the swarm
					float new_infected <- self.perc_sane*(inters.area / self.geom.area);
					self.perc_infected <- self.perc_infected + new_infected;
					self.perc_sane <- 1 - self.perc_infected;
					self.color <- #red;
					
					//add a link to the infection network
					infection_graph <- infection_graph add_edge(self.name :: myself.name);
					
					//log
					if do_log {
						save(""+cycle+";"+myself.name+";"+self.name+";"+new_infected)
							to: path_to_log_infection_network type: "text" rewrite: false;
					}
				}	
			}
		}
	}
	
	
	//reflex that defines the check process
	reflex check_date when: today.year >= start_cutting_year{
		if today.month = 4 or today.month = 8 or today.month = 12{
			
			//it is checking period
			if self.checked = false and flip(prob_checking){
				do cut;
				checked <- true;
				olive_checked <- olive_checked + 1;
			}
		}
	}
	
	//action that implements the cut procedure
	action cut {
		
		//cut hubs procedure
		if cut_hubs = true {
			
			if degrees_map[self.name] >= cut_off {
				//in lecce cut a fraction between 0% and 100% of infected
				perc_infected <- 0.0;
				float alpha <- 1 - (perc_sane + perc_infected);
				perc_sane <- 1.0 - perc_infected;
				area <- (1 - alpha)*area; 
				radius <- sqrt(area/#pi);
	   			geom <- circle(radius);
	   			olive_cutted <- olive_cutted + 1;
	   			color <- #gold;
			}
		}
		
		//normal cut	
		else {
			
			//percentage of infected is above 60%
			if perc_infected >= 0.6 and self.in_province != nil{
				
				//distinguish different provinces
				if self.in_province.name = "lecce" {
					
					//in lecce cut a fraction between 0% and 100% of infected
					perc_infected <- perc_infected - perc_infected*rnd(1.0);
					float alpha <- 1 - (perc_sane + perc_infected);
					perc_sane <- 1.0 - perc_infected;
					area <- (1 - alpha)*area; 
					radius <- sqrt(area/#pi);
		   			geom <- circle(radius);
		   			olive_cutted <- olive_cutted + 1;
		   			color <- #black;
				}
				
				else {
					
					//in other provinces cut a fraction between 0% and 50% of infected
					perc_infected <- perc_infected - perc_infected*rnd(0.5);
					float alpha <- 1 - (perc_sane + perc_infected);
					perc_sane <- 1.0 - perc_infected;
					area <- (1 - alpha)*area; 
					radius <- sqrt(area/#pi);
		   			geom <- circle(radius);	
		   			olive_cutted <- olive_cutted + 1;
		   			color <- #black;
				}
			}
		  
		  	//percentage of infected is between 30% and 60%
	    	else if perc_infected >= 0.3 and perc_infected < 0.6 and self.in_province != nil {
	    		
	    		//distinguish different provinces
	    		if self.in_province.name = "lecce" {
	    			
	    			//in lecce cut a fraction between 0% and 50% of infected
					perc_infected <- perc_infected - perc_infected*rnd(0.5);
					float alpha <- 1 - (perc_sane + perc_infected);
					perc_sane <- 1.0 - perc_infected;
					area <- (1 - alpha)*area; 
					radius <- sqrt(area/#pi);
		   			geom <- circle(radius);
		   			olive_cutted <- olive_cutted + 1;
		   			color <- #yellow;
				}
				
				else {
					
					//in other provinces cut a fraction between 0% and 25% of infected
					perc_infected <- perc_infected - perc_infected*rnd(0.25);
					float alpha <- 1 - (perc_sane + perc_infected);
					perc_sane <- 1.0 - perc_infected;
					area <- (1 - alpha)*area; 
					radius <- sqrt(area/#pi);
		   			geom <- circle(radius);	
		   			olive_cutted <- olive_cutted + 1;
		   			color <- #yellow;
				}	   	
	    	}    
		}
	}	 
}


/*
 * 
 * Experiment
 * 
 */
	

experiment xylella_exp type: gui {
	
	output {
		
		display region_display type:opengl {
			species province aspect: base;
			species xylella aspect: base;
			species olive aspect: base;
		}
		
	display infection_display refresh: every(1) { 
			chart "% of infected - total" type: series {
				data "Mean of % infected - Xylella" value: mean (xylella collect each.perc_infected) style: line color: #blue;
				data "Mean of % infected - Olive" value: mean (olive collect each.perc_infected) style: line color: #yellow;
			}
		}	
		
	display infection_display_xylella refresh: every(1) { 
			chart "% of infected - xylella" type: series {
				data "lecce" value: perc_infected_lecce_x style: line color: #blue;
				data "taranto" value: perc_infected_taranto_x  style: line color: #red;
				data "brindisi" value: perc_infected_brindisi_x  style: line color: #yellow;
			}
		}
		
	display infection_display_olive refresh: every(1) { 
			chart "% of infected - olive" type: series {
				data "lecce" value: perc_infected_lecce_o style: line color: #blue;
				data "taranto" value: perc_infected_taranto_o  style: line color: #red;
				data "brindisi" value: perc_infected_brindisi_o  style: line color: #yellow;
			}
		}
		
	display cut_display refresh: every(1) { 
			chart "olives checked and cutted" type: series {
				data "Olives cutted" value: olive_cutted style: line color: #red;
				data "Olives checked" value: olive_checked style: line color: #blue;
				data "Total olives" value: length(olive) style: line color: #black;
			}
		}
		
	display radius_display refresh: every(1) { 
			chart "Mean Radius - Province" type: series {
				data "Mean Radius Lecce" value: mean_radius_lecce style: line color: #blue;
				data "Mean Radius Brindisi" value: mean_radius_brindisi style: line color: #red;
				data "Mean Radius Taranto" value: mean_radius_taranto style: line color: #yellow;
			}
		}
	}
}