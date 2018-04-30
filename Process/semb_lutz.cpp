  omp_set_num_threads(_sct);
  #pragma omp parallel shared 
           (traveltime, trace, nstep, nostat, nsamp, ntimes, sembmax,latv,lonv,sembmaxX,sembmaxY) 
           private(i,j,k,l,x,y) reduction(+:semb,nomin,denom,sum,relstart_samples)   
  {
	cout <<"Thread Number "<<omp_get_thread_num()<<",from nthreads "<<omp_get_num_threads()<<endl; 
	#pragma omp for schedule (static,1) 
	
	for (i=0;i<ntimes;i++) 
      	{
	//+++++++++display data in files ++++ for each timestep one file +++++folders name same as event time++++
		timestep = boost::lexical_cast <double> (_origin->time().value().toString("%s"))+i*nstep; 
		cout << "TIMESTEP " << timestep <<endl;

		sembpath = path+"/"+boost::lexical_cast<std::string>(i)+".ASC";
		cout << "SEMB " << sembpath << endl;

		sembout.open(sembpath.c_str(),ios::out);
	
	//	if (!sembout) {	assert (false);}
	//++++++++++++++++++++++++write calculated data into file with header+++++++++++++++++++++++++++++++++++++++++++

		sembout << "# "<<timestep <<" , "<<recordstarttime <<endl<< "# CalcWindowStart "
                      <<_calcTimeWindow.startTime().toString("%F %T")<<endl<< "# CalcWindowEnd   "
                      <<_calcTimeWindow.endTime().toString("%F %T")<<endl<< "# step: "<<_step<<"s| ntimes: "
                      <<ntimes << "| winlen: "<<_winlen<<"s"<<endl;

		sembout << "# ";

//		for (unsigned int u = 0; u < sembstreams.size(); u++)  	{sembout <<sembstreams[u]<< " ";}
		sembout << endl<< "# southwestlat: "<<Latul<<" dlat: "<<_gridspacing<<
			              " nlat: "<<_dimX<<endl<< "# southwestlon: "<<Lonul<<" dlon: "<<_gridspacing<<
				          " nlat: "<<_dimY<<endl<< "#ddepth: "<<"0"<<" ndepth: "<<"1"<<endl;

		// loop over grid points
		sembmax  = 0; sembmaxX = 0;	sembmaxY = 0;

		for (j=0; j < _dimX * _dimY; j++){
			semb = 0;	nomin = denom = 0;
				
			for (l=0;l<nsamp;l++){
		 	    sum=0;	

			    for (k=0;k<nostat;k++) {
                    relstart = (int)((traveltime[k][j]-mint) * _new_frequence + 0.5) + i*nstep;
			        sum   += trace[k][relstart+l];
		            denom += trace[k][relstart+l] * trace[k][relstart+l];
		         			 
			    } //for nostat

	 		    nomin += sum*sum;
			} //for nsamp
	 		
	 		x    = latv[j];
			y    = lonv[j];
			semb = nomin /((double)nostat*denom);

		 	sembout <<x <<" "<< y<<" "<< semb << endl;
				
			if ( semb > sembmax ){
				sembmax  = semb;
                        // search for maximum and position of maximum on semblance grid for given time step	   
				sembmaxX = latv[j];
				sembmaxY = lonv[j];				
			}//if
		} //for _dimX*_dimY

		sembout << "# maximum semblance: "<<sembmax<<" at latitude: "
                      << sembmaxX<<" at longitude: "<< sembmaxY<<endl;
		sembout.close();
		
	} // for ntimes
	
  } // pragma end	
