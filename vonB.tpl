//><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>//
//Programer:													 //
//Date:															 //
//Purpose:Fit a simple vonBertalanffy model to length and age 	 //
//Notes: LOCAL_CALC caused issues with ADMB10_IDE				 //
//could not reproduce this behaviour							 //
//																 //
//><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>//

DATA_SECTION
	int sim;
	int rseed;
	LOCAL_CALCS //set up a command line option to generate data with randon seed option
		sim=0;
		rseed=0;
		int on,opt;
		if((on=option_match(ad_comm::argc,ad_comm::argv,"-sim",opt))>-1)
		{
			sim=1;
			rseed=atoi(ad_comm::argv[on+1]);
		}
	END_CALCS
	init_int nages; //read in the max age from the data file
	init_int nobs; //number of observations
	init_vector a(1,nobs); // vector for ages
	init_vector l(1,nobs); //vector for lengths
	init_int eofdat;//end of tile flag
	LOCAL_CALCS
		if(eofdat!=999)
		{
			cout<<"Error reading data file.\n Fix it."<<endl;
			ad_exit(1);
		}
	END_CALCS
	
	!!ad_comm::change_datafile_name("vonB.ctl");
	init_number tlinf;//these are initial parameter values and needed to simulate data
	init_number tvbk;
	init_number tto;
	init_number tcv;
	init_number tM;//These are needed to simulate the data
	init_number tselh;
	init_number tselsd;
	init_int eofctl;//end of tile flag
	int iter;
	!!iter=0;
	LOCAL_CALCS
		if(eofctl!=999)
		{
			cout<<"Error reading control file.\n Fix it."<<endl;
			ad_exit(1);
		}
	END_CALCS


PARAMETER_SECTION
	init_number log_linf;//we are going to search on the log transformed so that the value in regular space is >0
	init_number log_vbk;
	init_number log_to;
	init_number log_cv;
	likeprof_number linf; //by defining a likeprof_number we can use the -lprof command line option 
	likeprof_number vbk;//to generate a likelihood profile we also create a parameter value to transform the log_value
	number to;
	number cv;//
	//!!log_linf=log(ilinf);// set initial values can also use the .pin file for this. Need to add this later
	//!!log_vbk=log(ivbk);
	//!!log_cv=log(icv);
	objective_function_value nll;//must define the objective function
	vector lhat(1,nobs);// we are not searching on this vector but do need to keep track of the derivatives of these values

PRELIMINARY_CALCS_SECTION
	if(sim)
	{
		run_data_simulation();
	}

PROCEDURE_SECTION
	initialization();
	calculations();
	objectivefunction();
	if(mceval_phase()) mcmc_output();

FUNCTION initialization
	linf=mfexp(log_linf);//transfor the log parameter that are being searched on to regular space.
	vbk=mfexp(log_vbk);
	to=mfexp(log_to);
	cv=mfexp(log_cv);

FUNCTION calculations
	//calcuate mean length at age 
	lhat=linf*(1.-mfexp(-vbk*(a+to)));

FUNCTION objectivefunction 
	//calculate objective function
	dvar_vector epsilon=log(elem_div(l,lhat));//calculate deviations from the mean
	dvar_vector sigma=cv*lhat;//create a vector of standard deviations for each observation depending on the estiamted mean length at age
	nll=dnorm(epsilon,sigma);//use the statslib dnorm function to caalculate the total negative log likelihood.

FUNCTION mcmc_output

	dvector iage(1,nages);
	dvector ilhat(1,nages);
	iage.fill_seqadd(1.,1.);
	ilhat=value(linf)*(1.-mfexp(-value(vbk)*(iage+value(to))));
	if(iter==0)
	{
		ofstream ofs("lpars.mcmc");
		ofs<<"linf\t vbk\t to\t cv\t"<<endl;
		ofstream ofs1("lbar.mcmc");
	}
	iter++;
	ofstream ofs("lpars.mcmc",ios::app);
	ofs<<linf<<"\t"<<vbk<<"\t"<<to<<"\t"<<cv<<endl;
	ofstream ofs1("lbar.mcmc",ios::app);
	ofs1<<ilhat<<endl;
	

FUNCTION run_data_simulation
	random_number_generator rng(rseed);
	dvector iage(1.,nages);
	iage.fill_seqadd(1.,1.);
	dvector lx=pow(mfexp(-tM),iage-1.);			//survivorship
	dvector pa=lx/sum(lx);						//proportion-at-age
	dvector sa=plogis(iage,tselh,tselsd);				//logistic selectivity
	dvector p=elem_prod(pa,sa);				//probability of capturing a fish of age a
	a.fill_multinomial(rng,p);//draw ages from random multinomial distribution based on the probability of capture
	dvector err(1.,nobs); //create the err vector
	err.fill_randn(rng); //fill err with standard normal random deviates
	l=tlinf*(1.-mfexp(-tvbk*(a+tto))); //calculate expected mean lengths based on specificed vonB paramters
	err=elem_prod(err,tcv*l); // change standard normal error to lenght dependent standard deviations
	l+=err; //asdd error to mean length to create observations

REPORT_SECTION
	report<<"linf\n"<<linf<<endl;// write to the report file
	report<<"vbk\n"<<vbk<<endl;
	report<<"to\n"<<to<<endl;
	report<<"cv\n"<<cv<<endl;
	report<<"ages\n"<<a<<endl;
	report<<"lengths\n"<<l<<endl;
	report<<"meanlength\n"<<lhat<<endl;
GLOBALS_SECTION
	/**
	\def REPORT(object)
	Prints name and value of \a object on ADMB report %ofstream file.
	*/
	#undef REPORT
	#define REPORT(object) report << #object "\n" << object << endl;

	#include <admodel.h>
	#include <time.h>
	#include <contrib.h>//IF you have ADMB-11
	//#include<stats.cxx>//If you have ADMB-10 and make sure stats.cxx is in your working directory
	time_t start,finish;
	long hour,minute,second;
	double elapsed_time;

FINAL_SECTION
	time(&finish);
	elapsed_time=difftime(finish,start);
	hour=long(elapsed_time)/3600;
	minute=long(elapsed_time)%3600/60;
	second=(long(elapsed_time)%3600)%60;
	cout<<"*******************************************"<<endl;
	cout<<"--Start time: "<<ctime(&start)<<endl;
	cout<<"--Finish time: "<<ctime(&finish)<<endl;
	cout<<"--Runtime: ";
	cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
	cout<<"*******************************************"<<endl;


