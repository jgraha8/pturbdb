#include <mpi.h>

class Clock {
	double start_;
	double stop_;
public:
	Clock(): start_(0.0L), stop_(0.0L){};

	void start(){ this->start_ = MPI_Wtime(); }
	void stop(){ this->stop_ = MPI_Wtime(); }
	void resume(){ double _time = this->time(); this->start_ = MPI_Wtime() - _time; }
	double time(){ return this->stop_ - this->start_; }	
};
