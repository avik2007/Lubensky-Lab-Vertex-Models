// include guard
#ifndef __BESSEL_H_INCLUDED__
#define __BESSEL_H_INCLUDED__



class Bessel
{
public:
	double bessj(int order, double input);
	double bessy(int order, double input);
	double bessi(int order, double input);
	double bessk(int order , double input);
};
#endif