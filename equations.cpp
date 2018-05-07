#include <iostream>
#include <iomanip>
#include <vector>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <assert.h>
#include <getopt.h>

using namespace std;

#define RIGID    100
#define MOLDABLE 200
#define GRID     300
#define NBMODELS 3

static int models[NBMODELS] = { 100, 200, 300 };

struct Params
{
    int    model;
	int    P;
	double MTBF;
	double waittime;
	double C;
};

static string modelname(int model)
{
    if( MOLDABLE == model )
        return string("MOLDABLE");
    if( RIGID == model )
        return string("RIGID");
    if( GRID == model )
        return string("GRID");
    return string("");
}

static void printParams(Params& p)
{
    if( p.P != 0 )
        cout << "#nb-procs = " << p.P << std::endl;
    if( p.MTBF != 0 )
        cout << "#MTBF = " << p.MTBF << std::endl;
    if( p.waittime != 0 )
        cout << "#wait-time = " << p.waittime << std::endl;
    if( p.C != 0 ) 
        cout << "#ckpt-time = " << p.C << std::endl;
    if( p.model != 0 )
        cout << "#Model = " << modelname(p.model) << std::endl;
}

static double mu(int n, const Params& p)
{
	return p.MTBF/n;
}

static double C(int i,const Params& p)
{
	return p.C;
}

static double rigidTime(int F, const Params& p)
{
	double result = 0;
	double result2 = 0;
	for (int i=p.P; i>=p.P-F; i--)
	{
		result += mu(i,p);
		result2 += (double)(p.P-F)/(double)(i)*(C(p.P-F,p)+sqrt(2*C(p.P-F,p)*mu(p.P-F,p))/2);
	}
	//cout << result << " " << result2 << "\n";
	result += p.waittime + result2;
	return result;
}

static double rigidWork(int F, const Params& p)
{
	double result = 0;
	for (int i=p.P; i>=p.P-F; i--)
	{
		result += mu(i,p)*(1-C(p.P-F,p)/sqrt(2*C(p.P-F,p)*mu(p.P-F,p)));
	}
	return result*(p.P-F);
}

static double rigidYield(int Q, const Params& p)
{
	double w = rigidWork(Q,p);
	double t = rigidTime(Q,p);
	return w/(p.P*t);
}

/* Cheap-ass differentiation for simple functions whose domain is integer
 * returns (a bad approx of) f'_p(x) */
static double differ(double (*f)(int i, const Params &p), int x, const Params &p)
{
    return f(x + 1, p) - f(x, p);
}

/* Assume f_p is a convex function in [0, p.P] that returns a value in [0, 1]
 * This function finds the index and value for which f_p is max in [0, p.P]
 *   It uses the fact that 
 *   if f'_p(a)*f'_p((a+b)/2) < 0.0, then the max is in [a, (a+b)/2]
 *   else it is in [(a+b/2), b]
 */
static pair<int, double> findMax(const Params& p, double (*yield_fn)(int i, const Params &p))
{
    double m = 0.0;
    int mp = -1;
    int lower = 0, upper = p.P - 1; // P-1 to avoid dividing by 0

    /* Manage the case where there is no trade-off: yield_fn is strictly
     * decreasing */
    if( differ(yield_fn, 0, p) < 0.0 ) {
        if( yield_fn(0, p) < 0.0 )
            return make_pair(-1, 0.0);
        return make_pair(0, yield_fn(0, p));
    }
    
    while( upper - lower > 32 ) {
        int mid = (lower + upper) / 2;
        double s = differ(yield_fn, lower, p) * differ(yield_fn, mid, p);
        if( s > 0.0 ) {
            lower = mid;
        } else {
            upper = mid;
        }
    }

    for (int i = lower; i < upper; i++) {
        double v = yield_fn(i, p);
        assert( v >= 0.0 );
        assert( v <= 1.0 );
        if( v > m ) {
            m = v;
            mp = i;
        } else if( v <= m*0.95 ) { // Stop computation when reaching 5% less than actual max
            return make_pair(mp, m);
        }
    }
    return make_pair(mp, m);
}

static double findMaxY(double target, const Params& orig, double (*yield_fn)(int i, const Params &p))
{
    Params p = orig;

    double m = 0.0;
    int mp = -1;
    double lower = 0, upper = 1000*365.65*24*3600;
    p.waittime = lower;
    double yieldLow = yield_fn(0,p);
    p.waittime = upper;
    double yieldUp = yield_fn(0,p);

    while( (yieldLow - yieldUp > 1e-5 || upper - lower > 60) && (upper > lower+1) ) {
	//cout << lower << " " << upper << ", " << yieldLow << " " << yieldUp << "\n";
        double mid = (lower + upper) / 2;
	p.waittime = mid;
	double val = yield_fn(0,p);
        if( val > target ) {
            lower = mid;
	    yieldLow = val;
        } else {
            upper = mid;
	    yieldUp = val;
        }
    }

    return lower; //ensure being better than target
}

static double findDoubleMaxY(double target, const Params& orig, double (*yield_fn)(int i, const Params &p))
{
    Params p = orig;

    double m = 0.0;
    int mp = -1;
    double lower = 0, upper = 1000*365.65*24*3600;
    p.waittime = lower;
    double yieldLow = yield_fn(findMax(p,yield_fn).first,p);
    p.waittime = upper;
    double yieldUp = yield_fn(findMax(p,yield_fn).first,p);

    while( (yieldLow - yieldUp > 1e-5 || upper - lower > 60) && (upper > lower+1) ) {
	//cout << lower << " " << upper << ", " << yieldLow << " " << yieldUp << "\n";
        double mid = (lower + upper) / 2;
	p.waittime = mid;
	double val = yield_fn(findMax(p,yield_fn).first,p);
        if( val > target ) {
            lower = mid;
	    yieldLow = val;
        } else {
            upper = mid;
	    yieldUp = val;
        }
    }

    return lower; //ensure being better than target
}

static void plot(const Params& p, double (*yield_fn)(int i, const Params &p))
{
    for (int i = 0; i < 0.1 * p.P; i++) {
        double v = yield_fn(i, p);
        cout << i << "\t" << v << "\t" << differ(yield_fn, i, p) << std::endl;
    }
}

static double moldableTime(int F, const Params& p)
{
	double result = 0;
	double result2 = 0;
	for (int i=p.P; i>p.P-F; i--)
	{
		result += mu(i,p);
		result2 += C(i,p) + (double)i/(double)(i-1)*sqrt(2*C(i,p)*mu(i,p))/2;
	}
	//cout << result << " " << result2 << "\n";
	//last iteration is different:
	result += mu(p.P-F,p) + C(p.P,p) + (double)(p.P-F)/(double)p.P*sqrt(2*C(p.P-F,p)*mu(p.P-F,p))/2 + p.waittime + result2;
	return result;
}

static double moldableWork(int F, const Params& p)
{
	double result = 0;
	for (int i=p.P; i>=p.P-F; i--)
	{
		result += mu(i,p)*(1-C(i,p)/sqrt(2*C(i,p)*mu(i,p)))*i;
	}
	return result;
}

static double moldableYield(int L, const Params& p)
{
	double w = moldableWork(L,p);
	double t = moldableTime(L,p);
	return w/(p.P*t);
}

static double C(int P, int Q, const Params& p)
{
	return p.C;
}

static double gridWork(int F, const Params& p)
{
	int p1 = sqrt(p.P);
	int p2 = sqrt(p.P);
	int spares = 0;
	
	double result = 0;
	for (int i=0; i<=F; i++)
	{
		result += mu(p.P-i,p)*(1-C(p1,p2,p)/sqrt(2*C(p1,p2,p)*mu(p1*p2,p)))*(p1*p2);
		if (spares == 0)
		{
			if (p1 < p2)
			{
				spares = p1-1;
				p2--;
			} else {
				spares = p2-1;
				p1--;
			}
		} else {
			spares--;
		}
	}

	return result;
}

static double gridTime(int F, const Params& p)
{
	int p1 = sqrt(p.P);
	int p2 = sqrt(p.P);
	int spares = 0;
	
	double result = 0;
	for (int i=0; i<F; i++)
	{
		result += mu(p.P-i,p);
		if (spares == 0)
		{
			if (p1 < p2)
			{
				result += C(p1,p2-1,p) + (double)p2/(double)(p2-1)*sqrt(2*C(p1,p2,p)*mu(p1*p2,p))/2;
				spares = p1-1;
				p2--;
                                //cout << "u_2i+2 = " << i+1 << "\n";
			} else {
				result += C(p1-1,p2,p) + (double)p1/(double)(p1-1)*sqrt(2*C(p1,p2,p)*mu(p1*p2,p))/2;
				spares = p2-1;
				p1--;
                                //cout << "u_2i+1 = " << i+1 << "\n";
			}
		} else {
			result += (double)(p1*p2)/(double)(p.P-i)*(C(p1,p2,p)+sqrt(2*C(p1,p2,p)*mu(p1*p2,p))/2);
			spares--;
		}
	}
	//need to check
	result += mu(p.P-F,p) + p.waittime + C(sqrt(p.P),sqrt(p.P),p) + (double)(p1*p2)/(double)(p.P)*(sqrt(2*C(p1,p2,p)*mu(p1*p2,p))/2);

	return result;
}

static double gridYield(int L, const Params& p)
{
	double w = gridWork(L,p);
	double t = gridTime(L,p);
	return w/(p.P*t);
}

static void waittime_expe(Params& orig)
{
    Params p = orig;
	p.waittime = 0;
    p.model = 0;
	printParams(p);
    cout << "Wait Time";
    for(unsigned int j = 0; j < NBMODELS; j++) {
        cout << "," << modelname(models[j]) << "_F," << modelname(models[j]) << "_Yield," << modelname(models[j]) << "_Work," << modelname(models[j]) << "_Time";
    }
    cout << ",NO_Yield,NO_Work,NO_Time" << std::endl;
	for (int i=1; i<=120; i++)
	{
        cout << p.waittime;
        for(unsigned int j = 0; j < NBMODELS; j++ ) {
            p.model = models[j];
            
            pair<int,double> res;
            if (p.model == RIGID)
                res = findMax(p, rigidYield);
            else if (p.model == MOLDABLE)
                res = findMax(p, moldableYield);
            else if (p.model == GRID)
                res = findMax(p, gridYield);
            cout << "," << res.first << "," << res.second;

            if (p.model == RIGID)
                cout << "," << rigidWork(res.first, p) << "," << rigidTime(res.first, p);
            else if (p.model == MOLDABLE)
                cout << "," << moldableWork(res.first, p) << "," << moldableTime(res.first, p);
            else if (p.model == GRID)
                cout << "," << gridWork(res.first, p) << "," << gridTime(res.first, p);
        }
        assert( fabs(rigidWork(0, p) - moldableWork(0, p)) < 1e-10 );
        assert( fabs(rigidYield(0, p) - moldableYield(0, p)) < 1e-10);
        cout << "," << moldableYield(0, p) << "," << moldableWork(0, p) << "," << moldableTime(0, p) << std::endl;

		p.waittime += 600;
	}
}

static void proc_expe(Params& orig)
{
    Params p = orig;
	p.P = 0;
    p.model = 0;
	printParams(p);
    cout << "Number of procs";
    for(unsigned int j = 0; j < NBMODELS; j++) {
        cout << "," << modelname(models[j]) << "_F," << modelname(models[j]) << "_Yield," << modelname(models[j]) << "_Work," << modelname(models[j]) << "_Time";
    }
    cout << ",NO_Yield,NO_Work,NO_Time" << std::endl;
	for (int i=1; i<=100; i++)
	{
        p.P = sqrt(p.P+5000);
        p.P = p.P*p.P;
        cout << p.P;
        for(unsigned int j = 0; j < NBMODELS; j++ ) {
            p.model = models[j];
            
            pair<int,double> res;
            if (p.model == RIGID)
                res = findMax(p, rigidYield);
            else if (p.model == MOLDABLE)
                res = findMax(p, moldableYield);
            else if (p.model == GRID)
                res = findMax(p, gridYield);
            cout << "," << res.first << "," << res.second;

            if (p.model == RIGID)
                cout << "," << rigidWork(res.first, p) << "," << rigidTime(res.first, p);
            else if (p.model == MOLDABLE)
                cout << "," << moldableWork(res.first, p) << "," << moldableTime(res.first, p);
            else if (p.model == GRID)
                cout << "," << gridWork(res.first, p) << "," << gridTime(res.first, p);
        }
        assert( fabs(rigidWork(0, p) - moldableWork(0, p)) < 1e-10 );
        assert( fabs(rigidYield(0, p) - moldableYield(0, p)) < 1e-10);
        cout << "," << moldableYield(0, p) << "," << moldableWork(0, p) << "," << moldableTime(0, p) << std::endl;

	}
}

static void ckpt_expe(Params& orig)
{
    Params p = orig;
	p.C = 0;
    p.model = 0;
	printParams(p);
    cout << "Checkpoint Time";
    for(unsigned int j = 0; j < NBMODELS; j++) {
        cout << "," << modelname(models[j]) << "_F," << modelname(models[j]) << "_Yield," << modelname(models[j]) << "_Work," << modelname(models[j]) << "_Time";
    }
    cout << ",NO_Yield,NO_Work,NO_Time" << std::endl;
	for (int i=1; i<=120; i++)
	{
		p.C += 30;
        cout << p.C;
        for(unsigned int j = 0; j < NBMODELS; j++ ) {
            p.model = models[j];
            
            pair<int,double> res;
            if (p.model == RIGID)
                res = findMax(p, rigidYield);
            else if (p.model == MOLDABLE)
                res = findMax(p, moldableYield);
            else if (p.model == GRID)
                res = findMax(p, gridYield);
            cout << "," << res.first << "," << res.second;

            if (p.model == RIGID)
                cout << "," << rigidWork(res.first, p) << "," << rigidTime(res.first, p);
            else if (p.model == MOLDABLE)
                cout << "," << moldableWork(res.first, p) << "," << moldableTime(res.first, p);
            else if (p.model == GRID)
                cout << "," << gridWork(res.first, p) << "," << gridTime(res.first, p);
        }
        assert( fabs(rigidWork(0, p) - moldableWork(0, p)) < 1e-10 );
        assert( fabs(rigidYield(0, p) - moldableYield(0, p)) < 1e-10);
        cout << "," << moldableYield(0, p) << "," << moldableWork(0, p) << "," << moldableTime(0, p) << std::endl;

	}
}

static void mtbf_expe(Params& orig)
{
    Params p = orig;
	p.MTBF = 0;
    p.model = 0;
	printParams(p);
    cout << "MTBF";
    for(unsigned int j = 0; j < NBMODELS; j++) {
        cout << "," << modelname(models[j]) << "_F," << modelname(models[j]) << "_Yield," << modelname(models[j]) << "_Work," << modelname(models[j]) << "_Time";
    }
    cout << ",NO_Yield,NO_Work,NO_Time" << std::endl;
	for (int i=1; i<=120; i++)
	{
		p.MTBF += 3600*24*365.25;
        cout << p.MTBF;
        for(unsigned int j = 0; j < NBMODELS; j++ ) {
            p.model = models[j];
            
            pair<int,double> res;
            if (p.model == RIGID)
                res = findMax(p, rigidYield);
            else if (p.model == MOLDABLE)
                res = findMax(p, moldableYield);
            else if (p.model == GRID)
                res = findMax(p, gridYield);
            cout << "," << res.first << "," << res.second;

            if (p.model == RIGID)
                cout << "," << rigidWork(res.first, p) << "," << rigidTime(res.first, p);
            else if (p.model == MOLDABLE)
                cout << "," << moldableWork(res.first, p) << "," << moldableTime(res.first, p);
            else if (p.model == GRID)
                cout << "," << gridWork(res.first, p) << "," << gridTime(res.first, p);
        }
        assert( fabs(rigidWork(0, p) - moldableWork(0, p)) < 1e-10 );
        assert( fabs(rigidYield(0, p) - moldableYield(0, p)) < 1e-10);
        cout << "," << moldableYield(0, p) << "," << moldableWork(0, p) << "," << moldableTime(0, p) << std::endl;

	}
}

static void downtime_opt(Params& orig, bool maxmax)
{
	Params p = orig;
	double target = 0.5;
	for (int i=1; i<=200; i++)
	{
		target += 0.002;
		double dOpt = 0;
		cout << target << " ";
		if (maxmax)
		    dOpt = findDoubleMaxY(target,p,rigidYield);
		else
		    dOpt = findMaxY(target,p,rigidYield);
		p.waittime = dOpt;
		if (maxmax)
		{
			if (fabs(findMax(p,rigidYield).second - target) < 1e-4)
				cout << dOpt << " ";
			else
				cout << "0 ";
		} else {
			if (fabs(rigidYield(0,p) - target) < 1e-4)
				cout << dOpt << " ";
			else
				cout << "0 ";
		}
		if (maxmax)
		    dOpt = findDoubleMaxY(target,p,moldableYield);
		else
		    dOpt = findMaxY(target,p,moldableYield);
		p.waittime = dOpt;
		if (maxmax)
		{
			if (fabs(findMax(p,moldableYield).second - target) < 1e-4)
				cout << dOpt << " ";
			else
				cout << "0 ";
		} else {
			if (fabs(moldableYield(0,p) - target) < 1e-4)
				cout << dOpt << " ";
			else
				cout << "0 ";
		}
		if (maxmax)
		    dOpt = findDoubleMaxY(target,p,gridYield);
		else
		    dOpt = findMaxY(target,p,gridYield);
		p.waittime = dOpt;
		if (maxmax)
		{
			if (fabs(findMax(p,gridYield).second - target) < 1e-4)
				cout << dOpt << " ";
			else
				cout << "0 ";
		} else {
			if (fabs(gridYield(0,p) - target) < 1e-4)
				cout << dOpt << " ";
			else
				cout << "0 ";
		}
		cout << "\n";	
	}
}

static void spectrum(void)
{
    Params p;

    cout << "#P,MTBF,WaitTime,CkptTime,RF,RYield,RWork,MF,MYield,MWork" << std::endl;
    for(int P = 1024; P <= 512*1024; P = P*2) {
        p.P = P;
        for(double mtbf = 5.0*365*24*3600;
            mtbf < 100.0*365*24*3600;
            mtbf += 10.0*365*24*3600) {
            p.MTBF = mtbf;
            for(double d = 3600; d <= 240.0*3600; d = d * 2) {
                p.waittime = d;
                for(double C = 60.0; C <= 1920.0; C = C * 2.0) {
                    p.C = C;
                    pair<int, double>res1 = findMax(p, rigidYield);
                    pair<int, double>res2 = findMax(p, moldableYield);
                    if( res1.first != -1 && res2.first != -1 ) {
                        cout << P << "," << mtbf << "," << d << "," << C << ","
                             << res1.first << "," << res1.second << ","
                             << rigidWork(res1.first, p) << ","
                             << res2.first << "," << res2.second << ","
                             << moldableWork(res2.first, p) << std::endl;;
                    }
                }
            }
        }
    }
}

static int plot(const Params &p)
{
    if( MOLDABLE == p.model )
        plot(p, moldableYield);
    else if( RIGID == p.model ) 
        plot(p, rigidYield);
    else
        plot(p, gridYield);
    return 0;
}

static bool startsWith(const char *pre, const char *str)
{
    size_t lenpre = strlen(pre),
           lenstr = strlen(str);
    return lenstr < lenpre ? false : strncmp(pre, str, lenpre) == 0;
}

int main(int argc, char **argv)
{
    static struct option longopts[] = {
        { "Model",      required_argument,            NULL,           'M' },
        { "pLot",       no_argument,                  NULL,           'L' },
        { "Proc-expe",  no_argument,                  NULL,           'P' },
        { "Ckpt-expe",  no_argument,                  NULL,           'C' },
        { "Waittime-expe",  no_argument,                  NULL,           'W' },
        { "MTBF-expe",  no_argument,                  NULL,           'Y' },
        { "Spectrum",   no_argument,                  NULL,           'S' },
        { "Test",       no_argument,                  NULL,           'T' },
        { "nb-procs",      required_argument,            NULL,           'n' },
        { "mtbf",          required_argument,            NULL,           'm' },
        { "wait-time",     required_argument,            NULL,           'w' },
        { "ckpt-time",     required_argument,            NULL,           'c' },
        { NULL,               0,                      NULL,           0 }
    };
	Params p;
    int ch;
    
	p.P = 150*150;
	p.MTBF = 20*365.25*24*3600;
	p.waittime = 3600*12;
	p.C = 120;
    p.model = RIGID;

    while ((ch = getopt_long(argc, argv, "n:m:w:c:CSTPLWYM:", longopts, NULL)) != -1) {
        switch (ch) {
        case 'M':
            if( startsWith(optarg, "moldable") ) {
                p.model = MOLDABLE;
            } else if ( startsWith(optarg, "rigid") ) {
                p.model = RIGID;
            } else if ( startsWith(optarg, "grid") ) {
                p.model = GRID;
            } else {
                cerr << "Unknown model type '" << optarg << "' -- ignored" << std::endl;
            }
            break;
        case 'P':
            proc_expe(p);
            break;
        case 'C':
            ckpt_expe(p);
            break;
        case 'W':
            waittime_expe(p);
            break;
        case 'Y':
            mtbf_expe(p);
            break;
        case 'L':
            plot(p);
            break;
        case 'S':
            spectrum();
            break;
        case 'T':
            {
                int op = p.P;

		/*double d = findMaxY(0.5,p,rigidYield);
		p.waittime = d;
		cout << d << " " << rigidYield(findMax(p,rigidYield).first,p) << " " << findMax(p,rigidYield).second << "\n";*/

		downtime_opt(p,true);
	
		/*pair<int,double> res = findMax(p,rigidYield);
		cout << rigidYield(res.first,p) << " " << rigidWork(res.first,p) << " " << rigidTime(res.first,p) << "\n";
		res = findMax(p,moldableYield);
		cout << moldableYield(res.first,p) << " " << moldableWork(res.first,p) << " " << moldableTime(res.first,p) << "\n";
		res = findMax(p,gridYield);
		cout << gridYield(res.first,p) << " " << gridWork(res.first,p) << " " << gridTime(res.first,p) << "\n";*/
                p.P = op;
            }
            break;
        case 'n':
            p.P = atoi(optarg);
            break;
        case 'm':
            p.MTBF = strtod(optarg, NULL);
            break;
        case 'w':
            p.waittime = strtod(optarg, NULL);
            break;
        case 'c':
            p.C = strtod(optarg, NULL);
            break;
        }
    }

	return 0;
}
